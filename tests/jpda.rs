//! JPDA update and tracker checks.

use fusetrack::ekf::Ekf;
use fusetrack::jpda::{pda_update, JpdaConfig};
use fusetrack::models::{Radar, RadarMeas, State};
use fusetrack::sim::{simulate_scenario, ScenarioConfig};
use fusetrack::tracker::{Track, TrackFilter, TrackStatus};
use fusetrack::{GnnTracker, ScoreConfig, TrackerConfig};
use nalgebra::{Matrix4, Vector2};
use rand::{rngs::StdRng, SeedableRng};

fn radar() -> Radar {
    Radar {
        pos: Vector2::new(0.0, 0.0),
        sigma_r: 15.0,
        sigma_b: 0.003,
    }
}

fn ekf() -> Ekf {
    let x = State::new(4000.0, 1500.0, 20.0, -10.0);
    let p = Matrix4::from_diagonal(&nalgebra::Vector4::new(400.0, 400.0, 25.0, 25.0));
    Ekf::new(x, p, 1.0)
}

fn meas(dr: f64, db: f64, e: &Ekf, r: &Radar) -> RadarMeas {
    let z = r.h(&e.x);
    RadarMeas {
        range: z.range + dr,
        bearing: z.bearing + db,
    }
}

#[test]
fn a_certain_single_measurement_gives_exactly_the_ekf_update() {
    let (r, a) = (radar(), ekf());
    let z = meas(8.0, 0.002, &a, &r);
    let inn = a.innovation(&z, &r).unwrap();
    let mut plain = a.clone();
    plain.update_with(&inn, &r);
    let mut pda = a.clone();
    pda_update(&mut pda, &r, &[(1.0, inn)], 0.0);
    assert!((plain.x - pda.x).norm() < 1e-9);
    assert!((plain.p - pda.p).norm() < 1e-6);
}

#[test]
fn no_detection_leaves_the_estimate_and_covariance_unchanged() {
    let (r, a) = (radar(), ekf());
    let inn = a.innovation(&meas(8.0, 0.002, &a, &r), &r).unwrap();
    let mut pda = a.clone();
    pda_update(&mut pda, &r, &[(0.0, inn)], 1.0);
    assert!((a.x - pda.x).norm() < 1e-12 && (a.p - pda.p).norm() < 1e-9);
}

#[test]
fn uncertainty_about_which_measurement_is_real_inflates_the_covariance() {
    let (r, a) = (radar(), ekf());
    let i1 = a.innovation(&meas(40.0, 0.0, &a, &r), &r).unwrap();
    let i2 = a.innovation(&meas(-40.0, 0.0, &a, &r), &r).unwrap();
    let mut split = a.clone();
    pda_update(&mut split, &r, &[(0.5, i1), (0.5, i2)], 0.0);
    // the same two measurements fused as if one were certain gives a smaller covariance: the spread-of-means
    // term must make the split case more uncertain (positive semi-definite difference)
    let mut certain = a.clone();
    pda_update(&mut certain, &r, &[(1.0, i1)], 0.0);
    let diff = split.p - certain.p;
    let min_eig = nalgebra::SymmetricEigen::new(0.5 * (diff + diff.transpose()))
        .eigenvalues
        .min();
    assert!(
        min_eig > -1e-9,
        "covariance should not shrink when the association is uncertain: {min_eig}"
    );
    // symmetric measurements about the prediction cancel in the combined innovation
    assert!((split.x - a.x).norm() < 1e-6, "{}", (split.x - a.x).norm());
}

fn imm_at(x: State) -> fusetrack::imm::Imm {
    use fusetrack::imm::{Imm, ImmConfig};
    let p = Matrix4::from_diagonal(&nalgebra::Vector4::new(400.0, 400.0, 25.0, 25.0));
    Imm::new(
        ImmConfig::cv_ct(0.10, 0.5, 0.5, 0.9),
        x,
        p,
        vec![0.8, 0.1, 0.1],
    )
}

#[test]
fn imm_pda_with_one_certain_measurement_equals_the_plain_imm_update() {
    let r = radar();
    let x = State::new(4000.0, 1500.0, 20.0, -10.0);
    let mut plain = imm_at(x);
    let mut pda = plain.clone();
    plain.predict(1.0);
    pda.predict(1.0);
    let z = RadarMeas {
        range: r.h(&x).range + 12.0,
        bearing: r.h(&x).bearing + 0.002,
    };
    assert!(plain.update(&z, &r));
    assert!(pda.pda_update(&r, &[(1.0, z, 1.0e-6)], 0.0, 1.0));
    assert!((plain.estimate().0 - pda.estimate().0).norm() < 1e-9);
    assert!(plain
        .mu
        .iter()
        .zip(&pda.mu)
        .all(|(a, b)| (a - b).abs() < 1e-9));
}

#[test]
fn imm_pda_with_no_measurement_changes_nothing() {
    let x = State::new(4000.0, 1500.0, 20.0, -10.0);
    let mut imm = imm_at(x);
    imm.predict(1.0);
    let before = (imm.estimate().0, imm.mu.clone());
    assert!(!imm.pda_update(&radar(), &[], 1.0, 0.9));
    assert_eq!(before.0, imm.estimate().0);
    assert_eq!(before.1, imm.mu);
}

#[test]
fn imm_pda_shifts_probability_to_the_model_that_explains_a_turn() {
    // a target turning at 0.10 rad/s: after some scans of noise-free measurements, the turning model should dominate
    use fusetrack::ekf::Ekf;
    let r = radar();
    let mut truth = State::new(4000.0, 1500.0, 25.0, 0.0);
    let mut imm = imm_at(truth);
    let f = Ekf::ct_transition(0.10, 1.0);
    for _ in 0..25 {
        truth = f * truth;
        imm.predict(1.0);
        let h = r.h(&truth);
        assert!(imm.pda_update(
            &r,
            &[(
                1.0,
                RadarMeas {
                    range: h.range,
                    bearing: h.bearing
                },
                1.0e-6
            )],
            0.0,
            0.9
        ));
    }
    let turn_mass: f64 = imm.mu[1] + imm.mu[2];
    assert!(turn_mass > 0.8, "mu = {:?}", imm.mu);
    assert!((imm.estimate().0 - truth).norm() < 15.0);
}

#[test]
fn jpda_with_an_imm_tracks_a_manoeuvring_target_in_clutter() {
    use fusetrack::imm::ImmConfig;
    use fusetrack::sim::simulate_maneuvering;
    let r_max = 15000.0;
    let clutter_rate = 20.0;
    let density = clutter_rate / (std::f64::consts::PI * r_max * r_max);
    let cfg = ScenarioConfig {
        radar: radar(),
        dt: 1.0,
        steps: 100,
        n_targets: 1,
        sigma_accel: 0.5,
        pd: 0.9,
        clutter_rate,
        r_max,
    };
    let score = Some(ScoreConfig {
        pd: 0.9,
        clutter_density: density,
        confirm: 12.0,
        delete: -5.0,
    });
    let tcfg = TrackerConfig {
        imm: Some(ImmConfig::cv_ct(0.10, 0.5, 0.5, 0.9)),
        score,
        jpda: Some(JpdaConfig::new(0.9, density)),
        ..TrackerConfig::default()
    };
    let mut near = 0usize;
    let mut scored = 0usize;
    for seed in 0..10u64 {
        let mut rng = StdRng::seed_from_u64(seed);
        let (states, _) = simulate_maneuvering(
            State::new(4000.0, 1500.0, 25.0, 0.0),
            1.0,
            &[(0.0, 30), (0.08, 30), (0.0, 40)],
            0.5,
            &mut rng,
        );
        let scans = fusetrack::sim::simulate_from_trajectories(&cfg, &[states], &mut rng);
        let mut tracker = GnnTracker::new(tcfg.clone());
        for (k, scan) in scans.iter().enumerate() {
            tracker.step(cfg.dt, &scan.meas, &cfg.radar);
            if k < 25 {
                continue;
            }
            scored += 1;
            let t = Vector2::new(scan.truth[0][0], scan.truth[0][1]);
            if tracker
                .tracks
                .iter()
                .any(|tr| tr.status == TrackStatus::Confirmed && (tr.position() - t).norm() < 100.0)
            {
                near += 1;
            }
        }
    }
    assert!(
        near as f64 / scored as f64 > 0.9,
        "tracked {near} of {scored} scored scans"
    );
}

/// Average NEES of the confirmed track nearest the single true target, in clutter.
fn mean_nees(jpda: bool, clutter_rate: f64, runs: u64) -> (f64, usize) {
    let r_max = 15000.0;
    let cfg = ScenarioConfig {
        radar: radar(),
        dt: 1.0,
        steps: 100,
        n_targets: 1,
        sigma_accel: 1.0,
        pd: 0.9,
        clutter_rate,
        r_max,
    };
    let density = clutter_rate / (std::f64::consts::PI * r_max * r_max);
    let score = Some(ScoreConfig {
        pd: 0.9,
        clutter_density: density,
        confirm: 12.0,
        delete: -5.0,
    });
    let tcfg = TrackerConfig {
        sigma_accel: 1.0,
        score,
        jpda: jpda.then(|| JpdaConfig::new(0.9, density)),
        ..TrackerConfig::default()
    };
    let (mut sum, mut n) = (0.0, 0usize);
    for seed in 0..runs {
        let mut rng = StdRng::seed_from_u64(seed);
        let scans = simulate_scenario(&cfg, &mut rng);
        let mut tracker = GnnTracker::new(tcfg.clone());
        for (k, scan) in scans.iter().enumerate() {
            tracker.step(cfg.dt, &scan.meas, &cfg.radar);
            if k < 20 {
                continue;
            }
            let truth = scan.truth[0];
            let nearest: Option<&Track> = tracker
                .tracks
                .iter()
                .filter(|t| t.status == TrackStatus::Confirmed)
                .min_by(|a, b| {
                    (a.position() - Vector2::new(truth[0], truth[1]))
                        .norm()
                        .total_cmp(&(b.position() - Vector2::new(truth[0], truth[1])).norm())
                });
            if let Some(t) = nearest {
                if (t.position() - Vector2::new(truth[0], truth[1])).norm() > 100.0 {
                    continue; // a false track, not an estimate of this target
                }
                if let TrackFilter::Single(e) = &t.filter {
                    let err = truth - e.x;
                    if let Some(pinv) = e.p.try_inverse() {
                        sum += (err.transpose() * pinv * err)[0];
                        n += 1;
                    }
                }
            }
        }
    }
    (sum / n.max(1) as f64, n)
}

#[test]
fn jpda_estimates_stay_consistent_in_clutter() {
    // A consistent 4-state filter has mean NEES near 4. JPDA is an approximation, so it is allowed to drift
    // a little; the band only has to catch a gross error such as a wrong covariance update.
    for clutter in [0.0, 50.0] {
        let (j, nj) = mean_nees(true, clutter, 30);
        let (g, ng) = mean_nees(false, clutter, 30);
        eprintln!("clutter {clutter}: mean NEES JPDA {j:.4} ({nj} scans), GNN {g:.4} ({ng} scans); consistent = 4");
        assert!(nj > 500 && ng > 500, "too few scored scans");
        assert!(
            (2.0..8.0).contains(&j),
            "JPDA mean NEES {j} is far from 4 at clutter {clutter}"
        );
    }
}

#[test]
fn shallow_crossing_has_the_requested_relative_speed_and_miss_distance() {
    use fusetrack::sim::{simulate_crossing, simulate_crossing_angle};
    let cfg = ScenarioConfig {
        radar: radar(),
        dt: 1.0,
        steps: 100,
        n_targets: 2,
        sigma_accel: 0.0,
        pd: 1.0,
        clutter_rate: 0.0,
        r_max: 15000.0,
    };
    // 90 degrees reproduces the original scenario exactly for the same seed
    let a = simulate_crossing(&cfg, 40.0, 25.0, &mut StdRng::seed_from_u64(3));
    let b = simulate_crossing_angle(&cfg, 90.0, 40.0, 25.0, &mut StdRng::seed_from_u64(3));
    assert_eq!(a.len(), b.len());
    for (x, y) in a.iter().zip(&b) {
        assert_eq!(x.truth, y.truth);
    }
    // with no acceleration noise the relative speed is 2 v sin(angle / 2) and the closest approach is the set miss
    let s = simulate_crossing_angle(&cfg, 10.0, 40.0, 25.0, &mut StdRng::seed_from_u64(3));
    let rel = |k: usize| s[k].truth[1] - s[k].truth[0];
    let speed = (rel(10)[2].powi(2) + rel(10)[3].powi(2)).sqrt();
    assert!((speed - 2.0 * 25.0 * (5.0f64).to_radians().sin()).abs() < 1e-9);
    let closest = (0..s.len())
        .map(|k| (rel(k)[0].powi(2) + rel(k)[1].powi(2)).sqrt())
        .fold(f64::INFINITY, f64::min);
    assert!((closest - 40.0).abs() < 1.0, "closest approach {closest}");
}

#[test]
fn unresolved_returns_merge_to_their_mean_and_distant_ones_do_not() {
    use fusetrack::sim::{merge_unresolved, Scan};
    let m = |range: f64, bearing: f64| RadarMeas { range, bearing };
    let mk = |meas| Scan {
        truth: vec![],
        meas,
    };
    // two returns 20 m apart in range and 10 m in cross-range at 5 km merge; a far third stays
    let mut scans = vec![mk(vec![m(5000.0, 0.0), m(5020.0, 0.002), m(9000.0, 1.0)])];
    merge_unresolved(&mut scans, 30.0, 40.0);
    assert_eq!(scans[0].meas.len(), 2);
    let near = scans[0].meas.iter().find(|z| z.range < 6000.0).unwrap();
    assert!((near.range - 5010.0).abs() < 1e-9 && (near.bearing - 0.001).abs() < 1e-12);
    // a chain merges until no pair is left, and a zero cell changes nothing
    let mut chain = vec![mk(vec![m(5000.0, 0.0), m(5025.0, 0.0), m(5050.0, 0.0)])];
    merge_unresolved(&mut chain, 30.0, 40.0);
    assert!(chain[0].meas.len() <= 2);
    let mut same = vec![mk(vec![m(5000.0, 0.0), m(5001.0, 0.0)])];
    merge_unresolved(&mut same, 0.0, 0.0);
    assert_eq!(same[0].meas.len(), 2);
}
