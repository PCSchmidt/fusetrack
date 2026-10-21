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

#[test]
#[should_panic(expected = "single-EKF")]
fn jpda_with_an_imm_is_refused() {
    use fusetrack::imm::ImmConfig;
    let imm = ImmConfig::cv_ct(0.05, 1.0, 1.0, 0.9);
    let _ = GnnTracker::new(TrackerConfig {
        imm: Some(imm),
        jpda: Some(JpdaConfig::new(0.9, 1e-7)),
        ..TrackerConfig::default()
    });
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
