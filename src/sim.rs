//! Scenario simulation: a target driven by piecewise-constant random
//! acceleration, matching the filter's process model.

use std::f64::consts::PI;

use rand::seq::SliceRandom;
use rand::Rng;
use rand_distr::{Distribution, Normal, Poisson};

use crate::ekf::Ekf;
use crate::models::{Radar, RadarMeas, State};

/// Advance a state by `dt` with a random acceleration `(ax, ay)` held constant
/// over the step.
pub fn step_truth(x: &State, dt: f64, ax: f64, ay: f64) -> State {
    State::new(
        x[0] + x[2] * dt + 0.5 * ax * dt * dt,
        x[1] + x[3] * dt + 0.5 * ay * dt * dt,
        x[2] + ax * dt,
        x[3] + ay * dt,
    )
}

/// Generate a truth trajectory of `steps + 1` states (including `x0`).
pub fn simulate_trajectory<R: Rng + ?Sized>(
    x0: State,
    dt: f64,
    steps: usize,
    sigma_accel: f64,
    rng: &mut R,
) -> Vec<State> {
    let noise = Normal::new(0.0, sigma_accel).expect("sigma_accel must be finite and >= 0");
    let mut out = Vec::with_capacity(steps + 1);
    out.push(x0);
    let mut x = x0;
    for _ in 0..steps {
        x = step_truth(&x, dt, noise.sample(rng), noise.sample(rng));
        out.push(x);
    }
    out
}

/// Generate a manoeuvring trajectory made of constant-turn-rate segments.
///
/// Each segment is `(omega, steps)`: the target turns at `omega` rad/s for
/// `steps` steps, with the same random-acceleration noise as
/// [`simulate_trajectory`]. Returns the states (including `x0`) and, for each
/// step, the turn rate in force.
pub fn simulate_maneuvering<R: Rng + ?Sized>(
    x0: State,
    dt: f64,
    segments: &[(f64, usize)],
    sigma_accel: f64,
    rng: &mut R,
) -> (Vec<State>, Vec<f64>) {
    let noise = Normal::new(0.0, sigma_accel).expect("sigma_accel must be finite and >= 0");
    let mut states = vec![x0];
    let mut omegas = Vec::new();
    let mut x = x0;
    for &(omega, steps) in segments {
        let f = Ekf::ct_transition(omega, dt);
        for _ in 0..steps {
            let (ax, ay) = (noise.sample(rng), noise.sample(rng));
            x = f * x + State::new(0.5 * ax * dt * dt, 0.5 * ay * dt * dt, ax * dt, ay * dt);
            states.push(x);
            omegas.push(omega);
        }
    }
    (states, omegas)
}

/// Parameters of a multi-target radar scenario with missed detections and
/// uniform clutter.
#[derive(Debug, Clone)]
pub struct ScenarioConfig {
    pub radar: Radar,
    pub dt: f64,
    pub steps: usize,
    pub n_targets: usize,
    /// Standard deviation of target acceleration (m/s^2).
    pub sigma_accel: f64,
    /// Probability that a target in range is detected on a scan.
    pub pd: f64,
    /// Mean number of clutter returns per scan (Poisson), uniform over the
    /// disc of radius `r_max` around the radar.
    pub clutter_rate: f64,
    pub r_max: f64,
}

/// One scan: the true target states and the (unlabelled, shuffled)
/// measurements the radar reported.
#[derive(Debug, Clone)]
pub struct Scan {
    pub truth: Vec<State>,
    pub meas: Vec<RadarMeas>,
}

/// Random starting states: 2-8 km from the radar at 15-35 m/s on random
/// headings.
fn random_starts<R: Rng + ?Sized>(cfg: &ScenarioConfig, rng: &mut R) -> Vec<State> {
    (0..cfg.n_targets)
        .map(|_| {
            let range = rng.gen_range(2000.0..8000.0);
            let bearing = rng.gen_range(-PI..PI);
            let speed = rng.gen_range(15.0..35.0);
            let heading = rng.gen_range(-PI..PI);
            State::new(
                cfg.radar.pos[0] + range * bearing.cos(),
                cfg.radar.pos[1] + range * bearing.sin(),
                speed * heading.cos(),
                speed * heading.sin(),
            )
        })
        .collect()
}

/// Turn truth trajectories (index 0 is the start state) into scans with
/// missed detections and clutter.
fn scans_from_trajectories<R: Rng + ?Sized>(
    cfg: &ScenarioConfig,
    trajectories: &[Vec<State>],
    rng: &mut R,
) -> Vec<Scan> {
    let clutter_count = (cfg.clutter_rate > 0.0)
        .then(|| Poisson::new(cfg.clutter_rate).expect("clutter_rate must be finite and > 0"));

    (1..=cfg.steps)
        .map(|k| {
            let truth: Vec<State> = trajectories.iter().map(|t| t[k]).collect();
            let mut meas = Vec::new();
            for x in &truth {
                let in_range = cfg.radar.h(x).range <= cfg.r_max;
                if in_range && rng.gen::<f64>() < cfg.pd {
                    meas.push(cfg.radar.measure(x, rng));
                }
            }
            let n_clutter = clutter_count.as_ref().map_or(0, |d| d.sample(rng) as usize);
            for _ in 0..n_clutter {
                // sqrt gives a uniform density over area.
                let range = cfg.r_max * rng.gen::<f64>().sqrt();
                meas.push(RadarMeas {
                    range,
                    bearing: rng.gen_range(-PI..PI),
                });
            }
            meas.shuffle(rng);
            Scan { truth, meas }
        })
        .collect()
}

/// Simulate `steps` scans of `n_targets` targets present throughout.
///
/// Targets start 2-8 km from the radar at 15-35 m/s on random headings and
/// are driven only by random acceleration noise (no deliberate manoeuvres).
pub fn simulate_scenario<R: Rng + ?Sized>(cfg: &ScenarioConfig, rng: &mut R) -> Vec<Scan> {
    let starts = random_starts(cfg, rng);
    let trajectories: Vec<Vec<State>> = starts
        .iter()
        .map(|x0| simulate_trajectory(*x0, cfg.dt, cfg.steps, cfg.sigma_accel, rng))
        .collect();
    scans_from_trajectories(cfg, &trajectories, rng)
}

/// Like [`simulate_scenario`], but each target flies a sequence of segments of
/// `segment_len` steps, each chosen at random from straight, a turn at
/// `+omega` or a turn at `-omega` (rad/s).
pub fn simulate_maneuvering_scenario<R: Rng + ?Sized>(
    cfg: &ScenarioConfig,
    omega: f64,
    segment_len: usize,
    rng: &mut R,
) -> Vec<Scan> {
    let starts = random_starts(cfg, rng);
    let n_segments = cfg.steps.div_ceil(segment_len);
    let trajectories: Vec<Vec<State>> = starts
        .iter()
        .map(|x0| {
            let segments: Vec<(f64, usize)> = (0..n_segments)
                .map(|_| {
                    let w = match rng.gen_range(0..3) {
                        0 => 0.0,
                        1 => omega,
                        _ => -omega,
                    };
                    (w, segment_len)
                })
                .collect();
            let (mut states, _) =
                simulate_maneuvering(*x0, cfg.dt, &segments, cfg.sigma_accel, rng);
            states.truncate(cfg.steps + 1);
            states
        })
        .collect();
    scans_from_trajectories(cfg, &trajectories, rng)
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::{rngs::StdRng, SeedableRng};

    #[test]
    fn zero_acceleration_is_constant_velocity() {
        let x0 = State::new(0.0, 0.0, 10.0, 5.0);
        let mut rng = StdRng::seed_from_u64(1);
        let traj = simulate_trajectory(x0, 1.0, 10, 0.0, &mut rng);
        assert_eq!(traj.len(), 11);
        assert!((traj[10][0] - 100.0).abs() < 1e-9);
        assert!((traj[10][1] - 50.0).abs() < 1e-9);
    }

    #[test]
    fn same_seed_same_trajectory() {
        let x0 = State::new(0.0, 0.0, 1.0, 1.0);
        let a = simulate_trajectory(x0, 1.0, 20, 2.0, &mut StdRng::seed_from_u64(7));
        let b = simulate_trajectory(x0, 1.0, 20, 2.0, &mut StdRng::seed_from_u64(7));
        assert_eq!(a, b);
    }

    fn cfg(pd: f64, clutter_rate: f64) -> ScenarioConfig {
        ScenarioConfig {
            radar: Radar {
                pos: nalgebra::Vector2::new(0.0, 0.0),
                sigma_r: 15.0,
                sigma_b: 0.003,
            },
            dt: 1.0,
            steps: 400,
            n_targets: 4,
            sigma_accel: 1.0,
            pd,
            clutter_rate,
            r_max: 50000.0,
        }
    }

    #[test]
    fn scenario_with_pd_one_and_no_clutter_has_one_measurement_per_target() {
        let scans = simulate_scenario(&cfg(1.0, 0.0), &mut StdRng::seed_from_u64(5));
        assert_eq!(scans.len(), 400);
        assert!(scans
            .iter()
            .all(|s| s.truth.len() == 4 && s.meas.len() == 4));
    }

    #[test]
    fn detection_and_clutter_rates_match_configuration() {
        let c = cfg(0.8, 20.0);
        let scans = simulate_scenario(&c, &mut StdRng::seed_from_u64(6));
        let mean = scans.iter().map(|s| s.meas.len() as f64).sum::<f64>() / scans.len() as f64;
        let expected = 0.8 * 4.0 + 20.0;
        assert!(
            (mean - expected).abs() < 1.0,
            "mean measurements {mean} vs expected {expected}"
        );
    }

    #[test]
    fn maneuvering_trajectory_follows_segments() {
        let x0 = State::new(0.0, 0.0, 40.0, 0.0);
        let segs = [(0.0, 10), (0.1, 20)];
        let (xs, om) = simulate_maneuvering(x0, 1.0, &segs, 0.0, &mut StdRng::seed_from_u64(1));
        assert_eq!(xs.len(), 31);
        assert_eq!(om.len(), 30);
        assert!(om[..10].iter().all(|&w| w == 0.0) && om[10..].iter().all(|&w| w == 0.1));
        // Straight leg: 10 steps at 40 m/s.
        assert!((xs[10][0] - 400.0).abs() < 1e-9 && xs[10][1].abs() < 1e-9);
        // Turn conserves speed and heading has rotated by 20 * 0.1 = 2 rad.
        let v = xs[30];
        assert!((v[2].hypot(v[3]) - 40.0).abs() < 1e-6);
        assert!((v[3].atan2(v[2]) - 2.0).abs() < 1e-6);
    }

    #[test]
    fn maneuvering_scenario_conserves_speed_and_actually_turns() {
        // Without acceleration noise, constant-turn-rate motion conserves speed
        // exactly, so every target must stay within its initial 15-35 m/s.
        let c = ScenarioConfig {
            sigma_accel: 0.0,
            steps: 100,
            ..cfg(1.0, 0.0)
        };
        let scans = simulate_maneuvering_scenario(&c, 0.1, 25, &mut StdRng::seed_from_u64(9));
        assert_eq!(scans.len(), 100);
        assert!(scans
            .iter()
            .all(|s| s.truth.len() == 4 && s.meas.len() == 4));
        for scan in &scans {
            for v in &scan.truth {
                let speed = v[2].hypot(v[3]);
                assert!(
                    (15.0 - 1e-6..=35.0 + 1e-6).contains(&speed),
                    "speed {speed}"
                );
            }
        }
        // At least one target changes heading over the run.
        let heading = |s: &Scan, i: usize| s.truth[i][3].atan2(s.truth[i][2]);
        let turned = (0..4).any(|i| {
            crate::models::wrap_angle(heading(&scans[99], i) - heading(&scans[0], i)).abs() > 0.1
        });
        assert!(turned, "no target turned");
    }
}
