//! End-to-end scenario evaluation: simulate, track, score with OSPA.

use nalgebra::{Matrix4, Vector2, Vector4};
use rand::{rngs::StdRng, SeedableRng};
use rand_distr::{Distribution, Normal};

use crate::ekf::Ekf;
use crate::imm::Imm;
use crate::metrics::ospa;
use crate::models::{Radar, RadarMeas, State};
use crate::sim::{simulate_maneuvering, simulate_scenario, ScenarioConfig};
use crate::tracker::{GnnTracker, TrackerConfig};

/// OSPA cutoff (m) and order used for all reported numbers.
pub const OSPA_CUTOFF: f64 = 200.0;
pub const OSPA_ORDER: f64 = 2.0;

/// Scans at the start of a run excluded from scoring while tracks form.
pub const BURN_IN: usize = 10;

#[derive(Debug, Clone, Copy)]
pub struct RunStats {
    /// Mean OSPA (m) over scored scans.
    pub mean_ospa: f64,
    /// Mean number of confirmed tracks per scored scan.
    pub mean_tracks: f64,
    /// Mean number of true targets per scored scan.
    pub mean_targets: f64,
}

/// Run one seeded scenario through the GNN tracker.
pub fn run_gnn(cfg: &ScenarioConfig, tracker_cfg: &TrackerConfig, seed: u64) -> RunStats {
    let mut rng = StdRng::seed_from_u64(seed);
    let scans = simulate_scenario(cfg, &mut rng);
    let mut tracker = GnnTracker::new(tracker_cfg.clone());

    let (mut ospa_sum, mut track_sum, mut target_sum, mut n) = (0.0, 0.0, 0.0, 0usize);
    for (k, scan) in scans.iter().enumerate() {
        tracker.step(cfg.dt, &scan.meas, &cfg.radar);
        if k < BURN_IN {
            continue;
        }
        let truth: Vec<Vector2<f64>> = scan
            .truth
            .iter()
            .map(|x| Vector2::new(x[0], x[1]))
            .collect();
        let est = tracker.confirmed_positions();
        ospa_sum += ospa(&truth, &est, OSPA_CUTOFF, OSPA_ORDER);
        track_sum += est.len() as f64;
        target_sum += truth.len() as f64;
        n += 1;
    }
    let n = n as f64;
    RunStats {
        mean_ospa: ospa_sum / n,
        mean_tracks: track_sum / n,
        mean_targets: target_sum / n,
    }
}

/// Average [`RunStats`] over `n_seeds` seeds starting at 0.
pub fn run_gnn_mean(cfg: &ScenarioConfig, tracker_cfg: &TrackerConfig, n_seeds: u64) -> RunStats {
    let runs: Vec<RunStats> = (0..n_seeds).map(|s| run_gnn(cfg, tracker_cfg, s)).collect();
    let k = runs.len() as f64;
    RunStats {
        mean_ospa: runs.iter().map(|r| r.mean_ospa).sum::<f64>() / k,
        mean_tracks: runs.iter().map(|r| r.mean_tracks).sum::<f64>() / k,
        mean_targets: runs.iter().map(|r| r.mean_targets).sum::<f64>() / k,
    }
}

/// A filter that tracks one target from range-bearing radar scans.
pub trait SingleTargetFilter {
    fn predict(&mut self, dt: f64);
    fn update(&mut self, z: &RadarMeas, radar: &Radar);
    fn position(&self) -> Vector2<f64>;
}

impl SingleTargetFilter for Ekf {
    fn predict(&mut self, dt: f64) {
        Ekf::predict(self, dt);
    }
    fn update(&mut self, z: &RadarMeas, radar: &Radar) {
        self.update_radar(z, radar);
    }
    fn position(&self) -> Vector2<f64> {
        Vector2::new(self.x[0], self.x[1])
    }
}

impl SingleTargetFilter for Imm {
    fn predict(&mut self, dt: f64) {
        Imm::predict(self, dt);
    }
    fn update(&mut self, z: &RadarMeas, radar: &Radar) {
        Imm::update(self, z, radar);
    }
    fn position(&self) -> Vector2<f64> {
        let (x, _) = self.estimate();
        Vector2::new(x[0], x[1])
    }
}

/// Manoeuvring-target benchmark: 40 steps straight, 40 turning at `+omega`,
/// 40 straight, 40 turning at `-omega`, at 40 m/s with one radar.
#[derive(Debug, Clone, Copy)]
pub struct ManeuverScenario {
    pub omega: f64,
    /// Truth acceleration noise (m/s^2).
    pub sigma_accel: f64,
}

/// Position RMSE (m) of a filter on the manoeuvre benchmark.
#[derive(Debug, Clone, Copy)]
pub struct ManeuverStats {
    pub rmse_all: f64,
    pub rmse_straight: f64,
    pub rmse_turn: f64,
}

pub const MANEUVER_RADAR: Radar = Radar {
    pos: Vector2::new(0.0, 0.0),
    sigma_r: 15.0,
    sigma_b: 0.003,
};

/// Run `make` (given an initial state and covariance) over `n_seeds` seeded
/// manoeuvre scenarios. Errors from the first 5 scans are excluded.
pub fn run_maneuver<F, M>(sc: ManeuverScenario, make: M, n_seeds: u64) -> ManeuverStats
where
    F: SingleTargetFilter,
    M: Fn(State, Matrix4<f64>) -> F,
{
    let radar = MANEUVER_RADAR;
    let dt = 1.0;
    let p0_std = Vector4::new(10.0, 10.0, 5.0, 5.0);
    let p0 = Matrix4::from_diagonal(&p0_std.component_mul(&p0_std));
    let segments = [(0.0, 40), (sc.omega, 40), (0.0, 40), (-sc.omega, 40)];

    let (mut all, mut straight, mut turn) = ((0.0, 0usize), (0.0, 0usize), (0.0, 0usize));
    for seed in 0..n_seeds {
        let mut rng = StdRng::seed_from_u64(seed);
        let x0 = State::new(5000.0, 2000.0, -28.0, 28.0);
        let (truth, omegas) = simulate_maneuvering(x0, dt, &segments, sc.sigma_accel, &mut rng);

        let unit = Normal::new(0.0, 1.0).unwrap();
        let init_err = Vector4::from_fn(|i, _| p0_std[i] * unit.sample(&mut rng));
        let mut filter = make(x0 + init_err, p0);

        for k in 1..truth.len() {
            filter.predict(dt);
            let z = radar.measure(&truth[k], &mut rng);
            filter.update(&z, &radar);
            if k <= 5 {
                continue;
            }
            let e = filter.position() - Vector2::new(truth[k][0], truth[k][1]);
            let sq = e.norm_squared();
            all.0 += sq;
            all.1 += 1;
            let bucket = if omegas[k - 1] == 0.0 {
                &mut straight
            } else {
                &mut turn
            };
            bucket.0 += sq;
            bucket.1 += 1;
        }
    }
    let rmse = |(sum, n): (f64, usize)| (sum / n as f64).sqrt();
    ManeuverStats {
        rmse_all: rmse(all),
        rmse_straight: rmse(straight),
        rmse_turn: rmse(turn),
    }
}
