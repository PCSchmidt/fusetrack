//! End-to-end scenario evaluation: simulate, track, score with OSPA.

use nalgebra::Vector2;
use rand::{rngs::StdRng, SeedableRng};

use crate::metrics::ospa;
use crate::sim::{simulate_scenario, ScenarioConfig};
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
