//! Sweep clutter rate and detection probability; print mean OSPA.
//!
//! Run with: `cargo run --release --example clutter_sweep`

use fusetrack::eval::{run_gnn_mean, OSPA_CUTOFF};
use fusetrack::sim::ScenarioConfig;
use fusetrack::{Radar, TrackerConfig};
use nalgebra::Vector2;

fn main() {
    let n_seeds = 20;
    println!(
        "GNN tracker, 5 targets, 100 scans, {n_seeds} seeds, OSPA cutoff {OSPA_CUTOFF} m (p=2)\n"
    );
    println!(
        "{:>8} {:>6} {:>10} {:>14}",
        "clutter", "Pd", "mean OSPA", "tracks/targets"
    );
    for pd in [1.0, 0.9, 0.7] {
        for clutter_rate in [0.0, 10.0, 50.0, 200.0] {
            let cfg = ScenarioConfig {
                radar: Radar {
                    pos: Vector2::new(0.0, 0.0),
                    sigma_r: 15.0,
                    sigma_b: 0.003,
                },
                dt: 1.0,
                steps: 100,
                n_targets: 5,
                sigma_accel: 1.0,
                pd,
                clutter_rate,
                r_max: 15000.0,
            };
            let s = run_gnn_mean(&cfg, &TrackerConfig::default(), n_seeds);
            println!(
                "{:>8.0} {:>6.1} {:>9.1}m {:>8.2}/{:.0}",
                clutter_rate, pd, s.mean_ospa, s.mean_tracks, s.mean_targets
            );
        }
    }
}
