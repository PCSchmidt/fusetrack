//! Multi-target GNN tracker on manoeuvring targets: single-model CV tracks at
//! two tunings vs. IMM tracks.
//!
//! Run with: `cargo run --release --example imm_tracker`

use fusetrack::eval::run_gnn_maneuvering_mean;
use fusetrack::sim::ScenarioConfig;
use fusetrack::{ImmConfig, Radar, TrackerConfig};
use nalgebra::Vector2;

fn main() {
    let n_seeds = 20;
    println!(
        "5 targets, 100 scans, random 25-scan segments (straight / +omega / -omega), Pd 0.9, {n_seeds} seeds"
    );
    println!("IMM: CV + CT(+/-0.10 rad/s), sigma 0.5, stay 0.9. M-of-N confirmation 3 of 5.\n");
    println!(
        "{:>6} {:>8} {:<20} {:>9} {:>6} {:>7} {:>8}",
        "omega", "clutter", "track filter", "mean OSPA", "false", "missed", "delay"
    );
    for omega in [0.10, 0.06, 0.14] {
        for clutter_rate in [0.0, 10.0] {
            let cfg = ScenarioConfig {
                radar: Radar {
                    pos: Vector2::new(0.0, 0.0),
                    sigma_r: 15.0,
                    sigma_b: 0.003,
                },
                dt: 1.0,
                steps: 100,
                n_targets: 5,
                sigma_accel: 0.5,
                pd: 0.9,
                clutter_rate,
                r_max: 15000.0,
            };
            let filters: Vec<(&str, TrackerConfig)> = vec![
                (
                    "CV EKF, sigma_a=0.5",
                    TrackerConfig {
                        sigma_accel: 0.5,
                        ..TrackerConfig::default()
                    },
                ),
                (
                    "CV EKF, sigma_a=3.0",
                    TrackerConfig {
                        sigma_accel: 3.0,
                        ..TrackerConfig::default()
                    },
                ),
                (
                    "IMM (CV, CT+, CT-)",
                    TrackerConfig {
                        imm: Some(ImmConfig::cv_ct(0.10, 0.5, 0.5, 0.9)),
                        ..TrackerConfig::default()
                    },
                ),
            ];
            for (name, t) in filters {
                let s = run_gnn_maneuvering_mean(&cfg, omega, 25, &t, n_seeds);
                println!(
                    "{:>6.2} {:>8.0} {:<20} {:>8.1}m {:>6.2} {:>7.2} {:>7.1}",
                    omega,
                    clutter_rate,
                    name,
                    s.mean_ospa,
                    s.mean_false,
                    s.mean_missed,
                    s.mean_confirm_delay
                );
            }
            println!();
        }
    }
}
