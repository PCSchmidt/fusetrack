//! Where do the extra tracks come from at high clutter, and what do the
//! fixes cost? Compares track-confirmation rules across clutter rate and Pd.
//!
//! Run with: `cargo run --release --example clutter_diagnose`

use std::f64::consts::PI;

use fusetrack::eval::run_gnn_mean;
use fusetrack::sim::ScenarioConfig;
use fusetrack::tracker::ScoreConfig;
use fusetrack::{Radar, TrackerConfig};
use nalgebra::Vector2;

const R_MAX: f64 = 15000.0;

fn main() {
    let n_seeds = 10;
    println!("5 targets, 100 scans, {n_seeds} seeds. false/missed = confirmed tracks / targets with no match within 100 m\n");
    println!(
        "{:>7} {:>4} {:<18} {:>9} {:>6} {:>7} {:>13}",
        "clutter", "Pd", "confirmation rule", "mean OSPA", "false", "missed", "confirm delay"
    );
    for pd in [0.9, 0.7] {
        for clutter_rate in [10.0, 50.0, 200.0] {
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
                r_max: R_MAX,
            };
            let density = clutter_rate / (PI * R_MAX * R_MAX);
            let score = |confirm: f64| ScoreConfig {
                pd,
                clutter_density: density,
                confirm,
                delete: -5.0,
            };
            let rules: Vec<(String, TrackerConfig)> = vec![
                ("3 of 5 (baseline)".into(), TrackerConfig::default()),
                (
                    "4 of 6".into(),
                    TrackerConfig {
                        confirm_hits: 4,
                        confirm_window: 6,
                        ..TrackerConfig::default()
                    },
                ),
                (
                    "5 of 7".into(),
                    TrackerConfig {
                        confirm_hits: 5,
                        confirm_window: 7,
                        ..TrackerConfig::default()
                    },
                ),
                (
                    "score >= 12".into(),
                    TrackerConfig {
                        score: Some(score(12.0)),
                        ..TrackerConfig::default()
                    },
                ),
                (
                    "score >= 20".into(),
                    TrackerConfig {
                        score: Some(score(20.0)),
                        ..TrackerConfig::default()
                    },
                ),
            ];
            for (name, t) in rules {
                let s = run_gnn_mean(&cfg, &t, n_seeds);
                println!(
                    "{:>7.0} {:>4.1} {:<18} {:>8.1}m {:>6.2} {:>7.2} {:>10.1} scans",
                    clutter_rate,
                    pd,
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
