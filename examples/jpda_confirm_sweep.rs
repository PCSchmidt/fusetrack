//! Is JPDA's extra false-track rate in clutter just a confirmation-threshold effect?
//!
//! Run with: `cargo run --release --example jpda_confirm_sweep`
//!
//! Clutter 100 returns/scan, two targets side by side 120 m apart; GNN at confirm = 12 as the reference,
//! JPDA (no merge) at several score thresholds.

use fusetrack::eval::run_quality_mean;
use fusetrack::jpda::JpdaConfig;
use fusetrack::sim::{simulate_parallel, ScenarioConfig};
use fusetrack::{Radar, ScoreConfig, TrackerConfig};
use nalgebra::Vector2;

fn main() {
    let (pd, r_max, clutter_rate) = (0.9, 15000.0, 100.0);
    let density = clutter_rate / (std::f64::consts::PI * r_max * r_max);
    let cfg = ScenarioConfig {
        radar: Radar {
            pos: Vector2::new(0.0, 0.0),
            sigma_r: 15.0,
            sigma_b: 0.003,
        },
        dt: 1.0,
        steps: 100,
        n_targets: 2,
        sigma_accel: 1.0,
        pd,
        clutter_rate,
        r_max,
    };
    let score = |confirm: f64| {
        Some(ScoreConfig {
            pd,
            clutter_density: density,
            confirm,
            delete: -5.0,
        })
    };
    println!("clutter {clutter_rate:.0}, formation 120 m, 40 seeds\n");
    let run = |label: &str, t: TrackerConfig| {
        let s = run_quality_mean(&cfg, &t, 40, |rng| {
            simulate_parallel(&cfg, 120.0, 25.0, rng)
        });
        println!(
            "  {label:<22} OSPA {:>5.1} m  tracks {:>5.2}  false {:>5.2}  ID switches {:>4.2}",
            s.mean_ospa, s.mean_tracks, s.mean_false, s.id_switches
        );
    };
    run(
        "GNN confirm 12",
        TrackerConfig {
            score: score(12.0),
            ..TrackerConfig::default()
        },
    );
    for confirm in [12.0, 16.0, 20.0, 25.0, 30.0] {
        let jp = JpdaConfig {
            merge_duplicates: false,
            ..JpdaConfig::new(pd, density)
        };
        run(
            &format!("JPDA confirm {confirm:.0}"),
            TrackerConfig {
                score: score(confirm),
                jpda: Some(jp),
                ..TrackerConfig::default()
            },
        );
    }
    // and the birth rule: spawn more readily
    for birth in [0.2, 0.5, 0.8] {
        let jp = JpdaConfig {
            merge_duplicates: false,
            birth_threshold: birth,
            ..JpdaConfig::new(pd, density)
        };
        run(
            &format!("JPDA birth {birth:.1}"),
            TrackerConfig {
                score: score(12.0),
                jpda: Some(jp),
                ..TrackerConfig::default()
            },
        );
    }
}
