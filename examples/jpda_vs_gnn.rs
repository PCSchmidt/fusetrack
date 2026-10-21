//! JPDA against GNN on two targets crossing with a controlled closest approach, over clutter.
//!
//! Run with: `cargo run --release --example jpda_vs_gnn`
//!
//! Both trackers use the same score-based confirmation (confirm at 12, delete at -5) and the true Pd and
//! clutter density. Each cell averages `N_SEEDS` seeded runs of 100 scans.

use fusetrack::eval::{run_crossing_mean, QualityStats};
use fusetrack::jpda::JpdaConfig;
use fusetrack::sim::ScenarioConfig;
use fusetrack::{Radar, ScoreConfig, TrackerConfig};
use nalgebra::Vector2;

const N_SEEDS: u64 = 40;
/// Heaviest clutter makes JPDA slow (large clusters), so it gets fewer seeds.
const N_SEEDS_HEAVY: u64 = 15;
const R_MAX: f64 = 15000.0;
const PD: f64 = 0.9;

fn row(label: &str, s: &QualityStats) {
    println!(
        "    {:<5} OSPA {:>6.1} m  GOSPA {:>6.1} m  ID switches {:>5.2}  fragmentations {:>5.2}  tracks {:>5.2}  false {:>5.2}  missed {:>5.2}  inexact scans {:>5.1}",
        label, s.mean_ospa, s.mean_gospa, s.id_switches, s.fragmentations, s.mean_tracks, s.mean_false, s.mean_missed, s.inexact_scans
    );
}

fn main() {
    println!("Two targets crossing at 90 degrees, 25 m/s, Pd {PD}, {N_SEEDS} seeds, 100 scans\n");
    for clutter_rate in [0.0, 20.0, 100.0, 300.0] {
        let density = clutter_rate / (std::f64::consts::PI * R_MAX * R_MAX);
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
            pd: PD,
            clutter_rate,
            r_max: R_MAX,
        };
        let score = (clutter_rate > 0.0).then_some(ScoreConfig {
            pd: PD,
            clutter_density: density,
            confirm: 12.0,
            delete: -5.0,
        });
        let gnn = TrackerConfig {
            score,
            ..TrackerConfig::default()
        };
        // with no clutter the PDA weight needs a density; use a tiny one so measurements dominate the miss hypothesis
        let jp_density = if clutter_rate > 0.0 { density } else { 1.0e-12 };
        let jpda = TrackerConfig {
            jpda: Some(JpdaConfig::new(PD, jp_density)),
            score,
            ..TrackerConfig::default()
        };
        let seeds = if clutter_rate >= 300.0 {
            N_SEEDS_HEAVY
        } else {
            N_SEEDS
        };
        println!("clutter {clutter_rate:.0} returns/scan");
        for miss in [25.0, 50.0, 100.0, 200.0, 500.0] {
            println!("  closest approach {miss:.0} m");
            row("GNN", &run_crossing_mean(&cfg, miss, 25.0, &gnn, seeds));
            row("JPDA", &run_crossing_mean(&cfg, miss, 25.0, &jpda, seeds));
        }
        println!();
    }
}
