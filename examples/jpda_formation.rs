//! JPDA against GNN on two targets flying side by side (same velocity), over clutter.
//!
//! Run with: `cargo run --release --example jpda_formation`
//!
//! Cross-range resolution at 6 km is about 18 m (range x sigma_b), so separations near that are hard for any
//! one-measurement-per-track rule. JPDA is shown with and without the duplicate-track merge, because the merge
//! rule treats two statistically indistinguishable tracks as one target and so fuses real formation members.

use fusetrack::eval::{run_quality_mean, QualityStats};
use fusetrack::jpda::JpdaConfig;
use fusetrack::sim::{simulate_parallel, ScenarioConfig};
use fusetrack::{Radar, ScoreConfig, TrackerConfig};
use nalgebra::Vector2;

const N_SEEDS: u64 = 40;
const R_MAX: f64 = 15000.0;
const PD: f64 = 0.9;

fn row(label: &str, s: &QualityStats) {
    println!(
        "    {:<12} OSPA {:>6.1} m  GOSPA {:>6.1} m  ID switches {:>5.2}  frag {:>5.2}  tracks {:>5.2}  false {:>5.2}  missed {:>5.2}",
        label, s.mean_ospa, s.mean_gospa, s.id_switches, s.fragmentations, s.mean_tracks, s.mean_false, s.mean_missed
    );
}

fn main() {
    println!("Two targets side by side at 25 m/s, Pd {PD}, {N_SEEDS} seeds, 100 scans\n");
    for clutter_rate in [0.0, 20.0, 100.0] {
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
        let jp = |merge: bool| JpdaConfig {
            merge_duplicates: merge,
            ..JpdaConfig::new(PD, density)
        };
        let gnn = TrackerConfig {
            score,
            ..TrackerConfig::default()
        };
        let jpda_merge = TrackerConfig {
            jpda: Some(jp(true)),
            score,
            ..TrackerConfig::default()
        };
        let jpda_free = TrackerConfig {
            jpda: Some(jp(false)),
            score,
            ..TrackerConfig::default()
        };
        println!("clutter {clutter_rate:.0} returns/scan");
        for sep in [30.0, 60.0, 120.0, 240.0] {
            println!("  separation {sep:.0} m");
            for (label, t) in [
                ("GNN", &gnn),
                ("JPDA", &jpda_free),
                ("JPDA+merge", &jpda_merge),
            ] {
                row(
                    label,
                    &run_quality_mean(&cfg, t, N_SEEDS, |rng| {
                        simulate_parallel(&cfg, sep, 25.0, rng)
                    }),
                );
            }
        }
        println!();
    }
}
