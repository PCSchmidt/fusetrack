//! JPDA and GNN, each with and without the IMM, on two targets crossing at a shallow angle where one of them turns
//! around the time they are closest. Every tracker sees the same scans on each seed, so differences are paired.
//!
//! Run with: `cargo run --release --example jpda_imm`
//!
//! Prints GOSPA and ID switches per run for each tracker and paired differences (mean, 95% interval of 1.96
//! standard errors) for the three comparisons that matter: IMM's effect under GNN, IMM's effect under JPDA, and
//! JPDA's effect given the IMM.

use fusetrack::eval::{evaluate_quality, QualityStats};
use fusetrack::imm::ImmConfig;
use fusetrack::jpda::JpdaConfig;
use fusetrack::sim::{simulate_crossing_turn, ScenarioConfig};
use fusetrack::{Radar, ScoreConfig, TrackerConfig};
use nalgebra::Vector2;
use rand::{rngs::StdRng, SeedableRng};

const N_SEEDS: u64 = 40;
const R_MAX: f64 = 15000.0;
const SPEED: f64 = 25.0;

fn mean_ci(v: &[f64]) -> (f64, f64) {
    let n = v.len() as f64;
    let m = v.iter().sum::<f64>() / n;
    let var = v.iter().map(|x| (x - m).powi(2)).sum::<f64>() / (n - 1.0);
    (m, 1.96 * (var / n).sqrt())
}

fn main() {
    println!("Two targets crossing, one turning for 20 s from step 40; {SPEED} m/s, {N_SEEDS} paired seeds, 100 scans");
    println!(
        "Differences are second minus first; negative is better for GOSPA and ID switches.
"
    );
    for (omega, label) in [
        (
            0.04,
            "gentle turn, 0.04 rad/s (about 1 m/s^2 sideways, 46 degrees)",
        ),
        (
            0.10,
            "hard turn, 0.10 rad/s (about 2.5 m/s^2 sideways, 115 degrees)",
        ),
    ] {
        println!(
            "=== {label}
"
        );
        for pd in [0.9, 0.7] {
            for clutter_rate in [0.0, 20.0, 60.0] {
                let density = (clutter_rate / (std::f64::consts::PI * R_MAX * R_MAX))
                    .max(1.0e-9 / (R_MAX * R_MAX));
                let cfg = ScenarioConfig {
                    radar: Radar {
                        pos: Vector2::new(0.0, 0.0),
                        sigma_r: 15.0,
                        sigma_b: 0.003,
                    },
                    dt: 1.0,
                    steps: 100,
                    n_targets: 2,
                    sigma_accel: 0.5,
                    pd,
                    clutter_rate,
                    r_max: R_MAX,
                };
                let score = ScoreConfig {
                    pd,
                    clutter_density: density,
                    confirm: 12.0,
                    delete: -5.0,
                };
                let imm = Some(ImmConfig::cv_ct(0.10, 0.5, 0.5, 0.9));
                let mk = |jpda: bool, imm: &Option<ImmConfig>| TrackerConfig {
                    score: Some(score),
                    imm: imm.clone(),
                    jpda: jpda.then(|| JpdaConfig::new(pd, density)),
                    ..TrackerConfig::default()
                };
                let trackers = [
                    ("GNN", mk(false, &None)),
                    ("GNN+IMM", mk(false, &imm)),
                    ("JPDA", mk(true, &None)),
                    ("JPDA+IMM", mk(true, &imm)),
                ];
                println!("Pd {pd}, clutter {clutter_rate:.0} returns/scan");
                for angle in [10.0, 20.0] {
                    for miss in [0.0, 40.0] {
                        let mut runs: Vec<Vec<QualityStats>> = vec![vec![]; trackers.len()];
                        for seed in 0..N_SEEDS {
                            let mut rng = StdRng::seed_from_u64(seed);
                            let scans = simulate_crossing_turn(
                                &cfg, angle, miss, SPEED, omega, 40, 20, &mut rng,
                            );
                            for (k, (_, t)) in trackers.iter().enumerate() {
                                runs[k].push(evaluate_quality(&cfg, &scans, t));
                            }
                        }
                        let m = |k: usize, f: fn(&QualityStats) -> f64| {
                            runs[k].iter().map(f).sum::<f64>() / N_SEEDS as f64
                        };
                        let d = |a: usize, b: usize, f: fn(&QualityStats) -> f64| -> (f64, f64) {
                            mean_ci(
                                &runs[a]
                                    .iter()
                                    .zip(&runs[b])
                                    .map(|(x, y)| f(y) - f(x))
                                    .collect::<Vec<_>>(),
                            )
                        };
                        let g = |s: &QualityStats| s.mean_gospa;
                        let i = |s: &QualityStats| s.id_switches;
                        print!("  angle {angle:>3.0} deg, closest {miss:>3.0} m | GOSPA");
                        for (k, (name, _)) in trackers.iter().enumerate() {
                            print!(" {name} {:>5.1}", m(k, g));
                        }
                        print!(" | ID sw");
                        for (k, (name, _)) in trackers.iter().enumerate() {
                            print!(" {name} {:>4.2}", m(k, i));
                        }
                        println!();
                        for (label, a, b) in [
                            ("IMM under GNN (GNN -> GNN+IMM)", 0, 1),
                            ("IMM under JPDA (JPDA -> JPDA+IMM)", 2, 3),
                            ("JPDA given IMM (GNN+IMM -> JPDA+IMM)", 1, 3),
                        ] {
                            let (dg, cg) = d(a, b, g);
                            let (di, ci) = d(a, b, i);
                            println!("      {label:<38} GOSPA {dg:>+6.1} ({cg:.1})   ID switches {di:>+5.2} ({ci:.2})");
                        }
                    }
                }
                println!();
            }
        }
    }
}
