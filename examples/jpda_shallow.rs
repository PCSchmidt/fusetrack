//! JPDA against GNN on two targets that cross at a shallow angle (similar velocities), at several detection
//! probabilities and clutter levels. Both trackers see identical scans each seed, so differences are paired.
//!
//! Run with: `cargo run --release --example jpda_shallow`
//!
//! Prints the mean paired difference (JPDA minus GNN) with a 95% interval (1.96 standard errors) over seeds.

use fusetrack::eval::{evaluate_quality, QualityStats};
use fusetrack::jpda::JpdaConfig;
use fusetrack::sim::{simulate_crossing_angle, ScenarioConfig};
use fusetrack::{Radar, ScoreConfig, TrackerConfig};
use nalgebra::Vector2;
use rand::{rngs::StdRng, SeedableRng};

const N_SEEDS: u64 = 60;
const R_MAX: f64 = 15000.0;
const SPEED: f64 = 25.0;

fn mean_ci(v: &[f64]) -> (f64, f64) {
    let n = v.len() as f64;
    let m = v.iter().sum::<f64>() / n;
    let var = v.iter().map(|x| (x - m).powi(2)).sum::<f64>() / (n - 1.0);
    (m, 1.96 * (var / n).sqrt())
}

fn main() {
    println!(
        "Two targets crossing at a shallow angle, {SPEED} m/s, {N_SEEDS} paired seeds, 100 scans"
    );
    println!("cells: mean (95% interval). Difference = JPDA - GNN; negative is better for every metric.\n");
    for pd in [0.9, 0.7] {
        for clutter_rate in [0.0, 20.0] {
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
                pd,
                clutter_rate,
                r_max: R_MAX,
            };
            let score = ScoreConfig {
                pd,
                clutter_density: density.max(1.0e-9 / (R_MAX * R_MAX)),
                confirm: 12.0,
                delete: -5.0,
            };
            let gnn = TrackerConfig {
                score: Some(score),
                ..TrackerConfig::default()
            };
            let jpda = TrackerConfig {
                jpda: Some(JpdaConfig::new(pd, score.clutter_density)),
                score: Some(score),
                ..TrackerConfig::default()
            };
            println!("Pd {pd}, clutter {clutter_rate:.0} returns/scan");
            for angle in [5.0, 10.0, 20.0, 40.0] {
                for miss in [0.0, 40.0] {
                    let (mut g, mut j): (Vec<QualityStats>, Vec<QualityStats>) = (vec![], vec![]);
                    for seed in 0..N_SEEDS {
                        let mut rng = StdRng::seed_from_u64(seed);
                        let scans = simulate_crossing_angle(&cfg, angle, miss, SPEED, &mut rng);
                        g.push(evaluate_quality(&cfg, &scans, &gnn));
                        j.push(evaluate_quality(&cfg, &scans, &jpda));
                    }
                    let d = |f: fn(&QualityStats) -> f64| -> Vec<f64> {
                        g.iter().zip(&j).map(|(a, b)| f(b) - f(a)).collect()
                    };
                    let gm =
                        |f: fn(&QualityStats) -> f64| g.iter().map(f).sum::<f64>() / N_SEEDS as f64;
                    let jm =
                        |f: fn(&QualityStats) -> f64| j.iter().map(f).sum::<f64>() / N_SEEDS as f64;
                    let (o, oc) = mean_ci(&d(|s| s.mean_gospa));
                    let (i, ic) = mean_ci(&d(|s| s.id_switches));
                    let (m, mc) = mean_ci(&d(|s| s.mean_missed));
                    println!(
                        "  angle {angle:>4.0} deg, closest {miss:>3.0} m | GOSPA GNN {:>6.1} JPDA {:>6.1} diff {o:>+6.1} ({oc:.1}) | ID sw GNN {:>4.2} JPDA {:>4.2} diff {i:>+5.2} ({ic:.2}) | missed diff {m:>+5.2} ({mc:.2}) | inexact {:.1}",
                        gm(|s| s.mean_gospa),
                        jm(|s| s.mean_gospa),
                        gm(|s| s.id_switches),
                        jm(|s| s.id_switches),
                        jm(|s| s.inexact_scans),
                    );
                }
            }
            println!();
        }
    }
}
