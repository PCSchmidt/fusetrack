//! JPDA against GNN when a finite-resolution radar merges close returns into one detection.
//!
//! Run with: `cargo run --release --example jpda_merged`
//!
//! Two scenarios (shallow crossing at 0 m closest approach, and side by side) with a resolution cell of none,
//! 30 m x 40 m, or 60 m x 80 m (range x cross-range). Both trackers see identical merged scans each seed;
//! differences are paired (JPDA minus GNN, mean with a 95% interval of 1.96 standard errors).

use fusetrack::eval::{evaluate_quality, QualityStats};
use fusetrack::jpda::JpdaConfig;
use fusetrack::sim::{
    merge_unresolved, simulate_crossing_angle, simulate_parallel, Scan, ScenarioConfig,
};
use fusetrack::{Radar, ScoreConfig, TrackerConfig};
use nalgebra::Vector2;
use rand::{rngs::StdRng, SeedableRng};

const N_SEEDS: u64 = 60;
const R_MAX: f64 = 15000.0;
const SPEED: f64 = 25.0;
const PD: f64 = 0.9;

fn mean_ci(v: &[f64]) -> (f64, f64) {
    let n = v.len() as f64;
    let m = v.iter().sum::<f64>() / n;
    let var = v.iter().map(|x| (x - m).powi(2)).sum::<f64>() / (n - 1.0);
    (m, 1.96 * (var / n).sqrt())
}

fn main() {
    println!("Finite resolution: close returns merge into one detection. Pd {PD}, {SPEED} m/s, {N_SEEDS} paired seeds, 100 scans");
    println!("Difference = JPDA - GNN; negative is better for every metric.\n");
    for clutter_rate in [0.0, 20.0] {
        let density =
            (clutter_rate / (std::f64::consts::PI * R_MAX * R_MAX)).max(1.0e-9 / (R_MAX * R_MAX));
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
        let score = ScoreConfig {
            pd: PD,
            clutter_density: density,
            confirm: 12.0,
            delete: -5.0,
        };
        let gnn = TrackerConfig {
            score: Some(score),
            ..TrackerConfig::default()
        };
        let jpda = TrackerConfig {
            jpda: Some(JpdaConfig::new(PD, density)),
            score: Some(score),
            ..TrackerConfig::default()
        };
        println!("clutter {clutter_rate:.0} returns/scan");
        type Gen = Box<dyn Fn(&mut StdRng) -> Vec<Scan>>;
        let mut scenarios: Vec<(String, Gen)> = Vec::new();
        for angle in [5.0, 10.0, 20.0] {
            let c = cfg.clone();
            scenarios.push((
                format!("crossing {angle:>3.0} deg, 0 m"),
                Box::new(move |r| simulate_crossing_angle(&c, angle, 0.0, SPEED, r)),
            ));
        }
        for sep in [30.0, 60.0, 120.0] {
            let c = cfg.clone();
            scenarios.push((
                format!("side by side {sep:>3.0} m   "),
                Box::new(move |r| simulate_parallel(&c, sep, SPEED, r)),
            ));
        }
        for (label, gen) in &scenarios {
            for (rr, rx) in [(0.0, 0.0), (30.0, 40.0), (60.0, 80.0)] {
                let (mut g, mut j): (Vec<QualityStats>, Vec<QualityStats>) = (vec![], vec![]);
                let mut merged_scans = 0.0;
                for seed in 0..N_SEEDS {
                    let mut rng = StdRng::seed_from_u64(seed);
                    let mut scans = gen(&mut rng);
                    let before: usize = scans.iter().map(|s| s.meas.len()).sum();
                    merge_unresolved(&mut scans, rr, rx);
                    let after: usize = scans.iter().map(|s| s.meas.len()).sum();
                    merged_scans += (before - after) as f64;
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
                println!(
                    "  {label} cell {rr:>2.0}x{rx:>2.0} m | merges/run {:>5.1} | GOSPA GNN {:>6.1} JPDA {:>6.1} diff {o:>+6.1} ({oc:.1}) | ID sw GNN {:>4.2} JPDA {:>4.2} diff {i:>+5.2} ({ic:.2}) | missed GNN {:>4.2} JPDA {:>4.2} | frag GNN {:>4.2} JPDA {:>4.2}",
                    merged_scans / N_SEEDS as f64,
                    gm(|s| s.mean_gospa),
                    jm(|s| s.mean_gospa),
                    gm(|s| s.id_switches),
                    jm(|s| s.id_switches),
                    gm(|s| s.mean_missed),
                    jm(|s| s.mean_missed),
                    gm(|s| s.fragmentations),
                    jm(|s| s.fragmentations),
                );
            }
        }
        println!();
    }
}
