//! Compare a constant-velocity EKF at two tunings with an IMM on a
//! manoeuvring target. Prints position RMSE in metres.
//!
//! Run with: `cargo run --release --example imm_maneuver`

use fusetrack::eval::{run_maneuver, ManeuverScenario};
use fusetrack::{Ekf, Imm, ImmConfig};

fn main() {
    let n_seeds = 200;
    println!("Manoeuvring target (40 straight, 40 turn, 40 straight, 40 turn), {n_seeds} seeds");
    println!("Position RMSE (m); IMM models: CV + CT(+/-0.10 rad/s), sigma 0.5, stay 0.9\n");
    println!(
        "{:>14} {:<22} {:>8} {:>10} {:>8}",
        "true turn rate", "filter", "overall", "straight", "turning"
    );
    for omega in [0.10, 0.06, 0.14] {
        let sc = ManeuverScenario {
            omega,
            sigma_accel: 0.5,
        };
        let cv_low = run_maneuver(sc, |x, p| Ekf::new(x, p, 0.5), n_seeds);
        let cv_high = run_maneuver(sc, |x, p| Ekf::new(x, p, 3.0), n_seeds);
        let imm = run_maneuver(
            sc,
            |x, p| {
                Imm::new(
                    ImmConfig::cv_ct(0.10, 0.5, 0.5, 0.9),
                    x,
                    p,
                    vec![0.8, 0.1, 0.1],
                )
            },
            n_seeds,
        );
        for (name, s) in [
            ("CV EKF, sigma_a=0.5", cv_low),
            ("CV EKF, sigma_a=3.0", cv_high),
            ("IMM (CV, CT+, CT-)", imm),
        ] {
            println!(
                "{:>14.2} {:<22} {:>8.1} {:>10.1} {:>8.1}",
                omega, name, s.rmse_all, s.rmse_straight, s.rmse_turn
            );
        }
        println!();
    }
}
