//! IMM vs. constant-velocity EKF on the manoeuvring-target benchmark.

use fusetrack::eval::{run_maneuver, ManeuverScenario};
use fusetrack::{Ekf, Imm, ImmConfig, State};
use nalgebra::Matrix4;

fn imm(x: State, p: Matrix4<f64>) -> Imm {
    Imm::new(
        ImmConfig::cv_ct(0.10, 0.5, 0.5, 0.9),
        x,
        p,
        vec![0.8, 0.1, 0.1],
    )
}

#[test]
fn imm_beats_both_cv_tunings_when_the_turn_rate_is_matched() {
    let sc = ManeuverScenario {
        omega: 0.10,
        sigma_accel: 0.5,
    };
    let cv_low = run_maneuver(sc, |x, p| Ekf::new(x, p, 0.5), 50);
    let cv_high = run_maneuver(sc, |x, p| Ekf::new(x, p, 3.0), 50);
    let imm = run_maneuver(sc, imm, 50);
    println!("low {cv_low:?}\nhigh {cv_high:?}\nimm {imm:?}");
    // A low-noise CV filter lags badly in turns; a high-noise one is noisy on
    // straight legs. The IMM should beat both overall.
    assert!(imm.rmse_all < cv_low.rmse_all);
    assert!(imm.rmse_all < cv_high.rmse_all);
    assert!(imm.rmse_turn < cv_low.rmse_turn);
    assert!(imm.rmse_straight < cv_low.rmse_straight);
}

#[test]
fn imm_stays_comparable_to_a_tuned_cv_when_the_turn_rate_is_mismatched() {
    let sc = ManeuverScenario {
        omega: 0.14,
        sigma_accel: 0.5,
    };
    let cv_high = run_maneuver(sc, |x, p| Ekf::new(x, p, 3.0), 50);
    let imm = run_maneuver(sc, imm, 50);
    println!("high {cv_high:?}\nimm {imm:?}");
    assert!(imm.rmse_all < cv_high.rmse_all * 1.15);
}
