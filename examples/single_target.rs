//! Track one simulated target with a single radar and print the error.
//!
//! Run with: `cargo run --release --example single_target`

use fusetrack::{sim, Ekf, Radar, State};
use nalgebra::{Matrix4, Vector2, Vector4};
use rand::{rngs::StdRng, SeedableRng};

fn main() {
    let mut rng = StdRng::seed_from_u64(42);
    let radar = Radar {
        pos: Vector2::new(0.0, 0.0),
        sigma_r: 15.0,
        sigma_b: 0.003,
    };
    let (dt, steps, sigma_accel) = (1.0, 60, 2.0);

    let x0 = State::new(6000.0, 3000.0, -25.0, 10.0);
    let truth = sim::simulate_trajectory(x0, dt, steps, sigma_accel, &mut rng);

    let p0 = Matrix4::from_diagonal(&Vector4::new(100.0, 100.0, 25.0, 25.0));
    let mut ekf = Ekf::new(x0, p0, sigma_accel);

    println!(
        "{:>4} {:>10} {:>10} {:>10}",
        "t", "truth_x", "est_x", "pos_err"
    );
    for (k, tk) in truth.iter().enumerate().skip(1) {
        ekf.predict(dt);
        let z = radar.measure(tk, &mut rng);
        ekf.update_radar(&z, &radar);
        if k % 10 == 0 {
            let err = ((ekf.x[0] - tk[0]).powi(2) + (ekf.x[1] - tk[1]).powi(2)).sqrt();
            println!("{:>4} {:>10.1} {:>10.1} {:>10.2}", k, tk[0], ekf.x[0], err);
        }
    }
}
