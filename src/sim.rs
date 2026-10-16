//! Scenario simulation: a target driven by piecewise-constant random
//! acceleration, matching the filter's process model.

use rand::Rng;
use rand_distr::{Distribution, Normal};

use crate::models::State;

/// Advance a state by `dt` with a random acceleration `(ax, ay)` held constant
/// over the step.
pub fn step_truth(x: &State, dt: f64, ax: f64, ay: f64) -> State {
    State::new(
        x[0] + x[2] * dt + 0.5 * ax * dt * dt,
        x[1] + x[3] * dt + 0.5 * ay * dt * dt,
        x[2] + ax * dt,
        x[3] + ay * dt,
    )
}

/// Generate a truth trajectory of `steps + 1` states (including `x0`).
pub fn simulate_trajectory<R: Rng + ?Sized>(
    x0: State,
    dt: f64,
    steps: usize,
    sigma_accel: f64,
    rng: &mut R,
) -> Vec<State> {
    let noise = Normal::new(0.0, sigma_accel).expect("sigma_accel must be finite and >= 0");
    let mut out = Vec::with_capacity(steps + 1);
    out.push(x0);
    let mut x = x0;
    for _ in 0..steps {
        x = step_truth(&x, dt, noise.sample(rng), noise.sample(rng));
        out.push(x);
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::{rngs::StdRng, SeedableRng};

    #[test]
    fn zero_acceleration_is_constant_velocity() {
        let x0 = State::new(0.0, 0.0, 10.0, 5.0);
        let mut rng = StdRng::seed_from_u64(1);
        let traj = simulate_trajectory(x0, 1.0, 10, 0.0, &mut rng);
        assert_eq!(traj.len(), 11);
        assert!((traj[10][0] - 100.0).abs() < 1e-9);
        assert!((traj[10][1] - 50.0).abs() < 1e-9);
    }

    #[test]
    fn same_seed_same_trajectory() {
        let x0 = State::new(0.0, 0.0, 1.0, 1.0);
        let a = simulate_trajectory(x0, 1.0, 20, 2.0, &mut StdRng::seed_from_u64(7));
        let b = simulate_trajectory(x0, 1.0, 20, 2.0, &mut StdRng::seed_from_u64(7));
        assert_eq!(a, b);
    }
}
