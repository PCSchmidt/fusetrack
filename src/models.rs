//! Shared types and sensor models.

use std::f64::consts::{PI, TAU};

use nalgebra::{Vector2, Vector4};
use rand::Rng;
use rand_distr::{Distribution, Normal};

/// Target state: `[x, y, vx, vy]` in metres and metres per second.
pub type State = Vector4<f64>;

/// Wrap an angle to `[-pi, pi)`.
pub fn wrap_angle(a: f64) -> f64 {
    (a + PI).rem_euclid(TAU) - PI
}

/// A radar measurement: range (m) and bearing (rad, from +x axis, counter-clockwise).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RadarMeas {
    pub range: f64,
    pub bearing: f64,
}

/// A 2D range-bearing radar at a known position.
#[derive(Debug, Clone, Copy)]
pub struct Radar {
    pub pos: Vector2<f64>,
    /// Range noise standard deviation (m).
    pub sigma_r: f64,
    /// Bearing noise standard deviation (rad).
    pub sigma_b: f64,
}

impl Radar {
    /// Noise-free measurement of a state.
    pub fn h(&self, x: &State) -> RadarMeas {
        let dx = x[0] - self.pos[0];
        let dy = x[1] - self.pos[1];
        RadarMeas {
            range: dx.hypot(dy),
            bearing: dy.atan2(dx),
        }
    }

    /// Noisy measurement of a true state.
    pub fn measure<R: Rng + ?Sized>(&self, truth: &State, rng: &mut R) -> RadarMeas {
        let m = self.h(truth);
        let nr = Normal::new(0.0, self.sigma_r).expect("sigma_r must be finite and >= 0");
        let nb = Normal::new(0.0, self.sigma_b).expect("sigma_b must be finite and >= 0");
        RadarMeas {
            range: m.range + nr.sample(rng),
            bearing: wrap_angle(m.bearing + nb.sample(rng)),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn wrap_angle_maps_into_range() {
        for k in -5..=5 {
            let a = 0.3 + k as f64 * TAU;
            assert!((wrap_angle(a) - 0.3).abs() < 1e-9);
        }
        assert!((wrap_angle(PI + 0.1) - (-PI + 0.1)).abs() < 1e-12);
        assert!(wrap_angle(-PI) >= -PI && wrap_angle(PI) < PI);
    }

    #[test]
    fn radar_h_known_geometry() {
        let radar = Radar {
            pos: Vector2::new(0.0, 0.0),
            sigma_r: 1.0,
            sigma_b: 0.01,
        };
        let m = radar.h(&State::new(3.0, 4.0, 0.0, 0.0));
        assert!((m.range - 5.0).abs() < 1e-12);
        assert!((m.bearing - 4.0_f64.atan2(3.0)).abs() < 1e-12);
    }
}
