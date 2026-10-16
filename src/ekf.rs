//! Extended Kalman filter with a constant-velocity motion model and a
//! range-bearing radar update.

use nalgebra::{Matrix2, Matrix2x4, Matrix4, Matrix4x2};

use crate::models::{wrap_angle, Radar, RadarMeas, State};

/// Single-target EKF.
#[derive(Debug, Clone)]
pub struct Ekf {
    pub x: State,
    pub p: Matrix4<f64>,
    /// Standard deviation of the unmodelled acceleration (m/s^2).
    pub sigma_accel: f64,
}

impl Ekf {
    pub fn new(x: State, p: Matrix4<f64>, sigma_accel: f64) -> Self {
        Self { x, p, sigma_accel }
    }

    /// Constant-velocity state transition for a step of `dt` seconds.
    pub fn transition(dt: f64) -> Matrix4<f64> {
        Matrix4::new(
            1.0, 0.0, dt, 0.0, //
            0.0, 1.0, 0.0, dt, //
            0.0, 0.0, 1.0, 0.0, //
            0.0, 0.0, 0.0, 1.0,
        )
    }

    /// Process noise for piecewise-constant acceleration over `dt` seconds.
    pub fn process_noise(dt: f64, sigma_accel: f64) -> Matrix4<f64> {
        let (dt2, dt3, dt4) = (dt * dt, dt * dt * dt, dt * dt * dt * dt);
        sigma_accel.powi(2)
            * Matrix4::new(
                dt4 / 4.0,
                0.0,
                dt3 / 2.0,
                0.0, //
                0.0,
                dt4 / 4.0,
                0.0,
                dt3 / 2.0, //
                dt3 / 2.0,
                0.0,
                dt2,
                0.0, //
                0.0,
                dt3 / 2.0,
                0.0,
                dt2,
            )
    }

    pub fn predict(&mut self, dt: f64) {
        let f = Self::transition(dt);
        let q = Self::process_noise(dt, self.sigma_accel);
        self.x = f * self.x;
        self.p = f * self.p * f.transpose() + q;
    }

    /// Fuse a radar measurement. Returns the normalised innovation squared
    /// (NIS), or `None` if the update was skipped (target on top of the
    /// sensor, or singular innovation covariance).
    pub fn update_radar(&mut self, z: &RadarMeas, radar: &Radar) -> Option<f64> {
        let dx = self.x[0] - radar.pos[0];
        let dy = self.x[1] - radar.pos[1];
        let r2 = dx * dx + dy * dy;
        if r2 < 1e-12 {
            return None;
        }
        let r = r2.sqrt();

        let pred = radar.h(&self.x);
        let y = nalgebra::Vector2::new(z.range - pred.range, wrap_angle(z.bearing - pred.bearing));

        let h = Matrix2x4::new(
            dx / r,
            dy / r,
            0.0,
            0.0, //
            -dy / r2,
            dx / r2,
            0.0,
            0.0,
        );
        let rm = Matrix2::new(
            radar.sigma_r.powi(2),
            0.0, //
            0.0,
            radar.sigma_b.powi(2),
        );

        let s = h * self.p * h.transpose() + rm;
        let s_inv = s.try_inverse()?;
        let k: Matrix4x2<f64> = self.p * h.transpose() * s_inv;

        self.x += k * y;
        // Joseph form keeps P symmetric positive semi-definite.
        let i_kh = Matrix4::identity() - k * h;
        self.p = i_kh * self.p * i_kh.transpose() + k * rm * k.transpose();

        Some((y.transpose() * s_inv * y)[0])
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::Vector2;

    fn radar() -> Radar {
        Radar {
            pos: Vector2::new(0.0, 0.0),
            sigma_r: 10.0,
            sigma_b: 0.005,
        }
    }

    #[test]
    fn predict_advances_position_and_grows_covariance() {
        let mut ekf = Ekf::new(
            State::new(100.0, 200.0, 10.0, -5.0),
            Matrix4::identity() * 100.0,
            1.0,
        );
        let trace_before = ekf.p.trace();
        ekf.predict(2.0);
        assert!((ekf.x[0] - 120.0).abs() < 1e-12);
        assert!((ekf.x[1] - 190.0).abs() < 1e-12);
        assert!(ekf.p.trace() > trace_before);
    }

    #[test]
    fn perfect_measurement_shrinks_covariance_and_keeps_it_symmetric() {
        let truth = State::new(3000.0, 4000.0, 0.0, 0.0);
        let mut ekf = Ekf::new(truth, Matrix4::identity() * 1.0e4, 1.0);
        let z = radar().h(&truth);
        let trace_before = ekf.p.trace();
        let nis = ekf.update_radar(&z, &radar()).expect("update applied");
        assert!(ekf.p.trace() < trace_before);
        assert!((ekf.p - ekf.p.transpose()).norm() < 1e-9);
        assert!(nis.abs() < 1e-9, "zero innovation expected, got {nis}");
    }

    #[test]
    fn update_skipped_when_target_on_sensor() {
        let mut ekf = Ekf::new(State::zeros(), Matrix4::identity(), 1.0);
        let z = RadarMeas {
            range: 0.0,
            bearing: 0.0,
        };
        assert!(ekf.update_radar(&z, &radar()).is_none());
    }

    #[test]
    fn bearing_innovation_wraps_across_the_branch_cut() {
        // Target just below the -x axis (bearing ~ -pi + eps); measurement just
        // above (bearing ~ pi - eps). The true angular error is tiny.
        let truth = State::new(-5000.0, -1.0, 0.0, 0.0);
        let mut ekf = Ekf::new(truth, Matrix4::identity() * 100.0, 1.0);
        let mut z = radar().h(&truth);
        z.bearing = wrap_angle(z.bearing + 0.002);
        let nis = ekf.update_radar(&z, &radar()).expect("update applied");
        assert!(nis < 10.0, "branch cut not handled, NIS = {nis}");
    }
}
