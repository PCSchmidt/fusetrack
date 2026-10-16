//! Extended Kalman filter with a constant-velocity motion model and a
//! range-bearing radar update.

use nalgebra::{Matrix2, Matrix2x4, Matrix4, Matrix4x2, Vector2};

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

    /// Linearised measurement model at the current estimate. The innovation
    /// covariance depends on the track but not on the measurement, so it can
    /// be computed once per track and reused for every candidate measurement.
    /// Returns `None` if the target is on top of the sensor or `S` is singular.
    pub fn measurement_model(&self, radar: &Radar) -> Option<MeasurementModel> {
        let dx = self.x[0] - radar.pos[0];
        let dy = self.x[1] - radar.pos[1];
        let r2 = dx * dx + dy * dy;
        if r2 < 1e-12 {
            return None;
        }
        let r = r2.sqrt();
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
        Some(MeasurementModel {
            pred: radar.h(&self.x),
            h,
            s,
            s_inv,
        })
    }

    /// Innovation for a candidate measurement, without changing the state.
    pub fn innovation(&self, z: &RadarMeas, radar: &Radar) -> Option<Innovation> {
        Some(self.measurement_model(radar)?.innovation(z))
    }

    /// Apply a previously computed innovation. Returns the NIS.
    pub fn update_with(&mut self, inn: &Innovation, radar: &Radar) -> f64 {
        let rm = Matrix2::new(
            radar.sigma_r.powi(2),
            0.0, //
            0.0,
            radar.sigma_b.powi(2),
        );
        let k: Matrix4x2<f64> = self.p * inn.h.transpose() * inn.s_inv;
        self.x += k * inn.y;
        // Joseph form keeps P symmetric positive semi-definite.
        let i_kh = Matrix4::identity() - k * inn.h;
        self.p = i_kh * self.p * i_kh.transpose() + k * rm * k.transpose();
        inn.nis()
    }

    /// Fuse a radar measurement. Returns the normalised innovation squared
    /// (NIS), or `None` if the update was skipped (target on top of the
    /// sensor, or singular innovation covariance).
    pub fn update_radar(&mut self, z: &RadarMeas, radar: &Radar) -> Option<f64> {
        let inn = self.innovation(z, radar)?;
        Some(self.update_with(&inn, radar))
    }
}

/// Per-track linearised measurement model; see [`Ekf::measurement_model`].
#[derive(Debug, Clone, Copy)]
pub struct MeasurementModel {
    pub pred: RadarMeas,
    pub h: Matrix2x4<f64>,
    pub s: Matrix2<f64>,
    pub s_inv: Matrix2<f64>,
}

impl MeasurementModel {
    pub fn innovation(&self, z: &RadarMeas) -> Innovation {
        Innovation {
            y: Vector2::new(
                z.range - self.pred.range,
                wrap_angle(z.bearing - self.pred.bearing),
            ),
            h: self.h,
            s: self.s,
            s_inv: self.s_inv,
        }
    }
}

/// Result of comparing a measurement with a filter's prediction.
#[derive(Debug, Clone, Copy)]
pub struct Innovation {
    pub y: Vector2<f64>,
    pub h: Matrix2x4<f64>,
    pub s: Matrix2<f64>,
    pub s_inv: Matrix2<f64>,
}

impl Innovation {
    /// Normalised innovation squared (squared Mahalanobis distance).
    pub fn nis(&self) -> f64 {
        (self.y.transpose() * self.s_inv * self.y)[0]
    }

    /// Natural log of the determinant of the innovation covariance.
    pub fn ln_det_s(&self) -> f64 {
        self.s.determinant().ln()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

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
