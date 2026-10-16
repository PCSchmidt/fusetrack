//! Interacting Multiple Model (IMM) filter (Blom & Bar-Shalom, 1988).
//!
//! Runs several motion-model EKFs in parallel on the same target and blends
//! them by model probability, so one filter can cover straight flight and
//! turns. See `docs/imm-derivation.md` for the mixing and combination steps.

use std::f64::consts::TAU;

use nalgebra::Matrix4;

use crate::ekf::Ekf;
use crate::models::{Radar, RadarMeas, State};

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum MotionModel {
    ConstantVelocity,
    /// Coordinated turn at a fixed turn rate (rad/s, counter-clockwise > 0).
    CoordinatedTurn {
        omega: f64,
    },
}

impl MotionModel {
    pub fn transition(&self, dt: f64) -> Matrix4<f64> {
        match *self {
            MotionModel::ConstantVelocity => Ekf::transition(dt),
            MotionModel::CoordinatedTurn { omega } => Ekf::ct_transition(omega, dt),
        }
    }
}

#[derive(Debug, Clone)]
pub struct ImmConfig {
    /// Each model with its process-noise acceleration sigma (m/s^2).
    pub models: Vec<(MotionModel, f64)>,
    /// Markov model-switching matrix: `pi[i][j] = P(model j at k | model i at k-1)`.
    /// Rows must sum to one and every entry must be positive.
    pub pi: Vec<Vec<f64>>,
}

impl ImmConfig {
    /// Constant velocity plus coordinated turns at `+omega` and `-omega`,
    /// with the given probability of staying in the same model each scan.
    pub fn cv_ct(omega: f64, sigma_cv: f64, sigma_ct: f64, stay: f64) -> Self {
        let off = (1.0 - stay) / 2.0;
        Self {
            models: vec![
                (MotionModel::ConstantVelocity, sigma_cv),
                (MotionModel::CoordinatedTurn { omega }, sigma_ct),
                (MotionModel::CoordinatedTurn { omega: -omega }, sigma_ct),
            ],
            pi: vec![
                vec![stay, off, off],
                vec![off, stay, off],
                vec![off, off, stay],
            ],
        }
    }
}

#[derive(Debug, Clone)]
pub struct Imm {
    cfg: ImmConfig,
    pub filters: Vec<Ekf>,
    /// Current model probabilities.
    pub mu: Vec<f64>,
}

impl Imm {
    /// All models start from the same state and covariance.
    pub fn new(cfg: ImmConfig, x: State, p: Matrix4<f64>, mu0: Vec<f64>) -> Self {
        let r = cfg.models.len();
        assert_eq!(mu0.len(), r, "one initial probability per model");
        assert_eq!(cfg.pi.len(), r, "pi must be r x r");
        assert!(cfg.pi.iter().all(|row| row.len() == r), "pi must be r x r");
        let filters = cfg
            .models
            .iter()
            .map(|&(_, sigma)| Ekf::new(x, p, sigma))
            .collect();
        Self {
            cfg,
            filters,
            mu: mu0,
        }
    }

    /// Mixing followed by per-model prediction.
    pub fn predict(&mut self, dt: f64) {
        let r = self.filters.len();
        // Predicted model probabilities c_j = sum_i pi_ij mu_i.
        let c: Vec<f64> = (0..r)
            .map(|j| (0..r).map(|i| self.cfg.pi[i][j] * self.mu[i]).sum())
            .collect();
        let prior: Vec<(State, Matrix4<f64>)> = self.filters.iter().map(|f| (f.x, f.p)).collect();

        for (j, &c_j) in c.iter().enumerate() {
            // Mixing weights mu_{i|j}.
            let w: Vec<f64> = (0..r)
                .map(|i| self.cfg.pi[i][j] * self.mu[i] / c_j)
                .collect();
            let x0: State = (0..r).map(|i| w[i] * prior[i].0).sum();
            let p0: Matrix4<f64> = (0..r)
                .map(|i| {
                    let d = prior[i].0 - x0;
                    w[i] * (prior[i].1 + d * d.transpose())
                })
                .sum();

            let (model, sigma) = self.cfg.models[j];
            let f = model.transition(dt);
            let q = Ekf::process_noise(dt, sigma);
            self.filters[j].x = x0;
            self.filters[j].p = p0;
            self.filters[j].predict_with(&f, &q);
        }
        self.mu = c;
    }

    /// Update every model with a measurement and reweight by likelihood.
    /// Returns `false` (and changes nothing) if no model could be updated.
    pub fn update(&mut self, z: &RadarMeas, radar: &Radar) -> bool {
        let r = self.filters.len();
        let inns: Vec<_> = self
            .filters
            .iter()
            .map(|f| f.innovation(z, radar))
            .collect();
        if inns.iter().all(|i| i.is_none()) {
            return false;
        }

        // Log-likelihoods of the innovation under each model's N(0, S).
        let log_w: Vec<f64> = (0..r)
            .map(|j| match &inns[j] {
                Some(inn) => {
                    self.mu[j].ln() - 0.5 * inn.nis() - TAU.ln() - 0.5 * inn.s.determinant().ln()
                }
                None => f64::NEG_INFINITY,
            })
            .collect();
        let max = log_w.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        let w: Vec<f64> = log_w.iter().map(|l| (l - max).exp()).collect();
        let total: f64 = w.iter().sum();

        for (j, inn) in inns.iter().enumerate() {
            if let Some(inn) = inn {
                self.filters[j].update_with(inn, radar);
            }
        }
        self.mu = w.iter().map(|x| x / total).collect();
        true
    }

    /// Moment-matched combined estimate over all models.
    pub fn estimate(&self) -> (State, Matrix4<f64>) {
        let x: State = self
            .filters
            .iter()
            .zip(&self.mu)
            .map(|(f, m)| *m * f.x)
            .sum();
        let p: Matrix4<f64> = self
            .filters
            .iter()
            .zip(&self.mu)
            .map(|(f, m)| {
                let d = f.x - x;
                *m * (f.p + d * d.transpose())
            })
            .sum();
        (x, p)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::{Vector2, Vector4};

    fn radar() -> Radar {
        Radar {
            pos: Vector2::new(0.0, 0.0),
            sigma_r: 15.0,
            sigma_b: 0.003,
        }
    }

    fn p0() -> Matrix4<f64> {
        Matrix4::from_diagonal(&Vector4::new(100.0, 100.0, 25.0, 25.0))
    }

    #[test]
    fn a_single_model_imm_equals_the_plain_ekf() {
        let cfg = ImmConfig {
            models: vec![(MotionModel::ConstantVelocity, 1.0)],
            pi: vec![vec![1.0]],
        };
        let x0 = State::new(4000.0, 2000.0, -20.0, 5.0);
        let mut imm = Imm::new(cfg, x0, p0(), vec![1.0]);
        let mut ekf = Ekf::new(x0, p0(), 1.0);
        let radar = radar();
        for k in 1..=20 {
            let truth = State::new(
                4000.0 - 20.0 * k as f64,
                2000.0 + 5.0 * k as f64,
                -20.0,
                5.0,
            );
            let z = radar.h(&truth);
            imm.predict(1.0);
            ekf.predict(1.0);
            assert!(imm.update(&z, &radar));
            ekf.update_radar(&z, &radar);
        }
        let (x, p) = imm.estimate();
        assert!((x - ekf.x).norm() < 1e-9);
        assert!((p - ekf.p).norm() < 1e-9);
        assert!((imm.mu[0] - 1.0).abs() < 1e-12);
    }

    #[test]
    fn mixing_identical_models_changes_nothing() {
        let cfg = ImmConfig::cv_ct(0.1, 0.5, 0.5, 0.9);
        let x0 = State::new(4000.0, 2000.0, -20.0, 5.0);
        let mut imm = Imm::new(cfg, x0, p0(), vec![0.5, 0.25, 0.25]);
        // Predicting with dt = 0 isolates the mixing step.
        imm.predict(0.0);
        for f in &imm.filters {
            assert!((f.x - x0).norm() < 1e-9);
        }
        // Probabilities remain a distribution.
        assert!((imm.mu.iter().sum::<f64>() - 1.0).abs() < 1e-12);
    }

    #[test]
    fn model_probabilities_always_sum_to_one() {
        let cfg = ImmConfig::cv_ct(0.1, 0.5, 0.5, 0.9);
        let x0 = State::new(4000.0, 2000.0, -20.0, 5.0);
        let mut imm = Imm::new(cfg, x0, p0(), vec![0.8, 0.1, 0.1]);
        let radar = radar();
        for k in 1..=30 {
            let truth = State::new(
                4000.0 - 20.0 * k as f64,
                2000.0 + 5.0 * k as f64,
                -20.0,
                5.0,
            );
            imm.predict(1.0);
            imm.update(&radar.h(&truth), &radar);
            assert!((imm.mu.iter().sum::<f64>() - 1.0).abs() < 1e-9);
            assert!(imm.mu.iter().all(|m| *m >= 0.0));
        }
    }

    #[test]
    fn straight_flight_favours_constant_velocity_and_a_turn_favours_the_matching_ct() {
        let cfg = ImmConfig::cv_ct(0.1, 0.5, 0.5, 0.9);
        let radar = radar();
        let x0 = State::new(5000.0, 2000.0, -30.0, 20.0);
        let mut imm = Imm::new(cfg, x0, p0(), vec![0.34, 0.33, 0.33]);

        // Noise-free straight flight.
        let mut truth = x0;
        for _ in 0..25 {
            truth = Ekf::transition(1.0) * truth;
            imm.predict(1.0);
            imm.update(&radar.h(&truth), &radar);
        }
        assert!(imm.mu[0] > 0.7, "CV should dominate: {:?}", imm.mu);

        // Noise-free counter-clockwise turn at the model's rate.
        for _ in 0..25 {
            truth = Ekf::ct_transition(0.1, 1.0) * truth;
            imm.predict(1.0);
            imm.update(&radar.h(&truth), &radar);
        }
        assert!(
            imm.mu[1] > imm.mu[0] && imm.mu[1] > imm.mu[2],
            "{:?}",
            imm.mu
        );
    }
}
