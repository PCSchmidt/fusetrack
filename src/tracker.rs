//! Global-nearest-neighbour (GNN) multi-target tracker.
//!
//! Per scan: predict every track, gate measurements with a chi-square test on
//! the innovation, solve the track-to-measurement assignment optimally with
//! the Hungarian algorithm, update assigned tracks, and start tentative
//! tracks from unassigned measurements. Tentative tracks are confirmed on
//! M hits within N scans and deleted otherwise, or alternatively by a
//! sequential likelihood-ratio score (see [`ScoreConfig`]); confirmed tracks
//! are deleted after too many consecutive misses.

use nalgebra::{Matrix4, Vector2};

use crate::assignment::{self, FORBIDDEN};
use crate::ekf::{Ekf, Innovation};
use crate::models::{Radar, RadarMeas, State};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TrackStatus {
    Tentative,
    Confirmed,
}

#[derive(Debug, Clone)]
pub struct Track {
    pub id: u64,
    pub ekf: Ekf,
    pub status: TrackStatus,
    /// Scans since the track was born (including the birth scan).
    pub age: usize,
    /// Scans on which a measurement was assigned (including birth).
    pub hits: usize,
    /// Consecutive scans without an assigned measurement.
    pub misses: usize,
    /// Log-likelihood-ratio score (track vs. clutter); only used when
    /// [`TrackerConfig::score`] is set.
    pub score: f64,
}

impl Track {
    pub fn position(&self) -> Vector2<f64> {
        Vector2::new(self.ekf.x[0], self.ekf.x[1])
    }
}

/// Sequential likelihood-ratio track scoring (Blackman).
///
/// Each scan adds `ln(Pd * N(y; 0, S) / (rho * range))` on a hit and
/// `ln(1 - Pd)` on a miss, where `rho` is the clutter density per square
/// metre and `rho * range` converts it to density per metre-radian of
/// range-bearing space. This assumes `pd` and `clutter_density` are known.
#[derive(Debug, Clone, Copy)]
pub struct ScoreConfig {
    pub pd: f64,
    /// Expected clutter returns per square metre per scan.
    pub clutter_density: f64,
    /// Confirm a tentative track when its score reaches this value.
    pub confirm: f64,
    /// Delete a tentative track when its score falls below this value.
    pub delete: f64,
}

#[derive(Debug, Clone)]
pub struct TrackerConfig {
    /// Chi-square gate on the innovation (2 DOF: 13.8 is the 99.9% point).
    pub gate: f64,
    /// Process-noise acceleration standard deviation (m/s^2).
    pub sigma_accel: f64,
    /// Velocity standard deviation assumed when initialising a track (m/s).
    pub init_speed_sigma: f64,
    /// Confirm a tentative track after `confirm_hits` hits in its first
    /// `confirm_window` scans.
    pub confirm_hits: usize,
    pub confirm_window: usize,
    /// Delete a confirmed track after this many consecutive misses.
    pub max_misses: usize,
    /// If set, tentative tracks are confirmed and deleted by score instead of
    /// the M-of-N rule.
    pub score: Option<ScoreConfig>,
}

impl Default for TrackerConfig {
    fn default() -> Self {
        Self {
            gate: 13.8,
            sigma_accel: 2.0,
            init_speed_sigma: 30.0,
            confirm_hits: 3,
            confirm_window: 5,
            max_misses: 5,
            score: None,
        }
    }
}

#[derive(Debug, Clone)]
pub struct GnnTracker {
    pub cfg: TrackerConfig,
    pub tracks: Vec<Track>,
    next_id: u64,
}

impl GnnTracker {
    pub fn new(cfg: TrackerConfig) -> Self {
        Self {
            cfg,
            tracks: Vec::new(),
            next_id: 0,
        }
    }

    /// Process one scan of measurements taken `dt` seconds after the last.
    pub fn step(&mut self, dt: f64, meas: &[RadarMeas], radar: &Radar) {
        for t in &mut self.tracks {
            t.ekf.predict(dt);
        }

        // Gated cost matrix: NIS + ln|S| (negative log-likelihood up to a
        // constant), forbidden outside the gate.
        let mut innovations: Vec<Vec<Option<Innovation>>> = Vec::with_capacity(self.tracks.len());
        let mut cost: Vec<Vec<f64>> = Vec::with_capacity(self.tracks.len());
        for t in &self.tracks {
            let mut inn_row = Vec::with_capacity(meas.len());
            let mut cost_row = Vec::with_capacity(meas.len());
            let model = t.ekf.measurement_model(radar);
            let ln_det_s = model.map_or(0.0, |m| m.s.determinant().ln());
            for z in meas {
                match model.map(|m| m.innovation(z)) {
                    Some(inn) if inn.nis() <= self.cfg.gate => {
                        cost_row.push(inn.nis() + ln_det_s);
                        inn_row.push(Some(inn));
                    }
                    _ => {
                        cost_row.push(FORBIDDEN);
                        inn_row.push(None);
                    }
                }
            }
            innovations.push(inn_row);
            cost.push(cost_row);
        }

        let assigned = assignment::solve(&cost);
        let mut meas_used = vec![false; meas.len()];
        let score_cfg = self.cfg.score;
        for (i, t) in self.tracks.iter_mut().enumerate() {
            t.age += 1;
            match assigned[i] {
                Some(j) => {
                    let inn = innovations[i][j].expect("assigned pairs are gated in");
                    if let Some(sc) = score_cfg {
                        let ln_n = -0.5 * inn.nis()
                            - std::f64::consts::TAU.ln()
                            - 0.5 * inn.s.determinant().ln();
                        let density = sc.clutter_density * meas[j].range.max(1.0);
                        t.score += sc.pd.min(0.999).ln() + ln_n - density.ln();
                    }
                    t.ekf.update_with(&inn, radar);
                    t.hits += 1;
                    t.misses = 0;
                    meas_used[j] = true;
                }
                None => {
                    if let Some(sc) = score_cfg {
                        t.score += (1.0 - sc.pd.min(0.999)).ln();
                    }
                    t.misses += 1;
                }
            }
        }

        // Status changes and deletions.
        let cfg = &self.cfg;
        for t in &mut self.tracks {
            let confirmed = match cfg.score {
                Some(sc) => t.score >= sc.confirm,
                None => t.hits >= cfg.confirm_hits,
            };
            if t.status == TrackStatus::Tentative && confirmed {
                t.status = TrackStatus::Confirmed;
            }
        }
        self.tracks.retain(|t| match (t.status, cfg.score) {
            (TrackStatus::Tentative, Some(sc)) => t.score >= sc.delete,
            (TrackStatus::Tentative, None) => {
                t.age < cfg.confirm_window || t.hits >= cfg.confirm_hits
            }
            (TrackStatus::Confirmed, _) => t.misses <= cfg.max_misses,
        });

        // New tentative tracks from unused measurements.
        for (j, z) in meas.iter().enumerate() {
            if !meas_used[j] {
                let track = self.spawn(z, radar);
                self.tracks.push(track);
            }
        }
    }

    fn spawn(&mut self, z: &RadarMeas, radar: &Radar) -> Track {
        let (s, c) = z.bearing.sin_cos();
        let x = State::new(
            radar.pos[0] + z.range * c,
            radar.pos[1] + z.range * s,
            0.0,
            0.0,
        );
        // Isotropic position uncertainty: the larger of range and cross-range
        // sigma, a simple conservative choice.
        let sp = radar.sigma_r.max(z.range * radar.sigma_b);
        let sv = self.cfg.init_speed_sigma;
        let p = Matrix4::from_diagonal(&nalgebra::Vector4::new(sp * sp, sp * sp, sv * sv, sv * sv));
        let id = self.next_id;
        self.next_id += 1;
        Track {
            id,
            ekf: Ekf::new(x, p, self.cfg.sigma_accel),
            status: TrackStatus::Tentative,
            age: 1,
            hits: 1,
            misses: 0,
            score: 0.0,
        }
    }

    /// Positions of confirmed tracks.
    pub fn confirmed_positions(&self) -> Vec<Vector2<f64>> {
        self.tracks
            .iter()
            .filter(|t| t.status == TrackStatus::Confirmed)
            .map(|t| t.position())
            .collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn radar() -> Radar {
        Radar {
            pos: Vector2::new(0.0, 0.0),
            sigma_r: 15.0,
            sigma_b: 0.003,
        }
    }

    #[test]
    fn a_noiseless_target_is_confirmed_after_three_scans() {
        let radar = radar();
        let mut tracker = GnnTracker::new(TrackerConfig::default());
        for k in 0..3 {
            let truth = State::new(5000.0 + 20.0 * k as f64, 2000.0, 20.0, 0.0);
            tracker.step(1.0, &[radar.h(&truth)], &radar);
        }
        assert_eq!(tracker.confirmed_positions().len(), 1);
        assert_eq!(tracker.tracks.len(), 1);
    }

    #[test]
    fn isolated_clutter_never_confirms_and_is_deleted() {
        let radar = radar();
        let mut tracker = GnnTracker::new(TrackerConfig::default());
        // One return on scan 0 at a place no later return is near.
        tracker.step(
            1.0,
            &[RadarMeas {
                range: 4000.0,
                bearing: 1.0,
            }],
            &radar,
        );
        assert_eq!(tracker.tracks.len(), 1);
        for _ in 0..6 {
            tracker.step(1.0, &[], &radar);
        }
        assert!(tracker.tracks.is_empty());
        assert!(tracker.confirmed_positions().is_empty());
    }

    #[test]
    fn confirmed_track_coasts_through_missed_detections_then_dies() {
        let radar = radar();
        let mut tracker = GnnTracker::new(TrackerConfig::default());
        for k in 0..4 {
            let truth = State::new(5000.0 + 20.0 * k as f64, 2000.0, 20.0, 0.0);
            tracker.step(1.0, &[radar.h(&truth)], &radar);
        }
        assert_eq!(tracker.confirmed_positions().len(), 1);
        for _ in 0..tracker.cfg.max_misses {
            tracker.step(1.0, &[], &radar);
        }
        assert_eq!(tracker.confirmed_positions().len(), 1, "still coasting");
        tracker.step(1.0, &[], &radar);
        assert!(tracker.tracks.is_empty(), "deleted after max_misses + 1");
    }

    #[test]
    fn two_separated_targets_get_two_tracks_with_stable_ids() {
        let radar = radar();
        let mut tracker = GnnTracker::new(TrackerConfig::default());
        let mut ids_at_3: Vec<u64> = Vec::new();
        for k in 0..15 {
            let a = State::new(5000.0 + 20.0 * k as f64, 2000.0, 20.0, 0.0);
            let b = State::new(-4000.0, -3000.0 + 15.0 * k as f64, 0.0, 15.0);
            tracker.step(1.0, &[radar.h(&a), radar.h(&b)], &radar);
            if k == 3 {
                ids_at_3 = tracker.tracks.iter().map(|t| t.id).collect();
            }
        }
        let mut ids_end: Vec<u64> = tracker.tracks.iter().map(|t| t.id).collect();
        ids_end.sort_unstable();
        ids_at_3.sort_unstable();
        assert_eq!(ids_end.len(), 2);
        assert_eq!(ids_end, ids_at_3, "no ID churn on clean data");
    }

    fn score_cfg() -> TrackerConfig {
        TrackerConfig {
            score: Some(ScoreConfig {
                pd: 0.9,
                clutter_density: 1.0e-8,
                confirm: 12.0,
                delete: -5.0,
            }),
            ..TrackerConfig::default()
        }
    }

    #[test]
    fn score_confirms_a_consistent_target_within_a_few_scans() {
        let radar = radar();
        let mut tracker = GnnTracker::new(score_cfg());
        for k in 0..4 {
            let truth = State::new(5000.0 + 20.0 * k as f64, 2000.0, 20.0, 0.0);
            tracker.step(1.0, &[radar.h(&truth)], &radar);
        }
        assert_eq!(tracker.confirmed_positions().len(), 1);
    }

    #[test]
    fn score_deletes_a_lone_return_after_a_few_misses() {
        let radar = radar();
        let mut tracker = GnnTracker::new(score_cfg());
        tracker.step(
            1.0,
            &[RadarMeas {
                range: 4000.0,
                bearing: 1.0,
            }],
            &radar,
        );
        assert_eq!(tracker.tracks.len(), 1);
        for _ in 0..3 {
            tracker.step(1.0, &[], &radar);
        }
        assert!(tracker.tracks.is_empty());
    }
}
