//! Global-nearest-neighbour (GNN) multi-target tracker.
//!
//! Per scan: predict every track, gate measurements with a chi-square test on
//! the innovation, solve the track-to-measurement assignment optimally with
//! the Hungarian algorithm, update assigned tracks, and start tentative
//! tracks from unassigned measurements. Tentative tracks are confirmed on
//! M hits within N scans and deleted otherwise, or alternatively by a
//! sequential likelihood-ratio score (see [`ScoreConfig`]); confirmed tracks
//! are deleted after too many consecutive misses.
//!
//! Each track runs either a single constant-velocity EKF or, if
//! [`TrackerConfig::imm`] is set, an IMM filter. IMM tracks are gated and
//! costed with the combined predicted state and covariance, then updated by
//! the IMM so each model is reweighted by its own likelihood.

use nalgebra::{Matrix4, Vector2};

use crate::assignment::{self, FORBIDDEN};
use crate::ekf::{Ekf, Innovation, MeasurementModel};
use crate::imm::{Imm, ImmConfig};
use crate::jpda::{self, JpdaConfig};
use crate::models::{Radar, RadarMeas, State};

/// Two tracks whose states differ by less than this (chi-square, 4 degrees of freedom, 95%) given their
/// combined covariance are treated as the same target by JPDA's duplicate merge.
const MERGE_CHI2: f64 = 9.49;

/// Lowest clutter density (returns per square metre per scan) JPDA will use.
const JPDA_MIN_DENSITY: f64 = 1.0e-15;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TrackStatus {
    Tentative,
    Confirmed,
}

/// The filter carried by a track.
#[derive(Debug, Clone)]
pub enum TrackFilter {
    /// Single constant-velocity EKF.
    Single(Ekf),
    /// Interacting multiple model filter.
    Multi(Imm),
}

impl TrackFilter {
    fn predict(&mut self, dt: f64) {
        match self {
            TrackFilter::Single(e) => e.predict(dt),
            TrackFilter::Multi(i) => i.predict(dt),
        }
    }

    /// Linearised measurement model at the (combined) predicted state.
    fn measurement_model(&self, radar: &Radar) -> Option<MeasurementModel> {
        match self {
            TrackFilter::Single(e) => e.measurement_model(radar),
            TrackFilter::Multi(i) => {
                let (x, p) = i.estimate();
                Ekf::new(x, p, 0.0).measurement_model(radar)
            }
        }
    }

    /// `inn` is the gating innovation of the combined prediction; the IMM
    /// recomputes a per-model innovation internally.
    fn update(&mut self, inn: &Innovation, z: &RadarMeas, radar: &Radar) {
        match self {
            TrackFilter::Single(e) => {
                e.update_with(inn, radar);
            }
            TrackFilter::Multi(i) => {
                i.update(z, radar);
            }
        }
    }

    /// Current state estimate.
    pub fn state(&self) -> State {
        match self {
            TrackFilter::Single(e) => e.x,
            TrackFilter::Multi(i) => i.estimate().0,
        }
    }
}

#[derive(Debug, Clone)]
pub struct Track {
    pub id: u64,
    pub filter: TrackFilter,
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
        let x = self.filter.state();
        Vector2::new(x[0], x[1])
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
    /// If set, tracks run this IMM instead of a single constant-velocity EKF
    /// (`sigma_accel` is then ignored). The first model should be constant
    /// velocity: new tracks start with probability 0.8 on it and the rest
    /// shared equally across the other models.
    pub imm: Option<ImmConfig>,
    /// If set, measurements are associated by JPDA instead of GNN assignment: every track is updated with all
    /// measurements in its gate, weighted by association probability. Single-EKF tracks only (`imm` must be
    /// `None`). Track confirmation then uses the PDA likelihood ratio when `score` is set, or hit counting
    /// (a hit is a scan on which the track's miss probability is at most one half) otherwise.
    pub jpda: Option<JpdaConfig>,
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
            imm: None,
            jpda: None,
        }
    }
}

#[derive(Debug, Clone)]
pub struct GnnTracker {
    pub cfg: TrackerConfig,
    pub tracks: Vec<Track>,
    next_id: u64,
    /// Scans on which JPDA had a cluster too large to enumerate and used the independent-PDA fallback.
    pub jpda_inexact_scans: usize,
}

impl GnnTracker {
    pub fn new(cfg: TrackerConfig) -> Self {
        assert!(
            cfg.jpda.is_none() || cfg.imm.is_none(),
            "JPDA is implemented for single-EKF tracks; combine it with the IMM in a later milestone"
        );
        Self {
            cfg,
            tracks: Vec::new(),
            next_id: 0,
            jpda_inexact_scans: 0,
        }
    }

    /// Process one scan of measurements taken `dt` seconds after the last.
    pub fn step(&mut self, dt: f64, meas: &[RadarMeas], radar: &Radar) {
        for t in &mut self.tracks {
            t.filter.predict(dt);
        }

        // Gated cost matrix: NIS + ln|S| (negative log-likelihood up to a
        // constant), forbidden outside the gate.
        let mut innovations: Vec<Vec<Option<Innovation>>> = Vec::with_capacity(self.tracks.len());
        let mut cost: Vec<Vec<f64>> = Vec::with_capacity(self.tracks.len());
        for t in &self.tracks {
            let mut inn_row = Vec::with_capacity(meas.len());
            let mut cost_row = Vec::with_capacity(meas.len());
            let model = t.filter.measurement_model(radar);
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

        let meas_used: Vec<bool> = if let Some(jp) = self.cfg.jpda {
            self.associate_jpda(jp, &innovations, meas, radar)
        } else {
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
                        t.filter.update(&inn, &meas[j], radar);
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
            meas_used
        };

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

        // JPDA lets two tracks share one target's measurements instead of starving the second, so duplicates
        // survive; drop the weaker of any pair that is statistically the same target.
        if self.cfg.jpda.is_some_and(|jp| jp.merge_duplicates) {
            self.merge_duplicates();
        }

        // New tentative tracks from unused measurements.
        for (j, z) in meas.iter().enumerate() {
            if !meas_used[j] {
                let track = self.spawn(z, radar);
                self.tracks.push(track);
            }
        }
    }

    /// JPDA association and update for one scan; returns which measurements are explained by existing tracks.
    fn associate_jpda(
        &mut self,
        jp: JpdaConfig,
        innovations: &[Vec<Option<Innovation>>],
        meas: &[RadarMeas],
        radar: &Radar,
    ) -> Vec<bool> {
        // detection weight Pd * N(y; 0, S) / lambda for every gated pair
        let weights: Vec<Vec<Option<f64>>> = innovations
            .iter()
            .map(|row| {
                row.iter()
                    .enumerate()
                    .map(|(j, inn)| {
                        inn.as_ref().map(|inn| {
                            let n_ij = (-0.5 * inn.nis()).exp()
                                / (std::f64::consts::TAU * inn.s.determinant().sqrt());
                            // a zero density would make every weight infinite; floor it (no clutter is then
                            // modelled as an extremely low density)
                            let lambda =
                                jp.clutter_density.max(JPDA_MIN_DENSITY) * meas[j].range.max(1.0);
                            jp.pd * n_ij / lambda
                        })
                    })
                    .collect()
            })
            .collect();
        let marg = jpda::marginals(&weights, 1.0 - jp.pd, jp.max_events);
        if !marg.exact {
            self.jpda_inexact_scans += 1;
        }

        let scoring = self.cfg.score.is_some();
        let mut claimed = vec![0.0; meas.len()];
        for (i, t) in self.tracks.iter_mut().enumerate() {
            t.age += 1;
            let mut pairs: Vec<(f64, Innovation)> = Vec::new();
            for (j, inn) in innovations[i].iter().enumerate() {
                if let Some(inn) = inn {
                    pairs.push((marg.beta[i][j], *inn));
                    claimed[j] += marg.beta[i][j];
                }
            }
            if !pairs.is_empty() {
                let TrackFilter::Single(ekf) = &mut t.filter else {
                    unreachable!("JPDA requires single-EKF tracks")
                };
                jpda::pda_update(ekf, radar, &pairs, marg.miss[i]);
            }
            if scoring {
                // PDA track likelihood ratio: a miss plus every gated measurement as a possible detection
                let lr = (1.0 - jp.pd) + weights[i].iter().flatten().sum::<f64>();
                t.score += lr.ln();
            }
            if !pairs.is_empty() && marg.miss[i] <= 0.5 {
                t.hits += 1;
                t.misses = 0;
            } else {
                t.misses += 1;
            }
        }
        claimed.iter().map(|c| *c >= jp.birth_threshold).collect()
    }

    /// Drop the weaker of each pair of tracks whose states are statistically indistinguishable. Only pairs
    /// with at least one confirmed track are considered, so two young tentative tracks born close together
    /// are left alone. The stronger track is the confirmed one, then the one with more hits, then the older.
    /// This also merges two *real* targets that move together, which is a known cost of the rule.
    fn merge_duplicates(&mut self) {
        let n = self.tracks.len();
        let mut drop = vec![false; n];
        let rank = |t: &Track| {
            (
                t.status == TrackStatus::Confirmed,
                t.hits,
                std::cmp::Reverse(t.id),
            )
        };
        for i in 0..n {
            for j in i + 1..n {
                if drop[i] || drop[j] {
                    continue;
                }
                let (ti, tj) = (&self.tracks[i], &self.tracks[j]);
                if ti.status == TrackStatus::Tentative && tj.status == TrackStatus::Tentative {
                    continue;
                }
                let (TrackFilter::Single(a), TrackFilter::Single(b)) = (&ti.filter, &tj.filter)
                else {
                    continue;
                };
                let d = a.x - b.x;
                if let Some(inv) = (a.p + b.p).try_inverse() {
                    if (d.transpose() * inv * d)[0] < MERGE_CHI2 {
                        if rank(ti) >= rank(tj) {
                            drop[j] = true;
                        } else {
                            drop[i] = true;
                        }
                    }
                }
            }
        }
        let mut k = 0;
        self.tracks.retain(|_| {
            k += 1;
            !drop[k - 1]
        });
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
            filter: match &self.cfg.imm {
                None => TrackFilter::Single(Ekf::new(x, p, self.cfg.sigma_accel)),
                Some(imm_cfg) => {
                    let r = imm_cfg.models.len();
                    let mu0 = if r == 1 {
                        vec![1.0]
                    } else {
                        let mut v = vec![0.2 / (r - 1) as f64; r];
                        v[0] = 0.8;
                        v
                    };
                    TrackFilter::Multi(Imm::new(imm_cfg.clone(), x, p, mu0))
                }
            },
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

    #[test]
    fn imm_tracks_confirm_and_follow_a_turning_target() {
        use crate::imm::ImmConfig;
        let radar = radar();
        let cfg = TrackerConfig {
            imm: Some(ImmConfig::cv_ct(0.1, 0.5, 0.5, 0.9)),
            ..TrackerConfig::default()
        };
        let mut tracker = GnnTracker::new(cfg);
        let mut truth = State::new(5000.0, 2000.0, -25.0, 15.0);
        for k in 0..40 {
            // Straight for 15 scans, then a constant turn.
            let omega = if k < 15 { 0.0 } else { 0.1 };
            truth = crate::ekf::Ekf::ct_transition(omega, 1.0) * truth;
            tracker.step(1.0, &[radar.h(&truth)], &radar);
        }
        assert_eq!(tracker.confirmed_positions().len(), 1);
        assert_eq!(tracker.tracks.len(), 1, "no track lost or duplicated");
        let est = tracker.tracks[0].filter.state();
        assert!((est[0] - truth[0]).hypot(est[1] - truth[1]) < 30.0);
        assert!(matches!(tracker.tracks[0].filter, TrackFilter::Multi(_)));
    }
}
