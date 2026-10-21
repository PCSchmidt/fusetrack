//! Joint probabilistic data association (JPDA).
//!
//! Instead of forcing each track onto one measurement (GNN), JPDA updates every track with *all*
//! measurements in its gate, weighted by the probability that each one is the track's own. Those
//! weights come from the joint events consistent with the rule that a measurement can come from at most
//! one track and a track produces at most one measurement.
//!
//! # Weights
//!
//! For track `i` and gated measurement `j` the detection weight is `w_ij = Pd * N(y_ij; 0, S_i) / lambda_j`,
//! where `lambda_j` is the clutter density at the measurement (per metre-radian of range-bearing space) and
//! a missed detection has weight `1 - Pd`. A joint event assigns each track in a cluster either a distinct
//! measurement or a miss; its probability is proportional to the product of those weights (measurements
//! left unassigned are clutter, which is what dividing by `lambda_j` accounts for). Marginals `beta_ij`
//! are sums over events. Tracks that share no gated measurement are independent, so events are enumerated
//! per cluster of connected tracks. The number of events grows quickly with cluster size; past
//! `max_events` the cluster falls back to independent per-track PDA (exclusivity ignored) and the result
//! is flagged inexact.
//!
//! # Update
//!
//! The PDA update of Bar-Shalom and Fortmann, with the combined innovation and the spread-of-the-means
//! term in the covariance:
//! `P = b0 P- + (1 - b0) Pc + K (sum_j b_j v_j v_j' - v v') K'`, `v = sum_j b_j v_j`, where `Pc` is the
//! covariance after an update with a correct measurement.

use nalgebra::{Matrix2, Matrix4, Matrix4x2, Vector2};

use crate::ekf::{Ekf, Innovation};
use crate::models::Radar;

/// JPDA settings. `pd` and `clutter_density` are assumed known.
#[derive(Debug, Clone, Copy)]
pub struct JpdaConfig {
    pub pd: f64,
    /// Expected clutter returns per square metre per scan (as in [`crate::ScoreConfig`]); converted to
    /// range-bearing space by multiplying by the measurement's range. Values below 1e-15 are raised to it,
    /// so a clutter-free scenario does not divide by zero.
    pub clutter_density: f64,
    /// Cap on joint events per cluster before falling back to independent PDA.
    pub max_events: usize,
    /// A measurement whose total association probability over all tracks is below this starts a new
    /// tentative track.
    pub birth_threshold: f64,
    /// Merge tracks that are statistically the same target (see `GnnTracker`); this prevents duplicate tracks
    /// but also fuses real targets that fly together, so switch it off for formations.
    pub merge_duplicates: bool,
}

impl JpdaConfig {
    pub fn new(pd: f64, clutter_density: f64) -> Self {
        Self {
            pd,
            clutter_density,
            max_events: 5_000,
            birth_threshold: 0.5,
            merge_duplicates: true,
        }
    }
}

/// Association probabilities.
#[derive(Debug, Clone)]
pub struct Marginals {
    /// `beta[i][j]`: probability that measurement `j` belongs to track `i` (0 when not gated).
    pub beta: Vec<Vec<f64>>,
    /// `miss[i]`: probability that track `i` was not detected (or its detection is outside the gate).
    pub miss: Vec<f64>,
    /// False if any cluster exceeded the event cap and used the independent fallback.
    pub exact: bool,
}

/// Joint association probabilities from detection weights.
///
/// `weights[i][j]` is `Some(w_ij)` for a gated pair and `None` otherwise; `miss_weight` is `1 - Pd`.
pub fn marginals(weights: &[Vec<Option<f64>>], miss_weight: f64, max_events: usize) -> Marginals {
    let n = weights.len();
    let m = weights.first().map_or(0, |r| r.len());
    let mut beta = vec![vec![0.0; m]; n];
    let mut miss = vec![1.0; n];
    let mut exact = true;

    for cluster in clusters(weights) {
        let tracks = &cluster;
        match enumerate(weights, tracks, miss_weight, max_events) {
            Some((b, mi)) => {
                for (k, &i) in tracks.iter().enumerate() {
                    beta[i].clone_from(&b[k]);
                    miss[i] = mi[k];
                }
            }
            None => {
                exact = false;
                for &i in tracks {
                    let total = miss_weight + weights[i].iter().flatten().sum::<f64>();
                    for (j, w) in weights[i].iter().enumerate() {
                        beta[i][j] = w.map_or(0.0, |w| w / total);
                    }
                    miss[i] = miss_weight / total;
                }
            }
        }
    }
    Marginals { beta, miss, exact }
}

/// Connected components of tracks that share at least one gated measurement.
// column-wise traversal of a row-major matrix, so index loops are the clear form
#[allow(clippy::needless_range_loop)]
fn clusters(weights: &[Vec<Option<f64>>]) -> Vec<Vec<usize>> {
    let n = weights.len();
    let mut parent: Vec<usize> = (0..n).collect();
    fn find(p: &mut [usize], x: usize) -> usize {
        let mut r = x;
        while p[r] != r {
            r = p[r];
        }
        let mut y = x;
        while p[y] != r {
            let next = p[y];
            p[y] = r;
            y = next;
        }
        r
    }
    let m = weights.first().map_or(0, |r| r.len());
    for j in 0..m {
        let mut first: Option<usize> = None;
        for i in 0..n {
            if weights[i][j].is_some() {
                match first {
                    None => first = Some(i),
                    Some(f) => {
                        let (a, b) = (find(&mut parent, f), find(&mut parent, i));
                        parent[a] = b;
                    }
                }
            }
        }
    }
    let mut groups: Vec<Vec<usize>> = Vec::new();
    let mut root_index = vec![usize::MAX; n];
    for i in 0..n {
        let r = find(&mut parent, i);
        if root_index[r] == usize::MAX {
            root_index[r] = groups.len();
            groups.push(Vec::new());
        }
        groups[root_index[r]].push(i);
    }
    groups
}

/// Exact marginals for one cluster by depth-first enumeration of joint events; `None` if the number of
/// events exceeds `max_events`.
fn enumerate(
    weights: &[Vec<Option<f64>>],
    tracks: &[usize],
    miss_weight: f64,
    max_events: usize,
) -> Option<(Vec<Vec<f64>>, Vec<f64>)> {
    let m = weights.first().map_or(0, |r| r.len());
    let k = tracks.len();
    let mut beta = vec![vec![0.0; m]; k];
    let mut miss = vec![0.0; k];
    let mut used = vec![false; m];
    let mut choice: Vec<Option<usize>> = vec![None; k];
    let mut total = 0.0;
    let mut events = 0usize;

    #[allow(clippy::too_many_arguments)]
    fn rec(
        t: usize,
        weight: f64,
        weights: &[Vec<Option<f64>>],
        tracks: &[usize],
        miss_weight: f64,
        used: &mut [bool],
        choice: &mut [Option<usize>],
        beta: &mut [Vec<f64>],
        miss: &mut [f64],
        total: &mut f64,
        events: &mut usize,
        max_events: usize,
    ) -> bool {
        if t == tracks.len() {
            *events += 1;
            if *events > max_events {
                return false;
            }
            *total += weight;
            for (k, c) in choice.iter().enumerate() {
                match c {
                    Some(j) => beta[k][*j] += weight,
                    None => miss[k] += weight,
                }
            }
            return true;
        }
        // track t missed
        choice[t] = None;
        if !rec(
            t + 1,
            weight * miss_weight,
            weights,
            tracks,
            miss_weight,
            used,
            choice,
            beta,
            miss,
            total,
            events,
            max_events,
        ) {
            return false;
        }
        // track t takes a free gated measurement
        for j in 0..weights[tracks[t]].len() {
            if let Some(w) = weights[tracks[t]][j] {
                if !used[j] {
                    used[j] = true;
                    choice[t] = Some(j);
                    let ok = rec(
                        t + 1,
                        weight * w,
                        weights,
                        tracks,
                        miss_weight,
                        used,
                        choice,
                        beta,
                        miss,
                        total,
                        events,
                        max_events,
                    );
                    used[j] = false;
                    choice[t] = None;
                    if !ok {
                        return false;
                    }
                }
            }
        }
        true
    }

    let ok = rec(
        0,
        1.0,
        weights,
        tracks,
        miss_weight,
        &mut used,
        &mut choice,
        &mut beta,
        &mut miss,
        &mut total,
        &mut events,
        max_events,
    );
    if !ok || total <= 0.0 || !total.is_finite() {
        return None;
    }
    for row in &mut beta {
        for b in row.iter_mut() {
            *b /= total;
        }
    }
    for x in &mut miss {
        *x /= total;
    }
    Some((beta, miss))
}

/// Update `ekf` with the gated measurements of one track. `inns[n]` is `(beta_n, innovation_n)`; `beta0` is
/// the probability that none of them is the track's measurement. All innovations share one `H` and `S` (one
/// linearisation per track), so the gain is computed once.
pub fn pda_update(ekf: &mut Ekf, radar: &Radar, inns: &[(f64, Innovation)], beta0: f64) {
    let Some((_, first)) = inns.first() else {
        return;
    };
    let rm = Matrix2::new(radar.sigma_r.powi(2), 0.0, 0.0, radar.sigma_b.powi(2));
    let k: Matrix4x2<f64> = ekf.p * first.h.transpose() * first.s_inv;

    let mut v = Vector2::zeros();
    let mut spread = Matrix2::zeros();
    for (b, inn) in inns {
        v += *b * inn.y;
        spread += *b * inn.y * inn.y.transpose();
    }
    let p_pred = ekf.p;
    // covariance after an update with the correct measurement (Joseph form, as in `Ekf::update_with`)
    let i_kh = Matrix4::identity() - k * first.h;
    let p_correct = i_kh * p_pred * i_kh.transpose() + k * rm * k.transpose();
    let p_tilde = k * (spread - v * v.transpose()) * k.transpose();

    ekf.x += k * v;
    let p = beta0 * p_pred + (1.0 - beta0) * p_correct + p_tilde;
    ekf.p = 0.5 * (p + p.transpose()); // keep it symmetric against rounding
}

#[cfg(test)]
mod tests {
    use super::*;

    fn close(a: f64, b: f64) -> bool {
        (a - b).abs() < 1e-12
    }

    #[test]
    fn one_track_one_measurement_is_the_pda_ratio() {
        let m = marginals(&[vec![Some(3.0)]], 0.2, 1000);
        assert!(close(m.beta[0][0], 3.0 / 3.2) && close(m.miss[0], 0.2 / 3.2) && m.exact);
    }

    #[test]
    fn two_tracks_sharing_a_measurement_compete_for_it() {
        // events: both miss (m*m), t0 takes z (w0*m), t1 takes z (m*w1)
        let (w0, w1, m) = (4.0, 1.0, 0.1);
        let r = marginals(&[vec![Some(w0)], vec![Some(w1)]], m, 1000);
        let total = m * m + w0 * m + m * w1;
        assert!(close(r.beta[0][0], w0 * m / total));
        assert!(close(r.beta[1][0], m * w1 / total));
        assert!(
            r.beta[0][0] + r.beta[1][0] <= 1.0,
            "a measurement cannot belong to two tracks"
        );
        assert!(close(r.miss[0], (m * m + m * w1) / total));
    }

    #[test]
    fn exclusivity_changes_the_answer_relative_to_independent_pda() {
        let weights = vec![vec![Some(5.0), None], vec![Some(5.0), Some(1.0)]];
        let joint = marginals(&weights, 0.1, 1000);
        let indep = marginals(&weights, 0.1, 0); // cap of 0 forces the independent fallback
        assert!(joint.exact && !indep.exact);
        assert!(
            joint.beta[0][0] < indep.beta[0][0],
            "sharing a measurement lowers each track's claim to it"
        );
    }

    #[test]
    fn disjoint_tracks_are_independent_pdas() {
        let weights = vec![vec![Some(2.0), None], vec![None, Some(3.0)]];
        let r = marginals(&weights, 0.5, 1000);
        assert!(close(r.beta[0][0], 2.0 / 2.5) && close(r.beta[1][1], 3.0 / 3.5));
        assert_eq!(clusters(&weights).len(), 2);
    }

    #[test]
    fn probabilities_of_each_track_sum_to_one_and_each_measurement_to_at_most_one() {
        let w = |x| Some(x);
        let weights = vec![
            vec![w(1.0), w(2.0), None],
            vec![w(3.0), None, w(0.5)],
            vec![None, w(1.5), w(2.5)],
        ];
        let r = marginals(&weights, 0.3, 10_000);
        for i in 0..3 {
            assert!(close(r.miss[i] + r.beta[i].iter().sum::<f64>(), 1.0));
        }
        for j in 0..3 {
            assert!((0..3).map(|i| r.beta[i][j]).sum::<f64>() <= 1.0 + 1e-12);
        }
    }

    #[test]
    fn an_ungated_track_has_miss_probability_one() {
        let r = marginals(&[vec![None, None]], 0.4, 100);
        assert!(close(r.miss[0], 1.0) && r.beta[0] == vec![0.0, 0.0]);
    }
}
