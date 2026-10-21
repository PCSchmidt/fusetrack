//! Multi-target evaluation metrics.

use nalgebra::Vector2;

use crate::assignment;

/// OSPA distance (Schuhmacher, Vo & Vo 2008) between two point sets.
///
/// `c` is the cutoff (m) that caps per-target localisation error and sets the
/// penalty for a missed or false track; `p` is the order. Returns 0 for two
/// empty sets.
pub fn ospa(truth: &[Vector2<f64>], est: &[Vector2<f64>], c: f64, p: f64) -> f64 {
    let (small, large) = if truth.len() <= est.len() {
        (truth, est)
    } else {
        (est, truth)
    };
    let (m, n) = (small.len(), large.len());
    if n == 0 {
        return 0.0;
    }

    let cost: Vec<Vec<f64>> = small
        .iter()
        .map(|a| {
            large
                .iter()
                .map(|b| (a - b).norm().min(c).powf(p))
                .collect()
        })
        .collect();
    let assigned = assignment::solve(&cost);
    let loc: f64 = assigned
        .iter()
        .enumerate()
        .filter_map(|(i, j)| j.map(|j| cost[i][j]))
        .sum();
    let card = c.powf(p) * (n - m) as f64;
    ((loc + card) / n as f64).powf(1.0 / p)
}

/// GOSPA distance (Rahmathullah, Garcia-Fernandez & Svensson 2017) with `alpha = 2`.
///
/// Unlike OSPA it penalises missed and false tracks separately: the result is
/// `(sum over assigned pairs of min(d, c)^p + c^p / 2 * (missed + false))^(1/p)`, where optimal assignment
/// pairs `min(|truth|, |est|)` targets with tracks. Returns 0 for two empty sets.
pub fn gospa(truth: &[Vector2<f64>], est: &[Vector2<f64>], c: f64, p: f64) -> f64 {
    let (small, large) = if truth.len() <= est.len() {
        (truth, est)
    } else {
        (est, truth)
    };
    if large.is_empty() {
        return 0.0;
    }
    let cost: Vec<Vec<f64>> = small
        .iter()
        .map(|a| {
            large
                .iter()
                .map(|b| (a - b).norm().min(c).powf(p))
                .collect()
        })
        .collect();
    let assigned = assignment::solve(&cost);
    let paired: f64 = assigned
        .iter()
        .enumerate()
        .filter_map(|(i, j)| j.map(|j| cost[i][j]))
        .sum();
    // each point left unpaired, in either set, costs c^p / 2: |X| + |Y| - 2k of them for k pairs
    let k = assigned.iter().flatten().count();
    let unpaired = small.len() + large.len() - 2 * k;
    (paired + 0.5 * c.powf(p) * unpaired as f64).powf(1.0 / p)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn pt(x: f64, y: f64) -> Vector2<f64> {
        Vector2::new(x, y)
    }

    #[test]
    fn identical_sets_have_zero_distance() {
        let a = vec![pt(0.0, 0.0), pt(100.0, 50.0)];
        assert!(ospa(&a, &a, 200.0, 2.0) < 1e-12);
        assert_eq!(ospa(&[], &[], 200.0, 2.0), 0.0);
    }

    #[test]
    fn one_empty_set_gives_the_cutoff() {
        let a = vec![pt(0.0, 0.0), pt(1.0, 1.0)];
        assert!((ospa(&a, &[], 200.0, 2.0) - 200.0).abs() < 1e-9);
        assert!((ospa(&[], &a, 200.0, 2.0) - 200.0).abs() < 1e-9);
    }

    #[test]
    fn cardinality_error_is_penalised() {
        let truth = vec![pt(0.0, 0.0)];
        let est = vec![pt(0.0, 0.0), pt(9999.0, 0.0)];
        // ((0 + c^2) / 2)^(1/2) = c / sqrt(2)
        let d = ospa(&truth, &est, 200.0, 2.0);
        assert!((d - 200.0 / 2.0_f64.sqrt()).abs() < 1e-9);
    }

    #[test]
    fn localisation_error_is_capped_at_the_cutoff() {
        let d = ospa(&[pt(0.0, 0.0)], &[pt(5000.0, 0.0)], 200.0, 2.0);
        assert!((d - 200.0).abs() < 1e-9);
        let d = ospa(&[pt(0.0, 0.0)], &[pt(30.0, 40.0)], 200.0, 2.0);
        assert!((d - 50.0).abs() < 1e-9);
    }

    #[test]
    fn gospa_basics() {
        let a = [Vector2::new(0.0, 0.0), Vector2::new(100.0, 0.0)];
        assert_eq!(gospa(&[], &[], 200.0, 2.0), 0.0);
        assert!(gospa(&a, &a, 200.0, 2.0) < 1e-12);
        // one missed target: c^p / 2 for it
        let one = [Vector2::new(0.0, 0.0)];
        let g = gospa(&a, &one, 200.0, 2.0);
        assert!((g - (0.5f64 * 200.0f64.powi(2)).sqrt()).abs() < 1e-9, "{g}");
        // far apart (beyond the cutoff): a pair costs c^p, the same as one miss plus one false
        let far = [Vector2::new(1.0e5, 0.0)];
        let g = gospa(&one, &far, 200.0, 2.0);
        assert!((g - 200.0).abs() < 1e-9, "{g}");
        // GOSPA separates cardinality error from localisation error more than OSPA does: both missed and false count
        assert!(gospa(&a, &[], 200.0, 2.0) > gospa(&a, &one, 200.0, 2.0));
    }
}
