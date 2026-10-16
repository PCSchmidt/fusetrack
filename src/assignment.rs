//! Optimal linear assignment (Hungarian / Kuhn-Munkres, O(n^2 m)).

/// Cost used to mark a pairing as forbidden (e.g. outside the gate).
pub const FORBIDDEN: f64 = 1.0e9;

/// Minimise total cost over one-to-one assignments of rows to columns.
///
/// `cost` is an `n x m` matrix (rows of equal length). Returns, for each row,
/// the assigned column or `None`. Every row is assigned when `n <= m`; with
/// `n > m` the surplus rows are left unassigned. Pairings whose cost is at
/// least [`FORBIDDEN`] are never returned.
pub fn solve(cost: &[Vec<f64>]) -> Vec<Option<usize>> {
    let n = cost.len();
    let m = cost.first().map_or(0, |r| r.len());
    if n == 0 || m == 0 {
        return vec![None; n];
    }

    let mut out = if n <= m {
        solve_wide(cost, n, m)
    } else {
        // Transpose so rows <= cols, solve, then invert the mapping.
        let t: Vec<Vec<f64>> = (0..m)
            .map(|j| (0..n).map(|i| cost[i][j]).collect())
            .collect();
        let col_to_row = solve_wide(&t, m, n);
        let mut res = vec![None; n];
        for (j, assigned) in col_to_row.iter().enumerate() {
            if let Some(i) = assigned {
                res[*i] = Some(j);
            }
        }
        res
    };

    for (i, a) in out.iter_mut().enumerate() {
        if let Some(j) = *a {
            if cost[i][j] >= FORBIDDEN {
                *a = None;
            }
        }
    }
    out
}

/// Potentials-based Hungarian algorithm, requires `n <= m`.
fn solve_wide(cost: &[Vec<f64>], n: usize, m: usize) -> Vec<Option<usize>> {
    debug_assert!(n <= m);
    let inf = f64::INFINITY;
    // 1-indexed; row 0 / column 0 are sentinels.
    let mut u = vec![0.0; n + 1];
    let mut v = vec![0.0; m + 1];
    let mut p = vec![0usize; m + 1];
    let mut way = vec![0usize; m + 1];

    for i in 1..=n {
        p[0] = i;
        let mut j0 = 0;
        let mut minv = vec![inf; m + 1];
        let mut used = vec![false; m + 1];
        loop {
            used[j0] = true;
            let i0 = p[j0];
            let mut delta = inf;
            let mut j1 = 0;
            for j in 1..=m {
                if !used[j] {
                    let cur = cost[i0 - 1][j - 1] - u[i0] - v[j];
                    if cur < minv[j] {
                        minv[j] = cur;
                        way[j] = j0;
                    }
                    if minv[j] < delta {
                        delta = minv[j];
                        j1 = j;
                    }
                }
            }
            for j in 0..=m {
                if used[j] {
                    u[p[j]] += delta;
                    v[j] -= delta;
                } else {
                    minv[j] -= delta;
                }
            }
            j0 = j1;
            if p[j0] == 0 {
                break;
            }
        }
        loop {
            let j1 = way[j0];
            p[j0] = p[j1];
            j0 = j1;
            if j0 == 0 {
                break;
            }
        }
    }

    let mut ans = vec![None; n];
    for j in 1..=m {
        if p[j] != 0 {
            ans[p[j] - 1] = Some(j - 1);
        }
    }
    ans
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::{rngs::StdRng, Rng, SeedableRng};

    fn total(cost: &[Vec<f64>], a: &[Option<usize>]) -> f64 {
        a.iter()
            .enumerate()
            .filter_map(|(i, j)| j.map(|j| cost[i][j]))
            .sum()
    }

    /// Exhaustive minimum over assignments of all rows (requires n <= m).
    fn brute_force(cost: &[Vec<f64>]) -> f64 {
        fn rec(cost: &[Vec<f64>], row: usize, used: &mut Vec<bool>) -> f64 {
            if row == cost.len() {
                return 0.0;
            }
            let mut best = f64::INFINITY;
            for j in 0..used.len() {
                if !used[j] {
                    used[j] = true;
                    best = best.min(cost[row][j] + rec(cost, row + 1, used));
                    used[j] = false;
                }
            }
            best
        }
        rec(cost, 0, &mut vec![false; cost[0].len()])
    }

    #[test]
    fn matches_brute_force_on_random_matrices() {
        let mut rng = StdRng::seed_from_u64(3);
        for _ in 0..200 {
            let n = rng.gen_range(1..=5);
            let m = rng.gen_range(n..=6);
            let cost: Vec<Vec<f64>> = (0..n)
                .map(|_| (0..m).map(|_| rng.gen_range(-5.0..20.0)).collect())
                .collect();
            let a = solve(&cost);
            assert!(a.iter().all(|x| x.is_some()));
            let mut cols: Vec<usize> = a.iter().map(|x| x.unwrap()).collect();
            cols.sort_unstable();
            cols.dedup();
            assert_eq!(cols.len(), n, "columns must be distinct");
            assert!((total(&cost, &a) - brute_force(&cost)).abs() < 1e-9);
        }
    }

    #[test]
    fn handles_more_rows_than_columns() {
        let mut rng = StdRng::seed_from_u64(4);
        for _ in 0..100 {
            let n = rng.gen_range(2..=6);
            let m = rng.gen_range(1..n);
            let cost: Vec<Vec<f64>> = (0..n)
                .map(|_| (0..m).map(|_| rng.gen_range(0.0..20.0)).collect())
                .collect();
            let a = solve(&cost);
            assert_eq!(a.iter().filter(|x| x.is_some()).count(), m);
            // Compare with the transposed problem solved directly.
            let t: Vec<Vec<f64>> = (0..m)
                .map(|j| (0..n).map(|i| cost[i][j]).collect())
                .collect();
            assert!((total(&cost, &a) - brute_force(&t)).abs() < 1e-9);
        }
    }

    #[test]
    fn forbidden_pairs_are_never_returned() {
        let cost = vec![vec![1.0, FORBIDDEN], vec![FORBIDDEN, FORBIDDEN]];
        let a = solve(&cost);
        assert_eq!(a, vec![Some(0), None]);
    }

    #[test]
    fn empty_inputs() {
        assert!(solve(&[]).is_empty());
        assert_eq!(solve(&[vec![], vec![]]), vec![None, None]);
    }
}
