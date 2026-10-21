# JPDA (milestone 0.4)

JPDA is implemented, tested and measured. **On the scenarios here it does not beat the GNN tracker; in clutter it is worse, and at close separation it makes more identity errors.** This page records the method, the checks, and the negative result with its likely causes, which I have only partly tested.

## What was built

- `jpda.rs`: joint association probabilities by exact enumeration of joint events within each cluster of tracks that share gated measurements, with a cap on events per cluster (default 5,000) beyond which a cluster falls back to independent per-track PDA and the scan is counted in `jpda_inexact_scans`; and the PDA state and covariance update (combined innovation plus the spread-of-the-means term).
- `tracker.rs`: `TrackerConfig::jpda` swaps the Hungarian assignment for JPDA. Single-EKF tracks only; combining it with the IMM is refused with an error. With score-based confirmation, a track's score uses the PDA likelihood ratio `(1 - Pd) + sum_j Pd N_j / lambda_j`; a hit is a scan on which the track's miss probability is at most one half; a measurement starts a new tentative track when its total association probability across tracks is below a threshold (0.5). An optional duplicate-track merge drops the weaker of two statistically identical tracks.
- Evaluation: GOSPA (alpha = 2), ID switches and fragmentations (`eval::evaluate_quality`), and two scripted scenarios, `simulate_crossing` (two targets crossing at 90 degrees with a set closest approach) and `simulate_parallel` (two targets side by side, same velocity, a set separation).

## Checks that pass

- The marginals match hand-computed values for small cases; each track's probabilities sum to one; no measurement is claimed by more than one track in total; disjoint tracks reduce to independent PDA; forcing the cap reproduces the independent fallback and differs from the joint answer when tracks share a measurement.
- With one certain measurement the PDA update equals the existing EKF update to 1e-9; with no detection it changes nothing; with two equally likely measurements the covariance does not shrink and symmetric innovations cancel.
- **Consistency in clutter.** For a single target, the mean normalised estimation error squared (a consistent 4-state filter gives 4), with the tracker's process noise matched to the scenario: JPDA 3.87 / 3.94 / 4.00 and GNN 3.87 / 3.94 / 4.08 at 0 / 50 / 300 clutter returns per scan (4,800 scored scans each; the 300 case from an earlier, longer run). Both are consistent and nearly identical, which is expected for one isolated target. With the tracker's default process noise (larger than the scenario's) both read 2.5: a pessimistic filter, not a bug.
- GOSPA properties (zero for identical sets, one missed target costs `c / sqrt(2)`, a pair beyond the cutoff costs `c`).

## The experiments

Two targets, Pd 0.9, a 15 km disc, 100 scans, score-based confirmation for both trackers (confirm 12, delete -5) with the true Pd and clutter density, 15 to 40 seeds per cell (15 at the heaviest clutter, where JPDA is slow). Full output: `results/jpda_vs_gnn.txt`, `results/jpda_formation.txt`, `results/jpda_confirm_sweep.txt`.

**Crossing at 90 degrees** (closest approach 25 to 500 m). The two targets differ by 35 m/s in velocity, so a tracker can tell them apart from their motion alone, and ID switches are rare for both (about 0.0 to 0.4 per run, within the noise). OSPA, mean over closest approaches:

| Clutter per scan | GNN OSPA | JPDA OSPA | JPDA tracks per scan (GNN 2.00 to 2.02) |
| --- | --- | --- | --- |
| 0 | 12.8 m | 12.8 m | 2.00 |
| 20 | 12.7 m | 14.1 m | 2.01 to 2.02 |
| 100 | 12.6 m | 24.6 m | 2.10 to 2.14 |
| 300 | 14.9 m | 32.6 m | 2.16 to 2.25 |

**Side by side, same velocity** (separation 30 to 240 m; cross-range resolution at 6 km is about 18 m). ID switches per run, no clutter:

| Separation | GNN | JPDA | JPDA + duplicate merge (fragmentations per run) |
| --- | --- | --- | --- |
| 30 m | 1.23 | 2.67 | 1.30 (2.90) |
| 60 m | 0.55 | 1.93 | 1.48 (1.73) |
| 120 m | 0.28 | 0.45 | 0.42 (0.20) |
| 240 m | 0.05 | 0.23 | 0.15 (0.12) |

At 100 clutter returns per scan and 30 m separation: GNN OSPA 12.3 m, 1.15 ID switches; JPDA OSPA 25.4 m, 2.77 switches, 0.09 false tracks per scan.

## What this says

- **In clutter JPDA confirms extra false tracks that GNN rejects.** Its OSPA is about twice GNN's at 100 returns per scan. Raising JPDA's confirmation threshold narrows the gap but does not close it (clutter 100, formation 120 m): confirm 12, OSPA 26.7 m and 0.12 false tracks; confirm 20, 21.5 m and 0.08; confirm 30, 17.1 m and 0.04; GNN at 12: 12.6 m and 0.00. So part of the excess is the PDA score being more permissive than the best-measurement score at the same threshold, which is a calibration issue; the rest I have not isolated. Delay to confirmation at the higher thresholds was not measured.
- **At close separation JPDA has more identity switches than GNN, not fewer.** This is the opposite of the textbook motivation. Possible reasons, not tested: coalescence of the two tracks (the classic JPDA failure when targets stay together), the PDA update's covariance inflation letting tracks drift onto the neighbour, and the merge rule, which fixes switches by fusing and fragmenting the real pair.
- **The duplicate merge is a blunt tool.** It lowers ID switches at 30 to 60 m but splits the formation (about 1.5 to 3.4 fragmentations per run at 30 to 60 m) because two real, nearby targets with similar velocity are statistically indistinguishable.
- **The crossing scenario was too easy to show a benefit** (velocity separates the targets), and the side-by-side one shows a cost. A scenario with targets that cross with similar velocities, lower Pd, or a longer time in a shared gate is where JPDA would be expected to help, and I have not run it.

## Limits

Simulated, two targets, Pd 0.9, a known clutter density and Pd, 15 to 40 seeds per cell without confidence intervals, single-EKF tracks only. JPDA's exact enumeration is expensive: the heaviest clutter case took minutes per cell and fell back to independent PDA on up to 0.1 scans per run (the cap is 5,000 events per cluster). Variants that address coalescence (for example JPDA*, or coupled-filter formulations), a tuned birth and confirmation rule for JPDA, and the IMM combination are not tried.
