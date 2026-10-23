# JPDA (milestone 0.4)

JPDA is implemented, tested and measured. **On the scenarios here it does not beat the GNN tracker; in clutter it is worse, and at close separation it makes more identity errors.** This page records the method, the checks, and the negative result with its likely causes, which I have only partly tested.

## What was built

- `jpda.rs`: joint association probabilities by exact enumeration of joint events within each cluster of tracks that share gated measurements, with a cap on events per cluster (default 5,000) beyond which a cluster falls back to independent per-track PDA and the scan is counted in `jpda_inexact_scans`; and the PDA state and covariance update (combined innovation plus the spread-of-the-means term).
- `tracker.rs`: `TrackerConfig::jpda` swaps the Hungarian assignment for JPDA. Works with single-EKF tracks and, since the last study below, with IMM tracks. With score-based confirmation, a track's score uses the PDA likelihood ratio `(1 - Pd) + sum_j Pd N_j / lambda_j`; a hit is a scan on which the track's miss probability is at most one half; a measurement starts a new tentative track when its total association probability across tracks is below a threshold (0.5). An optional duplicate-track merge drops the weaker of two statistically identical tracks.
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
- **The crossing scenario was too easy to show a benefit** (velocity separates the targets), and the side-by-side one shows a cost. The shallow-angle scenario below is the case where JPDA is expected to help; it did not.

## The harder scenario: shallow-angle crossings

Targets cross with headings 5 to 40 degrees apart (relative speed 2.2 to 17 m/s at 25 m/s each), so they share a gate for many scans, with closest approach 0 or 40 m, Pd 0.9 or 0.7, and 0 or 20 clutter returns per scan. Both trackers see identical scans on each of 60 seeds, so the table gives paired differences (JPDA minus GNN, mean with a 95% interval of 1.96 standard errors; positive is worse for JPDA). Full output: `results/jpda_shallow.txt`, from `examples/jpda_shallow.rs`.

| Pd | Clutter | Angle | GOSPA difference | ID-switch difference per run |
| --- | --- | --- | --- | --- |
| 0.9 | 0 | 5 deg | +1.5 to +2.3 m (+-1.3) | +0.22 (+-0.16) at 0 m, none at 40 m |
| 0.9 | 0 | 10 to 40 deg | +0.0 to +1.2 m (+-0.4 to 0.9) | +0.07 to +0.20, mostly within the interval |
| 0.9 | 20 | 5 to 40 deg | +2.5 to +4.7 m (+-1.0 to 1.7) | +0.00 to +0.37, mostly within the interval |
| 0.7 | 0 | 5 deg | +2.1 to +2.9 m (+-1.7 to 1.9) | +0.03 to +0.30 (+-0.2 to 0.3) |
| 0.7 | 0 | 10 to 40 deg | -0.1 to +1.3 m | +0.02 to +0.18 |
| 0.7 | 20 | 5 to 40 deg | +4.2 to +7.8 m (+-1.6 to 2.4) | -0.02 to +0.60, mostly within the interval |

**JPDA did not beat GNN in any of the 32 cells.** The GOSPA difference is never significantly negative and is significantly positive (interval excludes zero) in 23 of 32 cells, including every clutter cell. Identity switches are higher for JPDA in 28 of 32 cells and significantly higher in 11; none is significantly lower. With 32 cells at 95% a few chance exceedances are expected, but the direction is consistent. Where JPDA is worse in GOSPA the cause is probably the clutter result above (extra false tracks), not identity errors; I did not separate the two in this run.

Why the scenario may still not be hard enough, which I have not tested: GNN's own identity-switch rate is low here (0.03 to 0.6 per run), so there is little for JPDA to remove, and the simulator never merges two targets into one measurement (a real radar has a resolution cell below which two returns become one). Ambiguity here comes only from measurement noise (about 15 m in range, 15 m in cross-range at 5 km), so at 0 m closest approach the targets are ambiguous for a few scans and no more. A merged-measurement or unresolved-target model is the setting where JPDA's soft association is usually argued to matter, and it is not built.

### Finite resolution (merged returns)

The remaining untested explanation was that the simulator never merges two close targets into one detection. `sim::merge_unresolved` now does: any two returns closer than a resolution cell (range x cross-range) become one detection at their mean, repeatedly, clutter included. Six scenarios (shallow crossings at 5, 10 and 20 degrees with 0 m closest approach; side by side at 30, 60 and 120 m), cells of none, 30 x 40 m and 60 x 80 m, clutter 0 and 20, Pd 0.9, 60 paired seeds, 36 cells. Full output: `results/jpda_merged.txt`, from `examples/jpda_merged.rs`.

- **The merging bites.** At the 60 x 80 m cell, GNN's GOSPA rises from about 18 to 24 to 32 m in the close scenarios (3.5 to 4 fragmentations per run at 30 m separation), and the unresolved scans per run range from 1 to 15.
- **JPDA still never wins.** GOSPA difference significantly positive in 27 of 36 cells, never negative (not one cell has a negative mean). ID switches higher for JPDA in 30 of 36, significantly higher in 13, significantly lower in none.
- **The gap narrows where merging is heaviest.** Side by side at 30 m, no clutter: +8.0 m (+-2.2) with no merging, +7.0 with the 30 x 40 m cell, +2.1 (+-1.2) with 60 x 80 m. That is the direction one would expect if JPDA's soft association helps when returns are ambiguous, but it closes the gap only because GNN also gets worse; it does not reverse it.

This was the last scenario change I had in mind for the question "when does JPDA beat GNN here". The honest answer from this study is: not in any of about 100 cells across crossings, formations, two detection probabilities, clutter, and merged returns. That is a statement about this implementation and these scenarios (known Pd and clutter density, exact enumeration with a cap, track confirmation by score, two targets), not about JPDA in general; a better-tuned JPDA, a coalescence-avoiding variant, or a larger target count could differ.

### JPDA with the IMM

The one combination I had refused to build is now built: each mode filter of an IMM track gets a PDA update with its own innovations and the association probabilities from the joint events, and each mode's probability is reweighted by its PDA likelihood of the whole measurement set, `(1 - Pd) + sum Pd N / lambda`. The association weights themselves come from the track's combined (moment-matched) prediction, an approximation; a mixture over the modes would be more exact. Unit tests: with one certain measurement it equals the plain IMM update (state and mode probabilities); with no measurement it changes nothing; and on a noise-free turning target the turning modes take over (over 80% of the probability within 25 scans).

Study (`examples/jpda_imm.rs`, `results/jpda_imm.txt`): two targets crossing at 10 or 20 degrees with closest approach 0 or 40 m, one of them turning for 20 s as they meet (gentle: 0.04 rad/s, about 1 m/s^2 sideways; hard: 0.10 rad/s, about 2.5 m/s^2), Pd 0.9 and 0.7, clutter 0, 20 and 60 returns per scan, 40 paired seeds per cell, 48 cells per comparison. Cell counts below are "interval excludes zero".

| Comparison | GOSPA significantly better | GOSPA significantly worse | ID switches better / worse |
| --- | --- | --- | --- |
| IMM under GNN (gentle turn) | 4 of 24 | 0 | 0 / 0 |
| IMM under GNN (hard turn) | **24 of 24** | 0 | 7 / 0 |
| IMM under JPDA | 19 of 48 | 7 | 4 / 1 |
| **JPDA given the IMM** (GNN+IMM to JPDA+IMM) | **0 of 48** | **38** | 0 / 12 |

- **The IMM works.** With a hard turn it improves the GNN tracker in every cell (mean GOSPA difference -1.8 to -2.5 m), so the manoeuvre model is doing its job; with a gentle turn the single filter's process noise already absorbs it and the gain is about 0.2 m.
- **JPDA still does not beat GNN once both have the IMM.** The mean GOSPA penalty of JPDA+IMM against GNN+IMM is +0.5 to +0.6 m with no clutter, +3.4 to +3.7 m at 20 returns per scan and +17 m at 60; identity switches are more frequent in 12 cells and fewer in none. It is the same result as before, now with the best manoeuvre handling I have.
- **The IMM helps JPDA more often than not** (19 cells better, 7 worse), so some of JPDA's deficit in the earlier studies was missing manoeuvre handling, but not the part that matters here: the extra false tracks in clutter.

## Limits

Simulated, two targets, Pd 0.9, a known clutter density and Pd, 15 to 40 seeds per cell without confidence intervals, single-EKF tracks in the main studies (the IMM combination is studied above). JPDA's exact enumeration is expensive: the heaviest clutter case took minutes per cell and fell back to independent PDA on up to 0.1 scans per run (the cap is 5,000 events per cluster). The later scenarios add Pd 0.7 and merged returns but still use known Pd and clutter density and two targets. Variants that address coalescence (for example JPDA*, or coupled-filter formulations), and a tuned birth and confirmation rule for JPDA are not tried.
