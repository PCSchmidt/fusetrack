# FuseTrack: Evidence Notes

Companion to fusetrack.md. Labels: Confirmed, Ran it, Per the docs, My reading, Ask the author. Depth: LITE-plus (11 commits, results present). Repo state read: main at c5978a2 (2026-10-23). Read-only clone; builds were directed to a target folder in the scratchpad, so nothing was written to the clone.

## E1 Fact sheet

- Name: fusetrack; public GitHub repo PCSchmidt/fusetrack; MIT (Confirmed, LICENSE and GitHub metadata). 0 stars, 0 forks. Created 2026-10-16 17:44Z, last push 2026-10-23 15:04Z.
- Kind: Rust library plus 11 example programs (simulation study). Status per README: "v0.4, in progress".
- Size: about 3,480 lines of Rust in src and tests (excluding examples) (wc -l: src 2,906, tests 576; the guide now states both), plus 11 examples.
- Tests: 58 passed, 0 failed (Ran it: `cargo test --all-targets`; 39 unit tests in src, 2 + 5 + 2 + 10 in the four integration files; example binaries have none).
- Commits: 11, 2026-10-16 to 2026-10-23. All 11 carry an AI co-author trailer (Confirmed). No agent instruction files or settings folders in the repo (Confirmed).
- Runs from a clean copy: yes with a stable Rust toolchain (cargo 1.99.1 used). Dependencies: nalgebra 0.33, rand 0.8, rand_distr 0.4; Cargo.lock is committed.
- Needs: no paid services, no GPU, no downloads. Heavy examples take minutes (README: clutter_sweep about 100 s, clutter_diagnose about 5 min, jpda_imm about 3 min, jpda_vs_gnn tens of minutes). Writer did not run clutter_sweep, clutter_diagnose, jpda_*, imm_tracker; the independent review re-ran all of them except jpda_vs_gnn, jpda_confirm_sweep and the full-scale jpda_vs_gnn (saved files used), and all matched.
- Deployed: no. Outside users: none known (Ask the author).
- Related: part of physical-ai-portfolio (Per the docs).

## E2 Reading and running log

Read: README, docs/roadmap.md, docs/jpda.md (in full), the head of docs/imm-derivation.md listing only, ci.yml, Cargo.toml, .gitignore, src/tracker.rs (header, config, gating lines), parts of sim.rs, metrics.rs, ekf.rs/imm.rs signatures, results/jpda_vs_gnn.txt (parsed), git log with trailers. Not read in detail: eval.rs, jpda.rs body, tests bodies, other result files, imm-derivation.md body. Commands run: `cargo test --all-targets`; `cargo fmt --all -- --check` (exit 0); `cargo clippy --all-targets -- -D warnings` (finished, no warnings); `cargo run --release --example single_target`; `cargo run --release --example imm_maneuver`. No network use except `gh` queries (run list, repo metadata, code-scanning). No .env exists. `git status` in the clone after all runs: clean (no output).

## E3 Claims check

Levels: L1 doc matches saved file; L2 reproduced by re-running; L3 independently recomputed.

- CL-01 "58 tests" is the guide's own count; README gives none. Ran it, L2.
- CL-02 IMM turning table (README, 9 rows): re-run printed identical overall, straight and turning values for all nine rows (e.g. 64.0, 17.0, 12.6 at omega 0.10; 85.3, 20.3, 19.8 at 0.14). Supported, L2.
- CL-03 "IMM gains 26% at matched, 10% at 0.06, 2% at 0.14 against sigma_a = 3.0". Recomputed: (17.0-12.6)/17.0 = 25.9%; (13.2-11.9)/13.2 = 9.8%; (20.3-19.8)/20.3 = 2.5%. Supported, L3 (the 2% is rounded from 2.5%; the guide says "about 2%").
- CL-04 JPDA vs GNN OSPA means (jpda.md: 12.8/12.8, 12.7/14.1, 12.6/24.6, 14.9/32.6 at clutter 0/20/100/300). Recomputed from results/jpda_vs_gnn.txt: GNN 12.84, 12.66, 12.64, 14.86; JPDA 12.80, 14.14, 24.56, 32.64, five closest-approach rows per level. Supported, L3.
- CL-05 Side-by-side 30 m ID switches 2.67 (JPDA) vs 1.23 (GNN). Per jpda.md table; results/jpda_formation.txt not re-parsed. L1 not checked against the file. Per the docs.
- CL-06 "never better than GNN in 32 cells; significantly worse GOSPA in 23"; "36 more cells, never better, 27 significantly worse"; "JPDA+IMM worse in 38 of 48, better in none". Per the docs (jpda.md). Not recomputed from jpda_shallow.txt, jpda_merged.txt, jpda_imm.txt. Note: the guide quotes the 32-cell, 36-cell and 48-cell counts; it omits the 27-of-36 figure.
- CL-07 Clutter tables, 3-of-5 vs score rule (README): re-run in the independent review (clutter_diagnose, about 5 min); every row printed identical to the README. Ran it, L2. The guide's table rows (154.2, 15.8, 15.1 m; 7.99, 0.01, 0.00; 3.4, 5.7, 4.5 scans) are copied from the README table at clutter 200, Pd 0.9.
- CL-08 "OSPA of 200 for a tracker with no tracks". My reading of metrics.rs: with no estimates, n = number of truths, card = c^p x n, result = c. Matches the code; not run as a probe.
- CL-09 NEES 3.89, NIS 1.97, RMSE 14.2 vs 23.2 m (README). tests/ekf_consistency.rs exists; its assertions were not read. The single_target example output (a different table) was run. Per the docs.
- CL-10 "CI runs fmt, clippy -D warnings and tests" (README). ci.yml confirms. Supported, L1.
- CL-12 Clutter chart (images/fusetrack-clutter.png) values 12.5, 12.8, 22.3, 156.0 m at clutter 0, 10, 50, 200 (Pd 1.0, 5 targets, 20 seeds) are copied from the README table; README-reported and reproduced in the independent review (clutter_sweep re-run: all 12 rows identical). Ran it, L2.
- CL-13 Unit conversions added to the guide: 0.003 rad = 0.17 degrees, about 15 m sideways at 5 km; 0.10 rad/s = 5.7 degrees per second. My reading (arithmetic). "200 seeds" for the IMM chart is from README line 95 (re-run table). The five closest-approach distances are 25, 50, 100, 200 and 500 m (results/jpda_vs_gnn.txt; jpda.md). The 300-clutter point (14.9 vs 32.6) is now cited in the text.
- CL-11 README "gate: 99.9% chi-square, 2 dof = 13.8" (tracker.rs line 129, default 13.8). Confirmed. The guide does not quote the number.

Count: 11 claims checked; 0 contradicted; 5 only against the project's own docs.

## E4 Rule attacks

No bad-input probes were run: the program reads no files, no network and no command-line input (the examples take no arguments, as far as the README run list shows; not verified in each example). Edge behaviour listed in the guide comes from tests and code comments (My reading): measurements outside the gate are FORBIDDEN in the assignment cost (tracker.rs ~line 198 to 208); a confirmed track survives 5 missed scans and is deleted on the 6th (misses > max_misses, default 5; the test at about line 516 asserts "deleted after max_misses + 1"); JPDA caps joint events per cluster at 5,000 and falls back to independent PDA, counted in jpda_inexact_scans (jpda.md). Not attacked: negative or zero variances, empty scans, coincident targets, NaN measurements.

## E5 Code map

src/: sim.rs (truth, scenarios, crossing and parallel scenarios, merge_unresolved), models.rs (radar), ekf.rs (EKF, Joseph form, measurement model), imm.rs (IMM, per-mode PDA update), assignment.rs (Hungarian), tracker.rs (GNN tracker, confirmation by M-of-N or score, JPDA switch), jpda.rs, metrics.rs (OSPA, GOSPA), eval.rs (averaged runs, quality metrics), lib.rs. tests/: ekf_consistency, gnn_scenarios, imm_maneuver, jpda. examples/: 11 programs, one per README table. results/: 7 text files. docs/: roadmap, imm-derivation, jpda.

## E6 Model call sites

No learned models and no language-model or hosted-service calls anywhere. Nothing is downloaded. Randomness: seeded generators from rand and rand_distr (Per the docs; seeds in examples not individually inspected).

## E7 Tests and CI

- 58 tests (Ran it). Test types: unit (inline in src), scenario and statistical integration tests. Coverage not measured.
- Spot checks: not performed on three tests by line; README and jpda.md describe JPDA marginals against hand-computed values, PDA update equal to EKF update to 1e-9, and Hungarian against brute force. Per the docs, not read in the test source.
- CI (.github/workflows/ci.yml): on push to main and pull requests; Ubuntu; stable Rust; steps: cargo fmt check, cargo clippy --all-targets -D warnings, cargo test --all-targets. `gh run list`: all 11 runs (one per commit) completed with success, 24 to 35 s each (checked in the independent review). The README badge links to this workflow. CI does not run any example or compare result files.
- Local run of fmt check and clippy: both clean (Ran it).

## E8 Plan-versus-results and history

All 11 commits between 2026-10-16 and 2026-10-23: v0.1 to v0.4 on 16 October (five commits within about 1 h 15 min on the author-date clock), JPDA milestone and later studies 21 and 23 October. The roadmap states "Each milestone ends with a measured result in the README"; no locked protocol, no success criteria written before results. A commit (1205eaa) removes PROJECT_SUMMARY.md "added by mistake; it is a local file"; the file's content was not inspected (Ask the author whether it held anything private; history still contains it). Commit credits: 11 of 11 carry the same AI co-author trailer; the guide names no tool.

## E9 Gates and approvals

None.

## E10 Environment

Cargo.lock committed; versions pinned by lock. Results do not record versions. A fresh build worked with cargo 1.99.1 on Windows. The CI uses the stable channel, so a future compiler can change warnings under -D warnings (My reading).

## E11 Security and privacy

Secret pattern search over current files (Ran it, `git grep -nIiE "(api[_-]?key|secret|token|password)\s*[=:]"`): no matches. History scan not run. .env is git-ignored; none present. No outside services. GitHub code scanning: "no analysis found". No claim is made about live security state.

## E12 Shortcomings register

1. Simulated data from the filters' own motion model (Gap, High for any real-world claim; disclosed in Limitations). Options: (A) replay recorded or public radar data; (B) add mismatch scenarios (non-constant-velocity, correlated clutter). Recommend B first (cheap), then A. Effort medium; payoff high.
2. Clutter fix thresholds tuned on reported scenarios; score rule given true Pd and clutter density (Gap, Medium, disclosed). Options: (A) hold-out seeds and scenarios; (B) estimate Pd and density online. Recommend A. Effort small; payoff medium.
3. No confidence intervals in the main tables (Gap, Medium, disclosed). Option: report paired intervals as the JPDA studies already do. Effort small; payoff medium.
4. JPDA negative result has untested causes (Gap, Low, disclosed): coalescence-avoiding variants not tried.
5. Experiment outputs for the first four README tables are not saved in results/ (only the JPDA ones are). Reproducibility depends on the example commands (My reading). Effort small.
6. CI does not run examples or compare results (My reading). Payoff medium, effort small.
7. IMM cost not benchmarked; fixed turn rates (disclosed).

## E13 Method detail

- Kalman update toy: prior variance 400 (20 m), measurement variance 225 (15 m); gain 400/(400+225) = 0.64; 5000 + 0.64 x 30 = 5019.2 m; posterior variance (1 - 0.64) x 400 = 144 (12 m). The guide gives "about 5,019" and "more certain than either"; it treats the range as a direct linear measurement (the EKF's linearisation is not shown).
- Assignment toy: costs A1=1, A2=9, B1=2, B2=3; A1+B2 = 4; A2+B1 = 11. Invented toy numbers.
- OSPA toy: (10^2 + 20^2 + 200^2)/3 = 13,500; square root 116.2 m. Matches the metrics.rs form with p = 2, c = 200, n = 3 tracks, m = 2 truths.
- Score-rule formula (README): hit adds ln(Pd N(y;0,S) / (rho x range)); miss adds ln(1 - Pd); confirm at 12 or 20; delete a tentative track below -5.
- Consistency ratios: ideal 4 (state dimension) and 2 (measurement dimension); the guide now describes them in plain words without the acronyms NEES and NIS.

## E14 Discrepancies

- README intro says "(planned) head-to-head benchmark of classical vs. learned data association"; no learned component exists. Consistent with the roadmap (0.6 planned).
- README Limitations says "no confidence intervals yet"; docs/jpda.md reports paired 95% intervals for later JPDA studies. Not a contradiction; later work.
- README single-target table gives RMSE 14.2 versus 23.2 m; the single_target example prints a different table (error versus time). The README table's source is the ekf_consistency test (Per the docs).
- jpda.md "Limits" says "15 to 40 seeds per cell without confidence intervals" for the main studies, while the shallow, merged and IMM studies use paired intervals; both are in the same file.

## E15 Production readiness

Placement: Works locally (builds, passes 58 tests and lint, runs examples); nothing hosted. Checklist: clean start (Yes: lock file, README run commands); checks block a bad change (Partly: tests, fmt, clippy in CI; no result checks); quality against a baseline (Yes in simulation: tuned CV filter and GNN as baselines); safe failure behaviour (Not applicable, no service); human review (Not applicable); prompts and models versioning (Not applicable); logs and monitoring (No); data and secrets (Yes in the narrow sense: synthetic data, no secrets found). Gap to deployable: real data ingestion, multi-sensor timing, performance measurement, packaging, security review. Reusable: the OSPA and GOSPA scoring, the Hungarian solver, and the EKF consistency test are small and portable (My reading). Showcase: possible.

## E16 Review dispositions

Self-review of draft: removed first-person claims of checking; replaced a gain figure with "lower error"; trimmed length from about 4,750 to about 4,350 words by whitespace count (tables and code block inflate it; roughly 4,100 prose words, My estimate); removed two unverified claims about test contents. Reader-check pass (18 edits) applied: terms defined before first use, caveat stated once in the opening plus short pointers, self-praise and notes-to-self removed (moved to the questions file), clutter chart added, parts, IMM and GNN/JPDA charts regenerated. Edit 16 and 18 merged the Results and Strengths limits; the "what numbers support" paragraph now lives in the opening. Not yet reviewed by the author.

## E17 Completeness and backlog

Promised versus built: simulator, radar model, EKF, multi-target scenarios, Hungarian, GNN tracker, IMM, OSPA/GOSPA, JPDA (all built and tested, 0.1 to 0.4 done per roadmap); EO and RF sensors, out-of-sequence handling (0.5, planned); learned association comparison (0.6, planned); Cursor-on-Target output, Python bindings, replay viewer (0.7, planned). Backlog: Finish: none of the built items is unfinished. Fix: save outputs for the first four README tables. Improve: intervals and held-out seeds; mismatch scenarios. Add: coalescence-avoiding JPDA, second sensor, learned association. Next milestone per roadmap: 0.5. Whether to do these is for the author to decide.
