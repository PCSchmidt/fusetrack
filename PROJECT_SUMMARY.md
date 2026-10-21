# FuseTrack Project Summary

## Plain-Language Overview

FuseTrack is a Rust project for studying target tracking and sensor fusion in a realistic but controlled simulation environment. Its current focus is counter-UAS-style tracking: turning noisy radar detections into stable tracks for small moving targets in the presence of missed detections, clutter, and maneuvering motion.

In simple terms, the project asks:

> If a radar reports a messy cloud of measurements every second, which ones are real targets, which ones are false alarms, and how do we maintain a reliable track for each target over time?

FuseTrack builds the classical tracking stack needed to answer that question. It includes a simulator, radar model, Extended Kalman Filter, multi-target association, track confirmation and deletion logic, maneuver handling with an Interacting Multiple Model filter, and quantitative scoring.

The project is not yet a full deployed edge AI system. It is the tracking and evaluation foundation for one. It creates repeatable simulated scenarios and measures how well different tracking approaches perform. The planned longer-term direction is to compare these classical tracking methods against learned or AI-based association methods.

## Project Emphasis

The strongest theme in FuseTrack is disciplined baseline-building. Before adding learned models, the project establishes classical tracking methods that are understandable, testable, and measurable.

That matters because sensor fusion problems can look deceptively simple until clutter, missed detections, target maneuvers, and ambiguous measurements are introduced. FuseTrack emphasizes:

- physically meaningful simulation,
- explicit uncertainty modeling,
- reproducible seeded experiments,
- statistical consistency checks,
- clear failure diagnosis,
- side-by-side benchmark comparisons,
- honest documentation of limitations.

This makes the project especially relevant to physical AI and edge AI work. It is about perception under uncertainty, not just model inference. It shows how a system reasons from imperfect sensor measurements to a coherent world state.

## What FuseTrack Currently Does

FuseTrack currently supports:

- simulating one or more moving targets in 2D,
- generating noisy radar range-bearing measurements,
- randomly dropping real detections according to a detection probability,
- adding random clutter returns,
- tracking single targets with an Extended Kalman Filter,
- tracking multiple targets with a Global Nearest Neighbor tracker,
- assigning measurements to tracks using the Hungarian algorithm,
- confirming tentative tracks and deleting stale tracks,
- using score-based confirmation to reduce false tracks in heavy clutter,
- tracking maneuvering targets with an Interacting Multiple Model filter,
- evaluating multi-target tracking quality with OSPA,
- running reproducible examples and tests.

The current version is best understood as a simulation and benchmark framework for classical radar tracking.

## Repository Layout

The core project files are organized as follows:

- `src/models.rs`: shared state types, radar measurement types, radar model, angle wrapping.
- `src/sim.rs`: target motion simulation, multi-target scenarios, clutter, missed detections.
- `src/ekf.rs`: Extended Kalman Filter for single-target radar tracking.
- `src/assignment.rs`: Hungarian assignment algorithm for one-to-one matching.
- `src/tracker.rs`: multi-target Global Nearest Neighbor tracker and track management.
- `src/imm.rs`: Interacting Multiple Model filter for maneuvering targets.
- `src/metrics.rs`: OSPA multi-target tracking metric.
- `src/eval.rs`: scenario evaluation, benchmark helpers, repeated seeded runs.
- `examples/`: runnable experiments that reproduce README benchmark tables.
- `tests/`: regression and consistency tests.
- `docs/imm-derivation.md`: mathematical derivation of the IMM mixing step.
- `docs/roadmap.md`: planned milestones.

The public library entry point is `src/lib.rs`, which exposes the main modules and commonly used types.

## Simulation Model

FuseTrack simulates targets in a two-dimensional plane. Each target state is:

```text
[x, y, vx, vy]
```

That means:

- `x`, `y`: target position in meters,
- `vx`, `vy`: target velocity in meters per second.

The basic motion model is constant velocity plus random acceleration. This is a common tracking model: the target is assumed to keep moving mostly as it was moving before, while random acceleration accounts for uncertainty and small changes in motion.

The simulator can also create maneuvering trajectories. In those scenarios, a target flies through segments that may be straight or turning at a fixed positive or negative turn rate.

For multi-target scenarios, the simulator creates several targets and then generates radar scans. Each scan contains:

- the true target states,
- noisy radar measurements,
- missed detections,
- clutter measurements,
- measurements shuffled so they are unlabeled.

The tracker does not receive the truth labels. It must infer which measurement belongs to which target.

## Radar Model

The radar model is a 2D range-bearing sensor. For each target, it measures:

```text
range, bearing
```

The model adds Gaussian noise to both range and bearing. This means the radar measurement is close to the true target location but not exact.

The radar model also handles angle wrapping. This is important because angles near `+pi` and `-pi` are actually close to each other, even though their raw numeric values look far apart. The EKF and radar model account for this so bearing errors are interpreted correctly.

At present, FuseTrack only models one radar sensor type. EO/IR, RF bearing-only sensors, asynchronous sensor timing, and out-of-sequence measurements are planned but not implemented yet.

## Extended Kalman Filter

The Extended Kalman Filter, or EKF, is the single-target tracking filter. It estimates a target's position and velocity from noisy radar measurements.

It has two main steps:

1. Predict where the target should be now based on its previous state.
2. Update that prediction using the latest radar measurement.

The EKF is needed because the radar measurement is nonlinear. The target state is Cartesian position and velocity, but the radar measures range and bearing. The EKF linearizes that nonlinear measurement model around the current estimate.

Important implementation details include:

- constant-velocity prediction,
- process noise for unmodeled acceleration,
- nonlinear range-bearing radar update,
- Joseph-form covariance update for numerical stability,
- correct bearing innovation wrapping across the angle branch cut.

The EKF has Monte Carlo consistency tests. These tests check whether the filter's reported uncertainty matches its actual error. The README reports:

- mean NEES close to 4 for a 4D state,
- mean NIS close to 2 for a 2D measurement,
- filtered position RMSE better than raw radar fixes.

This is a strong sign that the EKF is not merely producing plausible-looking tracks; its uncertainty model is also behaving sensibly under the simulated assumptions.

## Multi-Target Tracking

The multi-target tracker is a Global Nearest Neighbor tracker, often abbreviated GNN.

Each radar scan can contain multiple target detections and multiple clutter detections. The key challenge is data association: deciding which measurement belongs to which existing track.

For each scan, the tracker:

1. Predicts every existing track forward.
2. Compares every track to every measurement.
3. Rejects unlikely pairings using a chi-square gate.
4. Builds a cost matrix for plausible track-measurement pairings.
5. Solves the one-to-one assignment problem with the Hungarian algorithm.
6. Updates assigned tracks.
7. Marks unassigned tracks as missed.
8. Starts tentative tracks from unused measurements.
9. Confirms or deletes tracks based on track-management rules.

The tracker is "global" because it chooses the best assignment across all tracks and measurements together, rather than greedily assigning one track at a time.

## Assignment Algorithm

The assignment module implements the Hungarian algorithm. This solves the problem:

> Given many tracks and many measurements, choose the lowest-cost one-to-one pairing.

Costs are based on the measurement innovation. A measurement is cheaper for a track if it is close to what that track predicted. Pairings outside the statistical gate are marked forbidden.

The assignment implementation is tested against brute-force search on small random matrices, which gives confidence that it finds the true optimum.

## Track Management

Tracking is not only about updating existing tracks. The tracker also needs to decide when a possible target becomes a real track and when an old track should be removed.

FuseTrack uses tentative and confirmed tracks:

- A tentative track is a possible target that has not yet earned enough evidence.
- A confirmed track is trusted enough to report as a real target.

The default rule confirms a tentative track after 3 hits in 5 scans. Confirmed tracks are deleted after too many consecutive misses.

This simple rule works well in clean or moderate clutter. However, the project shows that it can fail in heavy clutter because random clutter detections can accidentally form track-like sequences.

## Heavy-Clutter Diagnosis

One of the most valuable parts of FuseTrack is its heavy-clutter analysis.

The default tracker works well with no clutter, light clutter, and moderate clutter. But at very high clutter, such as 200 clutter returns per scan, it confirms too many false tracks. The issue is not that real targets are missed. The issue is that false tracks are being created from clutter.

The README explains the mechanism clearly:

- every unused clutter return can spawn a tentative track,
- heavy clutter creates many tentative tracks,
- some tentative tracks receive enough accidental follow-up clutter hits,
- those tracks satisfy the default 3-of-5 confirmation rule,
- the tracker ends up reporting extra confirmed tracks.

FuseTrack then compares fixes. A stricter M-of-N rule helps, but it can delay confirmation or miss real targets when detection probability is lower. The better option is score-based confirmation.

## Score-Based Confirmation

The score-based confirmation rule uses a likelihood-ratio idea. Instead of simply counting hits, it asks whether the evidence looks more like a real target or more like clutter.

On each scan:

- a hit increases or decreases the score depending on how likely the measurement is under the track model compared with clutter,
- a miss penalizes the score based on detection probability,
- a tentative track is confirmed when its score gets high enough,
- a tentative track is deleted when its score falls too low.

This rule uses assumptions about detection probability and clutter density. In the current benchmarks, those values are known from the simulator. In a real system, they would need to be estimated or configured.

The README results show that score-based confirmation dramatically improves heavy-clutter behavior. It cuts false confirmations while avoiding the long confirmation delays that can come from overly strict hit-count rules.

## Interacting Multiple Model Filter

FuseTrack includes an Interacting Multiple Model filter, or IMM, for maneuvering targets.

A normal constant-velocity EKF can struggle when a target turns sharply. Increasing process noise can help it follow turns, but that also makes it noisier during straight flight.

The IMM handles this by running several motion models in parallel:

- constant velocity,
- coordinated turn at a positive fixed turn rate,
- coordinated turn at a negative fixed turn rate.

Each model has a probability. As measurements arrive, the IMM updates those probabilities. If the target is flying straight, the constant-velocity model should dominate. If the target turns, the matching turn model should gain probability.

The IMM then combines the models into one estimated state and covariance.

The derivation in `docs/imm-derivation.md` explains the key approximation: the exact model-history mixture grows exponentially, so the IMM repeatedly merges mixtures by matching their mean and covariance. This keeps computation bounded while preserving useful maneuver information.

## IMM Results

The README reports two IMM benchmark types:

1. Single-target maneuver tracking.
2. IMM used inside the multi-target GNN tracker.

The results show a practical and nuanced conclusion:

- A low-noise constant-velocity EKF performs poorly during turns.
- A high-noise constant-velocity EKF is a stronger baseline because it can tolerate maneuvers.
- The IMM performs best when its assumed turn rates are close to the true turn rates.
- If the true turn rate is far from the IMM's fixed turn models, the IMM advantage shrinks or disappears.

This is a good engineering result. The project does not claim the IMM is universally better. It shows when it helps and when model mismatch limits it.

## Evaluation Metrics

FuseTrack uses OSPA as its main multi-target metric. OSPA stands for Optimal Subpattern Assignment.

In plain language, OSPA measures how different the set of estimated tracks is from the set of true targets. It penalizes:

- position error,
- false tracks,
- missed targets.

This is better than only measuring position RMSE because multi-target trackers can fail by creating extra tracks or dropping real targets. OSPA captures those failures.

The evaluation module also reports:

- mean number of confirmed tracks,
- mean number of true targets,
- mean false tracks,
- mean missed targets,
- mean confirmation delay.

These additional measures help diagnose why a tracker's OSPA score changed.

## Examples

The `examples` directory contains runnable experiments:

- `single_target.rs`: tracks one simulated target with an EKF.
- `clutter_sweep.rs`: measures GNN tracker performance across clutter and detection probabilities.
- `clutter_diagnose.rs`: compares track confirmation rules in heavy clutter.
- `imm_maneuver.rs`: compares IMM against constant-velocity EKFs on a maneuvering target.
- `imm_tracker.rs`: runs the multi-target tracker with IMM-based tracks.

These examples are more than demos. They are the source of the README's benchmark tables and project claims.

## Tests And Verification

The project has a useful test suite covering the main algorithms and behaviors.

The unit tests cover:

- radar geometry,
- angle wrapping,
- simulation reproducibility,
- constant-velocity motion,
- maneuvering motion,
- EKF prediction and update behavior,
- covariance symmetry,
- bearing branch-cut handling,
- assignment correctness,
- OSPA behavior,
- IMM model probabilities,
- IMM model preference during straight and turning flight,
- tracker confirmation and deletion behavior,
- score-based track management,
- IMM tracking inside the GNN tracker.

The integration tests cover:

- EKF statistical consistency,
- EKF improvement over raw measurements,
- clean multi-target tracking,
- light clutter and missed detections,
- degradation under increased clutter,
- score-based confirmation in heavy clutter,
- IMM performance on maneuvering targets.

Local verification was run with:

```bash
cargo test
```

The result was:

```text
41 tests passed
0 failed
```

Breakdown:

- 32 unit tests passed,
- 2 EKF consistency tests passed,
- 5 GNN scenario tests passed,
- 2 IMM maneuver tests passed.

This means the current codebase supports the high-level project claims.

## Current Strengths

FuseTrack has several strong qualities:

- It is built around clear physical assumptions.
- It uses reproducible seeded simulations.
- It includes meaningful statistical checks, not just visual demos.
- It exposes important tracker failure modes.
- It compares fixes quantitatively.
- It documents benchmark results in the README.
- It keeps classical baselines explicit before adding learned methods.
- It has a compact, readable Rust implementation.
- It has tests that cover both algorithmic pieces and end-to-end behavior.

The project is especially strong as a portfolio demonstration of physical AI engineering judgment: uncertainty, dynamics, association, clutter, and evaluation all matter.

## Current Limitations

FuseTrack is still a simulation-only project. The most important limitations are:

- only one radar sensor is modeled,
- no real sensor data is used,
- no EO/IR or RF sensors are implemented yet,
- no asynchronous sensor fusion is implemented yet,
- no out-of-sequence measurement handling is implemented yet,
- no JPDA association is implemented yet,
- no learned association model is implemented yet,
- target births and deaths are simplified,
- clutter is uniform and independent between scans,
- score-based confirmation assumes known detection probability and clutter density,
- IMM turn rates are fixed rather than estimated,
- benchmark results currently have modest seed counts,
- confidence intervals are not reported yet,
- there is no visualization or replay viewer yet,
- Python bindings are planned but not implemented.

These limitations do not weaken the current project; they define its current scope. FuseTrack is best described as a classical tracking benchmark foundation, not yet a complete deployed sensor-fusion product.

## Roadmap Direction

The roadmap points toward a broader sensor-fusion and edge-AI benchmark platform.

Planned areas include:

- JPDA association,
- GOSPA and additional track quality metrics,
- ID-switch and fragmentation metrics,
- EO bearing-only sensors,
- RF bearing sensors,
- per-sensor latency,
- out-of-sequence measurements,
- learned association models,
- Cursor-on-Target output,
- Python bindings,
- replay visualization.

The learned-association milestone is especially important for the physical AI theme. FuseTrack is positioning itself to answer:

> When does a learned association model outperform classical GNN or JPDA baselines, and under what density, clutter, and maneuver conditions?

That is a stronger framing than simply adding AI for its own sake.

## Overall Assessment

FuseTrack is a well-scoped and technically credible tracking project. It does not try to be everything at once. Instead, it builds the foundation pieces that a real sensor-fusion system needs:

- motion models,
- measurement models,
- filtering,
- uncertainty,
- association,
- track management,
- clutter handling,
- maneuver handling,
- scoring,
- reproducible evaluation.

The project's current emphasis is classical radar tracking under realistic simulation stress. Its future emphasis is comparing those classical methods against learned association models for edge/physical AI use cases.

The most important takeaway is that FuseTrack is not just a code exercise. It is a measurement-driven investigation of how tracking systems behave when the sensor picture is uncertain, incomplete, and cluttered. That makes it a strong foundation for later multi-sensor fusion and learned-association work.

