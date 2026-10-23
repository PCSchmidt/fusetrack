# FuseTrack roadmap

| Milestone | Deliverable | Status |
| --- | --- | --- |
| 0.1 | Simulator, radar model, single-target EKF, NEES/NIS consistency tests, CI | Done |
| 0.2 | Clutter and missed detections; multi-target scenarios; GNN association with gating; OSPA | Done |
| 0.3 | IMM (CV + coordinated-turn) for a manoeuvring target; write-up deriving the mixing step | Done (single target) |
| 0.3a | Diagnose the heavy-clutter breakdown; score-based track confirmation | Done |
| 0.3b | Integrate the IMM into the multi-target GNN tracker | Done |
| 0.4 | JPDA; OSPA/GOSPA, track fragmentation and ID-switch metrics | Done as a measured negative result: JPDA matches GNN without clutter and is worse in clutter and at close separation (see [jpda.md](jpda.md)); harder shallow-angle and merged-return scenarios gave the same answer; JPDA+IMM built and measured (same answer: no better than GNN+IMM); coalescence-avoiding variants not tried |
| 0.5 | EO bearing-only and RF bearing sensors; per-sensor latency and out-of-sequence measurements | Planned |
| 0.6 | Learned association model vs. classical baselines at rising density and clutter | Planned |
| 0.7 | Cursor-on-Target output; Python bindings (PyO3); replay viewer | Planned |

Each milestone ends with a measured result in the README, not just code.
