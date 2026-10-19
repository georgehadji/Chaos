# Verification: experiments runner and run store (Amendments D53, D54)

Module-build stage 04 for the second scope of `experiments`. Branch `module/experiments-runner`, code at commit 9cd8df6. Raw output: `experiments-runner-verify-log.txt` (main run PASS 73, FAIL 4, INFO 53, 3240 s at 23 threads, plus addendum sections); device log `gpu-rtx5060-d53-verify-log.txt`; figure `experiments-runner-cost.png`. Scripts: `verify/experiments-runner.jl`, `verify/experiments-runner-reader.jl` (HDF5.jl and SHA only, separate process). Verify environment: `HDF5` and `HDF5_jll` added as direct dependencies.

## Evidence per criterion

| Criterion | Evidence | Result |
|-----------|----------|--------|
| D53a to D53v, D54a to D54e (unit level) | package suite 4649 of 4649; full `Pkg.test` 54270 of 54270; Aqua, JET, layering clean | PASS |
| D53a, D53u | wrong thread count, changed preregistration byte, wrong panel hash each refused, nothing written; exit codes 0, 2, 3, 3, 4 | PASS |
| D53b [full] | real panel: per-arm items and steps before and after merge equal an independent count (Float64 2.9331e11 after merge, Double64 1.308e10, BigFloat 2.76e8, merge saves 1.032e10); synthetic capped panel reproduces all 11 Budget rows exactly | PASS |
| O14 plan | independent count 15741 items, 1387 chunks; 400 seeded item texts and hashes equal the plan and `sha256sum` | PASS |
| D53w reduced campaign | real panel, horizons 1/100: kill after 277 chunks and during the write of chunk 278, resume equals the uninterrupted campaign (0 differences outside provenance), injected failure and rerun recorded, rerun on a finished store changes nothing | PASS as mechanics; FAIL as a study run (findings 1 and 2) |
| D53x | `data/p1-le-benchmark/` never created; study manifest folder unchanged | PASS |
| D53m / O13 reader | independent reader on 6 stores x 1387 files: 0 failed checks; module reader equals it on 204228 attributes and 187688 datasets; `h5dump` agrees on 36 datasets | PASS (finding 5 on axis order) |
| D53i threads | `bin/chaoslab.jl run` at `-t 1`, `-t 4`, `-t 23` (898 s, 298 s, 191 s): all 110359 datasets identical, combined SHA-256 e320fad1... | PASS |
| D53t CPU rerun (O12) | 186 items, bitwise equal | PASS |
| D53f | 630 stored values against the verified estimator's direct call over Float64, Double64, BigFloat, lifted and wrapped: 0 differences | PASS |
| D53s | Y5 at 1e7 on CPU, 8 methods, every grid step: max error / bound 0.0625; other directions follow ln mu + ln(c)/n | PASS |
| D53h, D54f device | record at 23 threads, all 8 recorder entries `gpu_faster = false`; `auto` resolves every Float64 method to CPU; forced device arm completes with device provenance | PASS |
| D53t device agreement | ensemble means within 2.39 SE (bound 4) for all 8 methods; per-orbit worst 5.63 sigma for 7 methods (bound 6), L2 16.2 sigma (finding 3) | FAIL for L2 only |
| D53q pilot | 10 of 11 arms within factor 2 of the campaign cost; tau_sweep 2.19 | FAIL (marginal) |
| D53r | no study value in `run_*.jl` | PASS |

## Findings for the owner (study inputs and criteria)

1. Controls stop the campaign. Controls 878, 880, 881, 882 are periodic points whose binary64 start leaves the orbit; their lambda_1 at 1e5 differs from the analytic value by 0.10 to 0.65, so the control bound (1e-8) stops the campaign at the first control chunk. The other 6 controls and Y5 pass (projected bias at most 5.3e-9).
2. The failed-chunk stop rule fires at real horizons. L2 (main, long, eps_sweep) and L4, L5 at tau 8 and 16 alone give 103 chunks failed twice at real horizons (limit 0.05 x 1387 = 69): 85 L2 separation errors (mostly eps 1e-10, 1e-12 and the long arm) and 18 orthonormalisation breakdowns at large tau. Cause of the L2 errors not established.
3. L2 per-orbit device agreement 16.2 sigma: heavy-tailed L2 at long horizons (CPU-only perturbation spread 0.028 against L1 0.0047). No study effect while `auto` resolves to CPU.
4. Budget lines: at 23 threads only `ladder_bigfloat` pauses at its first group (193638 thread-ns against the 109200 cap); tau_sweep and long sit at 0.90 of the cap. Thread cost per step grows about 1.7x from 4 to 23 workers. Campaign estimate 58.7 thread-h, about 2.6 to 4.5 h wall on 23 threads.
5. Schema axis order: HDF5.jl stores Julia (n, k) arrays with file extents (k, n); `schema/run-store.md` says "n x k".
6. Failed chunk files blank rows after the failure point (t = 0, NaN); the schema does not say so.
7. The plan accepts L6 at Double64 although L6 has no Double64 support; not exercised by A4.

## Independent review

Read only the spec, A4, the schema and the code. No CRITICAL.

| ID | Severity | Finding |
|----|----------|---------|
| H1 | HIGH | A control-bound stop is undone by a plain resume: the chunk is stored as a final failed file and skipped, the stop reason is lost |
| H2 | HIGH | A resumed run with a different resolved backend table is not refused (spec R4 item 3); GPU and CPU chunks mix silently |
| M1 | MEDIUM | `mv(...; force = true)` is not atomic on Windows (falls back to delete and copy); manifest write can throw on a concurrent reader |
| M2 | MEDIUM | No fsync before rename; root attributes not in the read-back check |
| M3 | MEDIUM | Wall cap under-counts killed sessions (`s.now` stale at manifest writes) |
| M4 | MEDIUM | Missing or corrupt manifest silently resets wall time, pauses and overrides |
| M5 | MEDIUM | `code_version` hashes every file under `src/` (CONTEXT.md, coverage files) and omits `ext/` and `Project.toml` |
| M6 | MEDIUM | Panel-derived values (design class, control exponent, setup class) are in no hash; resume never compares `panel_sha256` |
| M7 | MEDIUM | `threads` argument trusted instead of `Threads.nthreads()` |
| M8 | MEDIUM | Cost-cap pause reaches disk only after the first chunk file; a kill in between loses it |
| M9 | MEDIUM | Chunk files do not record a dirty-tree session |
| L | LOW | control rule on attempt 1 only; budget line from the first precision of an arm; `csize = 0` not mapped to out-of-memory; launch budget calibrated on the AC13 chain; small TOCTOU on `sha_start`; `rerun_attempts` not in the hash; no store lock; stop check after a group; default manifest path not git-ignored |

## Audit

| Check | Result |
|-------|--------|
| Evidence | every criterion has a number |
| Review closed | FAIL: H1 and H2 open |

Outcome: back to stage 03 for H1, H2 and the MEDIUM items; findings 1, 2 and 4 need owner decisions before any study chunk (amendments are allowed until data exist).
