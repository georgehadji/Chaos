# Verification: experiments runner and run store (Amendments D53, D54, D55, D56)

Module-build stage 04 for the second scope of `experiments`. Branch `module/experiments-runner`, code of D55 and D56 (commits d787d25 and 4219a42), verify scripts at fa84391 plus the edits committed with this report. The first verification (code at 9cd8df6) failed the audit on H1 and H2; this report is the repeat after D55, preregistration A5 and D56.

Raw output: `experiments-runner-verify-log.txt` (main run: PASS 80, FAIL 6, INFO 58, 3329 s at 23 threads; addendum: rerun of sections Y5, REALFAIL and DEVICE after three script corrections); device log `gpu-rtx5060-d53-verify-log.txt`; figure `experiments-runner-cost.png`. Scripts: `verify/experiments-runner.jl`, `verify/experiments-runner-reader.jl` (HDF5.jl and SHA only, separate process).

The script had been written for A4. Under A5 it needed five corrections, all in the verification, none in the runner: reduce every `[run]` block, the A5 Budget rows (main and long split into non-L2 and L2 arms), the A5 control scope, the `breakdown_tau` flag in the independent reader, the `session` attribute in the resume comparison. A new section CTRL runs the real control rule at the real horizon.

## Evidence per criterion

| Criterion | Evidence | Result |
|-----------|----------|--------|
| D53a to D53v, D54a to D54e, D55, D56 (unit level) | package suite 4807 of 4807; full `Pkg.test` 54428 of 54428 (at commit 4219a42) | PASS |
| D53a, D53u | wrong thread count, changed preregistration byte, wrong panel hash each refused, nothing written; exit codes 0, 2, 3, 3, 4 | PASS |
| D53b [full], O14 | real panel: 15569 items, 1363 chunks equal an independent count; per-arm steps before and after merge equal; Float64 2.9159e11, Double64 1.308e10, BigFloat 2.76e8 after merge; synthetic capped panel reproduces the 13 Budget rows (A5 split) with relative difference 0; 400 seeded item texts and hashes equal the plan and `sha256sum` | PASS |
| D53w reduced campaign | real panel, horizons 1/100: kill after 272 chunks and during the write of chunk 273, resume equals the uninterrupted campaign (0 differences outside times, host fields and session), injected failure and rerun recorded, rerun on a finished store changes nothing | PASS |
| D53x | `data/p1-le-benchmark/` never created; study manifest folder unchanged | PASS |
| D53m / O13 reader | independent reader on 1363 files per store: 0 failed checks; module reader equals it on 206053 attributes and 188399 datasets; `h5dump` agrees on 36 datasets | PASS |
| D53i threads | `bin/chaoslab.jl run` at `-t 1`, `-t 4`, `-t 23` (745 s, 307 s, 204 s): all datasets identical, manifests equal outside provenance | PASS |
| D53t CPU rerun (O12) | 226 items, bitwise equal | PASS |
| D53f | 630 stored values against the verified estimator's direct call (Float64, Double64, BigFloat, lifted, wrapped): 0 differences | PASS |
| D53s | Y5 at 1e7 on CPU, 8 methods, every grid step: max error / bound 0.0625; other directions 0.000665; fixed-point controls on their orbit 0.0045 | PASS |
| D53w controls (A5 item 3) | off-orbit controls are exactly 878, 880, 881, 882 (INFO); no gated control has a projected bias above 1e-8 | PASS |
| D53k/REALFAIL (addendum) | real horizons, L2 wrapped in main_l2, long_l2, eps_sweep and L4, L5 in tau_sweep: 360 of 1363 chunks, 0 failed twice (limit 68.15) | PASS |
| D53h, D54f device | record at 23 threads: all 8 recorder entries `gpu_faster = false`; `auto` resolves every Float64 method to CPU; forced device arm completes with device provenance | PASS |
| D53q pilot | totals: 60.8 projected thread-hours against 174.1 allowed by the A4 lines with safety; no arm pauses; 11 of 12 arms within a factor 2 of the reduced-campaign cost; tau_sweep 2.46 (the pilot overestimates) | FAIL (marginal, conservative) |
| D53t device agreement | ensemble means within 2.39 SE (bound 4) for all 8 methods; per-orbit worst 5.63 sigma for 7 methods (bound 6), L2 16.2 sigma | FAIL for L2 only |
| D54f recorder retake | 7 of 8 methods within 1.5x of the stored entry on both CPU and device, L2 1.6x to 1.7x on both (twice), L5b CPU 5.6x once (transient); gpu slower than cpu for all 8 in every retake | FAIL (L2 record not reproducible; no decision changes) |
| CTRL (real control rule at N = 1e7) | gated controls 871, 873, 875, 877, 879, Y5: L1, L3, L4, L4b, L5, L5b, L6 at most 4.6e-9; L9 1.8e-8 (877) and 2.3e-8 (879) against the bound 1e-8 | FAIL for L9 |
| D53r | no study value in `run_*.jl` | PASS |

## Findings for the owner

1. **Control bound and L9 (new, needs a decision before any study chunk).** At the real horizon N = 1e7, L9 deviates from the exact exponent by 1.8e-8 on control 877 (K = 2 pi) and 2.3e-8 on control 879 (K = 10); the bound is 1e-8. The deviation falls as 1/n (1.8e-6 at 1e5, 1.8e-7 at 1e6, 1.8e-8 at 1e7), so it is the constant |ln c| = 0.18 and 0.23, not an error that shrinks away. The other seven methods are at most 4.6e-9 on the same orbits. The real campaign stops at the first L9 control chunk (ordinal 192 of the reduced plan). Options: a per-method bound (L9 at 5e-8), a single looser bound (3e-8), or L9 outside the control gate; each is an amendment, allowed until data exist.
2. **L2 recorder record.** The stored L2 entry (CPU 2.94 ns, device 11.19 ns per orbit-step) is not reproduced: retakes give 4.8 and 17.8 ns (1.6x to 1.7x, twice). The gpu stays slower than the CPU by 3.6x, so `auto` still resolves to CPU. A new record at the start of stage 03-run (K19) fixes it.
3. **L2 per-orbit device agreement** 16.2 sigma (heavy-tailed L2 at long horizons). No study effect while `auto` resolves to CPU.
4. **Pilot and budget.** The pilot overestimates tau_sweep by 2.46 (conservative). Reduced-campaign costs projected on the real plan: 59.1 thread-hours against 174.1 with safety (63.6 without); no arm's first run type exceeds the first-chunk cap. Estimated wall time at 23 threads about 2.6 to 4.5 hours.
5. L6 is accepted at Double64 by the plan although L6 has no Double64 support; not exercised by A5.

Closed by D55, A5 and D56 since the first verification: controls 878, 880, 881, 882 (now INFO), the failed-chunk stop at real horizons (L2 wrapped and the tau sweep breakdown flags: 0 failed), the BigFloat budget line (measured rates), the schema notes on axis order and blank rows of failed files (schema v2).

## Independent review

Read only the spec, A4, the schema and the code. No CRITICAL. H1 and H2 and the MEDIUM items M1 to M9 were fixed by D55 (tests in `test/experiments/ex_review.jl`); the LOW items are logged in D55. The verification above exercises the fixed paths: resume after kills, changed backend table, dirty tree, control stop (`bin` exit code 3), manifest reads.

## Audit

| Check | Result |
|-------|--------|
| Evidence | every criterion has a number |
| Review closed | PASS: H1, H2, M1 to M9 closed by D55 |
| Open FAIL rows | CTRL (L9 bound, owner decision), D54f L2 record, D53t L2 device agreement, D53q tau_sweep: none is a defect of the runner |

Outcome: stage 04 complete for the runner. Finding 1 needs an owner decision (amendment A6 or an accepted risk) before stage 03-run of the study.
