# Verification: experiments

Module-build stage 04 for `experiments` (panel scope, D49; amendments D50 withdrawn, D51 in force). Branch `module/experiments`, code at commit 3991155.

## Runs

- Package suite (stage 03): experiments 2597/2597; full `Pkg.test` 52218 pass, 0 fail, 0 error; Aqua, JET, layering clean; coverage of the module 98.9 percent.
- Oracle script: `julia --project=verify -t 23 verify/experiments.jl`, PASS 93, FAIL 0, INFO 60, 795 s. Raw output: `experiments-verify-log.txt`.
- Verify environment: `verify/Manifest.toml` resolved once (Chaoslab gained the `SHA` dependency; six stdlib jll patch bumps); `verify/Project.toml` unchanged.
- Figure: `experiments-portraits.png` (one panel per K, wrapped coordinates, starts by design class, controls, periodic points by class).

## Evidence per criterion

| Criterion | Evidence | Result |
|-----------|----------|--------|
| AC1 to AC18, AC20, AC22, AC24 to AC26, AC28, AC29 | package suite, per the test plan's criterion table | PASS |
| AC2 [full] hash | header 650b3558... equals the module value and an independent `sha256sum` under the SELF rule; a CRLF copy gives the same value; the panel manifest records it | PASS |
| AC12 [full] K = 0 screening | 65536 cells, closed form 1.117e-4, max abs error 4.5e-19 (bound 1.8e-15); T at K = 0: 40 drawn of 53248 qualifying, recount matches | PASS |
| AC19 thread independence | reduced input (grid side 64): panel files at `-t 1`, `-t 4`, `-t 23` byte-identical, SHA-256 f1b16818...; 143 s, 43 s, 29 s | PASS |
| AC21 [full] subsets | S_D 172 (cap 210), S_L 92 (cap 110), nesting holds, no cap applied, recomputation reproduces the file | PASS |
| AC23 cells | every cell populated, short or empty with its reason, or undefined at K = 0 (table below) | PASS |
| AC27 / O6 ChaosTools 3.5.4 | every ChaosTools orbit of orders 1 to 3 at K > 0 matched or accounted for, 0 unmatched, 0 class mismatches, trace agreement within 8.3e-7; at K = 2pi 907 ChaosTools points are D51 ghosts near (pi/2, 0) and (pi/2, pi) (scoring rule: the module counts ill-conditioned rejections at a dividing period and the o-step trace is within 1e-2 of 2); 2 points at K = 0.5 at (pi, pi) lie 1.06e-6 from the module point, just above the 1e-6 match tolerance, on the wrap boundary (INFO) | PASS |
| O7 independent re-evaluation | separate process, TOML and Base only, 5283 periodic rows at 1024 bits: max residual 1.7e-20, windings, classes, exponents, primitive periods agree, trace relative error 1.1e-16; minimum separation 3.14589e-4 reproduced; 10 fixed-point controls and Y5 (ln mu = 0.96242) match closed forms | PASS |
| D51a, D51b, D51d | package suite (bound identity on 200 seeded det-1 matrices within 1e-12) and the full build (no merged class, no `merged` table, no `sn_trace` input) | PASS |
| D51c [full] K = 2pi | build completes; 73963 ill-conditioned rejections, all at K = 2pi (period 1: 44032, 2: 21693, 3: 3475, 4: 3475, 5: 724, 6: 564); no accepted period-1 point within 1e-5 of (pi/2, 0) (nearest 1.571 away); smallest separation 3.146e-4 against 1e-7 | PASS |

## Panel (full study inputs, 23 threads, 401 s)

883 orbits: 870 drawn, 13 controls. Newton: 1966080 seeds, 1359242 accepted, 499944 diverged, 32931 at the cap, 73963 ill-conditioned.

| K | T | R | C | S | X |
|---|---|---|---|---|---|
| 0 | 40 | 40 | undefined | undefined | undefined |
| 0.5 | 40 | 40 | 40 | 25 short | 40 |
| 0.971635406 | 40 | 40 | 40 | 40 | 40 |
| pi/2 | 40 | 40 | 40 | 40 | 40 |
| 2pi | 40 | empty (no elliptic periodic point; accepted by the owner) | 40 | 5 short | 40 |
| 10 | empty (no qualifying grid cell) | empty (no elliptic periodic point) | 40 | empty (no qualifying grid cell) | 40 |

Controls: the fixed points (0, 0) and (pi, 0) at each K > 0, two period-1 accelerator modes at K = 10 (hyperbolic, exponents 2.2697 and 1.7229), and Y5. No ghost control.

## Independent review

A fresh reviewer read only the spec and the code. No CRITICAL or HIGH finding. Open:

| ID | Severity | Finding | Proposed fix |
|----|----------|---------|--------------|
| M1 | MEDIUM | A stored coordinate can be the binary64 2pi (6.283185307179586) where the true value is 0 (256-bit noise of either sign), breaking invariant I1 and making the control order depend on that noise; seen on a period-2 orbit at K = 5 | store x0, p0 through the same key rule (>= 2pi becomes 0); sort controls by the 256-bit key |
| M2 | MEDIUM | Deduplication drops a candidate within delta_dup of a kept orbit without a trace, so two distinct orbits closer than delta_dup merge silently and the E9 refusal never sees them | record the largest merge distance; refuse when it exceeds a small multiple of delta_dup / 10 |
| M3 | MEDIUM | Panel identity holds on one machine and Julia version only: Float64 elementary functions differ across versions and hosts, and a 1-ulp change can flip a cell near a threshold | state the limit in the spec; the manifest sidecar already records the Julia version |
| M4 | MEDIUM | Acceptance checks the residual at the representative only; I7 and AC8 require every point of the orbit | check the residual at every point when the orbit is formed (O7 found max 1.7e-20 on the real panel, so no current case) |
| L1 to L8 | LOW | K = 0 R decimal text not the exact 256-bit rational; negative K accepted; no check that N_scr / 10 is a multiple of tau; SELF rule fixed padding; conditioning at the searched, not primitive, period; StableRNGs version not recorded; path outside the repo not refused; minor items | listed for stage 03 |

## Audit

| Check | Result |
|-------|--------|
| Evidence | every acceptance criterion has a number or a plot (above) |
| Review closed | no open CRITICAL or HIGH finding |

## Notes for the owner

- The study panel and portraits were written to study/02 output by the verify script, but the panel's manifest records `git_dirty = true` (the verify Manifest was modified at build time), and fixing M1 changes stored coordinates. The study panel should be rebuilt after the stage 03 fixes from a clean commit.
- The K = 2pi and K = 10 portraits are dominated by the thousands of hyperbolic periodic points; drawing only elliptic and parabolic points and controls would make them readable.
- Julia 1.10 (CI `lts`) was not run; `Manifest.toml` for 1.10 has a stale `project_hash` until resolved on 1.10.
