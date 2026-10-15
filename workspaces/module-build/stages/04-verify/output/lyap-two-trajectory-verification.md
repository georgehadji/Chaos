# Verification: lyap-two-trajectory

Spec: `../../01-spec/output/lyap-two-trajectory-spec.md` (criteria 1 to 17, Amendment D35). Test plan: `../../02-tests/output/lyap-two-trajectory-test-plan.md` (Amendment D38). Package tests: `test/lyap-two-trajectory/` (`l2_*.jl`). Verify script: `verify/lyap-two-trajectory.jl`, run as `julia +release -t 2 --project=verify verify/lyap-two-trajectory.jl` from the repo root; full output in `lyap-two-trajectory-verify-log.txt` (one line per criterion: `AC<n> | quantity | value | bound | verdict`). Steps 1 and 3 of the stage contract; step 2 (independent review) is the orchestrator's, see "Review findings".

Run facts. Julia 1.13.1, 2 threads, base seed 20261012, verify offset 45000 (seeds `L2_SEED + 45000 + 50 ac + sub`; the package tests use `L2_SEED + 25000 .. 25999`, so all samples are fresh). Library state: working tree of branch `m2-p1-modules` at the time of the run (HEAD a2cc052; the L2 files uncommitted per impl notes). Wall time of the final run 1401 s, of which the package suite 369 s (4 min nominal, loaded machine), the two CairoMakie sections 385 s and 70 s (first-call compilation), the rest under 5 min; other agents were running concurrently, so timings are not benchmarks. The measurement helpers are those of `test/lyap-two-trajectory/l2_helpers.jl` (tolerance formulas `l2_BS`, `l2_D2T`, `l2_d0star`, `l2_valid`, Psi_n via `lb_psi.jl`, orbit and start helpers); the loops are written in the script (the testsets are top-level code and cannot be called), so a definition error in a helper would not be caught here. The independent pieces are the 256-bit closed forms, the L1 estimator and ChaosTools. Values are "worst ratio of error to bound" (pass at <= 1) unless stated. Final verdict tally: 57 PASS, 1 FAIL (AC13), 15 INFO or STAGE04-GPU lines.

## Evidence table

"Where": T = package test (`test/lyap-two-trajectory/`, run at full size by the `suite` section), V = `verify/lyap-two-trajectory.jl`.

| Criterion | Quantity | Value | Bound | Verdict | Where |
|---|---|---|---|---|---|
| suite | package suite at full size (L2_FULL_RUN), assertions | 568 pass, 0 fail, 0 error | all pass | PASS | T via V |
| 1 | linear map at origin, max error / B_S(n), tau in {1,5,10}, n in {20,100,1000}: Float64 / Double64 | 0.02266 / 0.004841 | 1 | PASS | T, V |
| 2 | K = 0 shear, max error / B_S(n), d0 in {2^-20, 1e-9}, tau in {1,10}, n up to 1e4: Float64 / Double64 | 0.06051 / 0.001159 | 1 | PASS | T, V |
| 3 | hyperbolic fixed point, max error / (B_S + 2 n E_nl), K in {1,2,6}, 100 u0 | 3.852e-08 | 1 | PASS | T, V |
| 4 | straddling pair and pair at x = 1: max abs(X2) / (2 E_ro); difference of the two / (2 E_ro) | 0; 0 | 1 | PASS | T, V |
| 4 | plain difference of the straddling pair, ln(norm / d0) (the T3 artifact) | 15.7 | >= 10 | PASS | T, V |
| 4 | `separation` of the pair: norm within 4 E_ro d0 of d0 (ratio); max entry 9.5e-7 <= pi; `separation(displaced(sep))` residual over 4 eps 2 pi | 0; 0 | 1 | PASS | T, V |
| 4 | plain-difference double on the straddling pair: rejected before any work | separation / step 0 / orbit 1 / cpu | same | PASS | T, V |
| 4 | double vs correct hook at x = 1, difference / (2 E_ro) | 0 | 1 | PASS | T, V |
| 4 | mirror case across p = 0, error / (B_S + 2 n E_ro) | 0.04865 | 1 | PASS | T, V |
| 5 | reference orbit differs bitwise from plain loop or L1 stream orbit, 40 cells per type (mode x K x tau x N), Float64 / Double64 | 0 / 0 | 0 | PASS | T, V |
| 6 | pathwise vs L1, max abs(X2 - X) / (Delta_2T Psi_n/n + B_S/n), K = 0.5 / 2 / 5 / 2 pi, d0 = 10^-m, m = 3..12, 0 points outside validity | 0.006326 / 0.02922 / 0.03372 / 0.0116 | 1 | PASS | T, V |
| 6 | old per-step ratio abs(X2 - X) / Delta_2T (reported), same four K | 6.902 / 1.345 / 1.004 / 0.2157 | reported | INFO | V |
| 7 | U-shape, K = 5, 200 starts: slope m in {2,3,4} | 0.9001 | [0.8, 2.2] | PASS (marginal) | T, V |
| 7 | slope m in {10,...,13} | -0.9512 | [-1.3, -0.7] | PASS | T, V |
| 7 | grid minimum at m = 8, d0/d0* (d0* = 5.959e-08) | 0.1678 | [0.01, 100] | PASS | T, V |
| 7 | grid minimum 1.661e-08 over Delta_2T(d0*) = 6.689e-06 | 0.002484 | 1 | PASS | T, V |
| 8 | Double64 T5a ratio, 50 starts, m = 3..28: K = 5 / 2 pi | 0.01376 / 0.01875 | 1 | PASS | T, V |
| 8 | Double64 grid minimum (K = 5), m = 16 | 1.287e-16 | <= 1e-9 | PASS | T, V |
| 8 | Float64 grid minimum / Double64 grid minimum | 1.291e+08 | >= 1e3 (T8 predicts 1.7e7) | PASS | T, V |
| 9 | tau sweep T6a ratio, d0 = 1e-8: K = 0.5 tau 1 / 2 / 4; K = 5 tau 1 / 2 | 0.002393 / 0.0001987 / 2.058e-05; 0.0009609 / 2.778e-05 | 1 | PASS | T, V |
| 9 | reported, not gated: K = 5 tau 4 / tau 8 | 5.747e-08 / 1.472e-12 | reported | INFO | V |
| 9 | control (T7), K = 10, d0 = 1e-3, tau = 16: max X2 | 0.4884 | <= 0.5249 | PASS (marginal) | T, V |
| 9 | control: min (X - X2) over the 50 starts with X >= 1.1 | 1.062 | >= 0.5 | PASS | T, V |
| 10 | (a) d0 = 1e-20 rejected before any work, kind separation, step 0, cpu; (b) collapse at second handoff, orbit 3 / step 2; (c) non-finite in ref / shadow / both: kind state; (d) run record keys; contract ArgumentErrors | all true | true | PASS | T, V |
| 11 | `propagate_shadow` equals serial map of `propagate_shadow_orbit` bitwise, 4096 orbits, chunk 64 / 1000 / default | true, true, true | true | PASS | T, V |
| 11 | outputs bitwise identical for 1, 2 and 4 threads (separate processes, hash of X2, sums, x, W) | true | true | PASS | V |
| 11 | bytes allocated after warm-up, n = 1e4, Float64 and Double64, wrapped and lifted | 0, 0, 0, 0 | 0 | PASS | T, V |
| 12 | K = 5, d0 = 1e-8, T5a ratio: wrapped n = 1e3 / lifted n = 1e3 (100 of 100 valid) | 0.001411 / 0.001392 | 1 | PASS | T, V |
| 12 | wrapped n = 1e5 (100 of 100 valid) | 0.002036 | 1 | PASS | T, V |
| 12 | lifted n = 1e5: 1 valid start tested, 99 with E_ro > 0.5 reported, 0 unexpected errors | 3.163e-05 | 1 | PASS (thin) | T, V |
| 12 | recorded: max abs(X2 - X) at n = 1e5, lifted / wrapped = 1.104e-06 / 1.5e-06 | 0.736 | recorded | INFO | V |
| 13 | ChaosTools oracle, K = 0.5 regular, 20 starts, n = 1e5, d0 = 1e-6, lifted: max abs(X2 - oracle) | 3.181e-04 | <= 1e-4 | **FAIL** | V |
| 14 | oracle ensemble, K = 2 pi, 1000 starts, n = 1e4: abs(diff of means) / combined SE (means 1.167 both, diff 1.66e-06) | 0.0004443 | <= 4 | PASS | V |
| 15 | K = 5 (200 starts): z; K = 10 (200 starts): z | 0.0004639; 0.009381 | <= 4 | PASS | V |
| 16 | K = 2 pi, 1000 starts, n = 1e4, X2 at d0* = 5.683e-08 vs L1 X: z (means 1.169, diff 3.67e-09) | 1.51e-06 | <= 4 | PASS | V |
| 17 | X2 (spec Amendment D42, L2 round-off units), regular K = 0.5 (S7: 1024 of 3072 regular, 17 replaced), lifted, n = 4e4, tau = 1, d0 = 1e-8: max over orbits of \|X2_cpu - X2_gpu\| / (4 E_ro Psi_n / n), all 1024 inside the T5 validity (E_ro max 0.052; Psi_n / n median 7359; max \|X2_cpu - X2_gpu\| 1.2e-3) | 1.104e-04 (bound median 563 against median \|X2\| 7.8e-4) | <= 1 | PASS, **uninformative** (bound > \|X2\| on every orbit; owner question) | `verify/gpu.jl` |
| 17 | X2, d0 = 1e-8, chaotic K = 10, Float64, 1024 orbits, n = 1e4, tau = 2: max per-orbit gap / sigma_FT; difference of means / SE | 4.893; 0.6222 | 6; 4 | PASS | `verify/gpu.jl` |
| 17 | Double64 kernel | does not compile for the device | compile or reason | FALLBACK(kernel does not compile for the element type) | `verify/gpu.jl` |
| R2 | wrapped oracle rule minus lifted rule, 40 starts, n = 1e5, d0 = 1e-6: K = 0.5 / K = 5 | 1.736e-04 / 2.288e-03 | recorded | INFO | V |

Figures: `lyap-two-trajectory-ushape.png` (left: AC7 and AC8, median and max of abs(X2 - X) against d0 for Float64 and Double64 with Delta_2T and d0*; right: AC9, max error against tau with the sigma^(2(tau-1)) Delta_2T line), `lyap-two-trajectory-oracle.png` (left: AC13 X2 against oracle per start; right: AC14 to AC16 ensemble means with one SE).

## Oracle comparison (D30): ChaosTools `lyapunov`

Installed ChaosTools 3.5.4 / DynamicalSystemsBase 3.19.3, source `~/.julia/packages/ChaosTools/0SbsT/src/chaosdetection/lyapunovs/lyapunov.jl`. Signature: `lyapunov(ds, T; u0, d0 = 1e-9, inittest = (u1, d0) -> u1 .+ d0/sqrt(D), Ttr = 0, Δt = 1, d0_upper = 1e3 d0, d0_lower = 1e-3 d0)`; it builds a `ParallelDynamicalSystem` with states `[u0, inittest(u0, d0)]`.

Convention matching.

| Item | ChaosTools | Chaoslab | Matched? |
|---|---|---|---|
| map | user rule; the standard-map rule used is `p += k sin(theta); theta += p` | `StandardMap.forward` | yes; lifted rule, no wrap. 1000 steps at K = 5 give a zero difference of the two orbits (same arithmetic) |
| coordinates | plain Euclidean distance, so lifted coordinates (P4) | lifted `StandardMap` in AC13 to AC15 (minimal image would be a no-op away from boundaries) | yes |
| start direction | `inittest` | seeded unit u_0 of the stream (default of `propagate_shadow`) | yes: `inittest = (u1, d) -> u1 .+ d .* u_0` with the same u_0 per start (orbit id = position) |
| d0 | 1e-6 | 1e-6 | yes (the spec's choice, criteria 13 to 15) |
| realised start separation | `d != d0 && rescale` at start | Algorithm 1 recomputes d_0^+ from stored states | effectively yes |
| transient, time step | `Ttr = 0`, `Δt = 1`, T = n | no transient, tau = 1, N = n | yes; both divide by n |
| rescaling schedule | only when distance leaves [1e-3 d0, 1e3 d0] = [1e-9, 1e-3] | every tau = 1 step | no (by design; hence ensemble comparison in AC14, AC15, and diagnostics for AC13) |
| logarithm reference | ln(d / d0) against the nominal d0 | ln(d / d_{j-1}^+) with the realised separation (T2) | no (placement rounding is in the oracle's sum, not ours) |
| rescale arithmetic | `u2 = u1 + (u2 - u1)/a` | `displaced(ref, delta d0/d)` | same idea |

Results. AC14 (K = 2 pi, 1000 starts): means 1.167 and 1.167, difference 1.66e-06 against a combined SE of 3.7e-03 (z = 0.0004). AC15: K = 5 means 0.9521 both, K = 10 means 1.619 both. AC16 (L1 against L2 at d0*): means 1.169 both. Pathwise (information): the orbits are bitwise the same on both sides and the start direction is the same, so the samples are almost perfectly paired. The combined-SE rule of the spec treats them as independent, hence z of 1e-3 to 1e-2 reflects the pairing, not a strong statistical margin. The informative numbers are the pathwise differences: K = 2 pi median 4.1e-05, max 5.4e-04; K = 5 median 3.9e-05, max 5.6e-04; K = 10 median 3.9e-05, max 6.1e-04; L2 vs L1 at d0* (AC16) median 1.5e-08, max 1.0e-05. The oracle's wider threshold window (up to 1e-3 relative to unit scale, nonlinear) is the likely reason for the 4e-05 pathwise scatter against 1.5e-08 for L2 vs L1; a tight-threshold oracle ([d0/2, 2 d0], information line) gives z = 0.0004 as well.

### AC13 FAIL and its diagnosis (numbers, no change to tolerances)

AC13 gates max abs(X2 - oracle) <= 1e-4 for the 20 regular K = 0.5 starts at n = 1e5, d0 = 1e-6, lifted. Measured 3.181e-04 (mean X2 1.472e-04, mean oracle 1.075e-04; the figure shows one start at X2 = 5.1e-04 against oracle 2.0e-04 and a cluster with X2 about 1.3e-04 above the oracle). The script follows the spec (default thresholds, lifted rule, aligned direction, d0 = 1e-6 both sides). Diagnostics (all in the log, not gates):

- Against the independent L1 tangent value X on the same orbit and direction: max abs(X2 - X) = 3.82e-04, max abs(oracle - X) = 6.40e-05. So the departure is on the Chaoslab side of this pair, in the lifted-mode round-off.
- Measured lifted M_state up to 6.283e+05 (a regular orbit at K = 0.5 reaches coordinate magnitude about n p); E_ro(d0 = 1e-6, M = 6.3e5) = 1.30e-03 per step (worst case); Psi_n/n up to 1.25e+05 (shear-like regular orbits, spec Amendment D35: "grows to the order of n"). The pathwise T5a bound with the measured M_state is then vacuous: the realised ratio is 5.6e-04, i.e. the observed 3.8e-04 is far inside the T5a bound but far outside the per-step-derived 1e-4 of criterion 13.
- Dependence on d0 (max abs(X2 - oracle) / max abs(X2 - X)): d0 = 1e-5: 2.26e-04 / 2.78e-04; d0 = 1e-4: 2.27e-03 / 9.4e-04.
- Oracle with narrower thresholds (max abs(X2 - oracle) / max abs(oracle - X)): [d0/2, 2 d0] 1.29e-04 / 2.57e-04; [d0/10, 10 d0] 1.18e-04 / 2.65e-04; [d0/100, 100 d0] 2.96e-04 / 8.7e-05. A frequent-rescale oracle departs from X by as much as X2 does and agrees with X2 to 1.2e-04, which supports the explanation (placement rounding at |x| about 6e5 amplified by Psi_n/n), not a bug in the oracle side or the schedule.

Reading: the 1e-4 of criterion 13 was derived from the per-step terms (spec: "algorithmic terms are Delta_2T(1e-6) = 1.5e-6") and the oracle bias at d0 = 1e-6; Amendment D35 replaced the pathwise bound of criteria 6, 8, 9, 12 by the Psi_n-amplified form but left criterion 13 unamended. The measured gap is consistent with the amplified form. This is not a failure of AC14 to AC16 (ensemble level, passing). Candidate for sending back to stage 01 (amend criterion 13 to the amplified bound or to an ensemble statement), for the human to decide; no src, test or tolerance was changed.

Recorded line R2 (wrapped rule, plain distance, boundary artifact of T3): K = 0.5 (regular): wrapped 2.784e-04 (SE 4.4e-05) against lifted 1.048e-04 (SE 3.3e-06), artifact +1.74e-04 (2.7 times the lifted value); K = 5: wrapped 0.9766 (SE 8.3e-04) against lifted 0.9743 (SE 7.4e-04), artifact +2.3e-03 (about 2 combined SE; chaotic orbits differ between the two rules after rounding, so only ensemble level).

## Margins and marginal criteria

- AC13: FAIL, see above.
- AC7(a): slope 0.9001 against the floor 0.8 (margin 0.10; the realised slope is near 1, the expected linear E_nl law). AC7(b) -0.9512 in [-1.3, -0.7]. AC7(c): grid minimum at d0 = 1e-8, 0.17 d0*, error 1.66e-08, 400 times below Delta_2T(d0*).
- AC9 control (T7): max X2 = 0.4884 against the cap 0.5249 (margin 7 percent; the cap is a hard bound, 0.488 is the saturation plateau of that d0 and tau).
- AC12 lifted n = 1e5: only 1 of 100 starts is valid (E_ro <= 0.5); 99 reported. The spec (P6) expected 2 of 100 for its start set. Lifted coverage therefore rests on the n = 1e3 line (100 of 100 valid, ratio 0.00139) plus that one start.
- AC8: grid-minimum ratio 1.29e+08 against the T8 prediction 1.7e7 and bound 1e3; Double64 minimum 1.29e-16 at m = 16.
- Realised pathwise ratios are small throughout: T5a ratio at most 0.034 (AC6), 0.0021 (AC12 wrapped, n = 1e5), 0.0010 (AC9 gated). The old per-step ratio abs(X2 - X) / Delta_2T exceeds 1 for K = 0.5 (6.90) and K = 2 (1.345), K = 5 (1.004), the reason for Amendment D35; T5a with Psi_n has margin of at least 29 over all cells here.
- AC11 thread invariance was checked across processes (1, 2, 4 threads), not in-process; AC10 (b) and (c) are the same fixed cases as in the package suite (deterministic maps, no seed dependence).
- AC17 (`verify/gpu.jl`): the first run used the tangent-core count form (tau = 10) and failed, 882 / 876 orbits above 1e-6 (difference 6, bound 3): the CPU itself is above 1e-6 from the Double64 reference on 882 of 1024 regular orbits at d0 = 1e-8 (the round-off of the difference, the AC13 mechanism). Owner's choice: the bound in L2 round-off units of Amendment D42 (tau = 1, where T5a holds directly). It passes, ratio 1.1e-4. Absolute size (CPU recomputation, same starts): bound 4 E_ro Psi_n / n min 2.1e-4, median 563, max 4.8e4; |X2| min 2.1e-5, median 7.8e-4, max 0.014; the bound exceeds |X2| on all 1024 orbits, so in this configuration the row cannot detect a GPU defect. The owner chose the L2 round-off units; the configuration (lifted mode, n = 4e4, tau = 1, factor 4) is the orchestrator's. Owner question at the checkpoint: accept, or rerun where the bound is below |X2| (wrapped mode, smaller n).

## Deviations from impl notes checked

- "Own per-orbit loop, bitwise equal to the serial map for any chunk and thread count" (impl notes Design): confirmed by AC11 (three chunk sizes, 1, 2, 4 threads, 4096 orbits).
- "Reference orbit is bitwise the stream's orbit" (I1): confirmed by AC5, 80 cells, 0 mismatches, both types.
- "0 bytes allocated per orbit after warm-up (AC11, Float64 and Double64, wrapped and lifted)": confirmed, 0 / 0 / 0 / 0.
- "The start check runs on the host for all orbits before any chunk (step 0, backend :cpu)": confirmed by AC10 (a) and the AC4 double (step 0, orbit 1, cpu).
- "d0 must widen exactly into T, finite and positive (ArgumentError)": confirmed in the AC10 contract line (0, negative, NaN, Inf, BigFloat narrowing).
- "Provenance: Tangent's keys with k = 0 plus estimator, d0::Float64, state_mode": confirmed for Float64 lifted, Float64 wrapped, Double64 wrapped.
- Impl notes say "None in behaviour" for deviations from the contract; none found. Not checked here: `state_mode` for systems other than the two standard-map variants.
- Impl note "suite time 127 s on 1.13.1, 194 s on 1.10.12 under load, budget 3 minutes (D35)", re-time item: see Runtime.

## Runtime

Final run (Julia 1.13.1, 2 threads, loaded machine, informative only): suite section 369 s (the package suite itself 4 min nominal at full size: AC1 17 s ... AC12 51 s in the first run), AC1 to 5 29 s, AC6 10 s, AC7 and 8 21 s, AC9 10 s, AC10 and 11 262 s (includes two child Julia processes for the thread check, each with compilation), AC12 59 s, oracle 76 s, figures 385 s and 70 s; total 1401 s. A first full run took about 25 min total on the same machine with the same sections.

Package suite on Julia 1.10: timed once after the final verify run, in a fresh process (`julia +1.10 -t 2`, Julia 1.10.12) with a scratch environment (copy of the prof-env dependency list, Chaoslab by `Pkg.develop`; Julia 1.10 does not read `[sources]`), the package precompile done beforehand and not timed (it took 945 s on the loaded machine). `@testset` around `test/lyap-two-trajectory/runtests.jl` in the package-suite size (no `L2_FULL_RUN`): 308 of 308 pass, 361.99 s (6m01s), compile of the tests included. That is above the 3 minute budget (D35), but the machine was heavily loaded by other agents' concurrent Julia processes (the same package precompile took 945 s here against under 2 min on a quiet machine, and the impl notes measured 194 s on 1.10.12 with two other suites running), so the number is not a verdict on the budget; a quiet re-time is needed before the plan's listed cuts are applied. The 1.10 run also reproduces the sweeps: AC6 T5a ratios 0.0186 / 0.0233 / 0.0146 / 0.0074, AC8 grid-minimum ratio 6.2e+07, AC9 gated ratios at most 0.0017, all inside their bounds.

## Files

- `verify/lyap-two-trajectory.jl` (rewritten; section 1 is the D38 full-size suite; sections AC1-5, AC6, AC7-8, AC9, AC10-11, AC12, oracle, figure, figure-oracle; `L2_ONLY` and `L2_LOG` switches).
- `verify/Manifest.toml`: modified by `Pkg.resolve()` under 1.13 (the stale manifest lacked PrecompileTools, a D39 dependency, so `Chaoslab` failed to precompile in `verify/`). Uncommitted; other verify agents share this file.
- Outputs: `lyap-two-trajectory-verification.md`, `lyap-two-trajectory-verify-log.txt`, `lyap-two-trajectory-ushape.png`, `lyap-two-trajectory-oracle.png`.

## Review findings

Independent review: fresh Sonnet subagent given the spec and the code, not the implementation notes. No CRITICAL and no HIGH finding. Disposition: closed (evidence named), fix (proposed change), record (decision row or spec amendment), owner (checkpoint question).

| Severity | Finding | Location | Disposition |
|---|---|---|---|
| MEDIUM | no per-handoff X2 series (Algorithm 4, L1 S8) | `src/Lyapunov/two_trajectory.jl` | owner: check against D34 "no series sink" and Q9 prefix invariance |
| MEDIUM | provenance stores d0 as Float64; a Double64 d0 cannot be replayed | `src/Lyapunov/two_trajectory.jl:219` | fix: store the exact value (string) |
| MEDIUM | d0 near the rounding floor: realized direction unchecked, E_ro <= 0.5 not enforced | `src/Lyapunov/two_trajectory.jl:33, 67` | owner: guard or record |
| LOW | GPU default literal 100.0 against Tangent `DEFAULT_GPU_NS` | `src/Lyapunov/two_trajectory.jl:174` | fix: use the constant |
| LOW | chunk budget below the real buffers | `src/Lyapunov/two_trajectory.jl:167-170` | fix: count the shadow buffers |
| LOW | default `separation` and `displaced` assume s == d | `src/Lyapunov/separation.jl:9, 17` | record |
| LOW | host-side seeding over all orbits before chunking | `src/Lyapunov/two_trajectory.jl:137-141` | record (memory only) |
| untested | thread count; post-reset window failure code 3; I3 over a run | `test/lyap-two-trajectory` | thread count closed by AC11 (1, 2, 4 threads, child processes, bitwise); the other two: fix, test cases |
