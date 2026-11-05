# Implementation notes: gali (catalog I3, GALI_k for k = 2..K)

Code (new): `src/Indicators/gali.jl` (189 lines), `bench/gali.jl` (28). Edited: `src/Indicators/Indicators.jl` (export list, include, module docstring), `src/Indicators/CONTEXT.md` (80 lines). Spec `gali-spec.md` with Amendment D107, test plan `gali-test-plan.md`, tests `test/gali/` (unchanged).

## Status

| Run | Result |
|-----|--------|
| RED (stage 02) | 5 pass, 46 error (all `UndefVarError`) |
| gali suite alone, `-t2`, testenv of stage 02, driver `ga02/myred.jl` | 2704 pass, 0 fail, 0 error, 2 min 26 s (51 testsets); the first draft (unrolled `SMatrix` algebra) passed 2696 and failed the 8 allocation tests of AC17c, then took 25 min because Double64 at d = 200 compiled for 901 s: rewritten, see Design |
| Whole suite `julia -t 4 --project -e "using Pkg; Pkg.test()"`, Julia 1.13.2 | 1 408 257 pass, 0 fail, 0 error, 1 broken (pre-existing, weighted-birkhoff), 59 min 24 s in the testset (60 min 47 s wall; the machine was not idle) |

Suites in the whole run: quality 12 (Aqua and JET), layering 2, analysis 2763, coupled-maps 643, diagnostics 245, experiments 9269, fli-ofli 263533, gali 2704, lyap-benettin 15009, lyap-gram-schmidt 18938, lyap-householder-qr 575, lyap-svd-stabilized 12915, lyap-two-trajectory 493, megno 1000962, naff-fma 900, precision 150, rli 64460, sali 9210, standard-map 474, stretching 3476, tangent-core 808, weighted-birkhoff 716 (1 broken). The gali suite inside `Pkg.test` took 4 min 58 s (4 threads, machine loaded by the run), against the 3 minute ceiling of D38 that the 2 min 26 s of the quiet run meets. The whole suite ran before the last one-line edit (the constant `ArgumentError` text in `reducer_init`, checked by hand); the gali suite ran after the algorithm was final.

## Change table

| File | Lines | Change |
|------|-------|--------|
| `src/Indicators/gali.jl` | +189 new | `GALI`, `GALIState`, `reducer_init/step/update/failed/error/final/provenance`, `gali_values`, `gali_value`, `gali_floor`, `check_frame(::GALI, ...)` |
| `src/Indicators/Indicators.jl` | +3 -2 | exports `GALI, gali_values, gali_value, gali_floor`; `include("gali.jl")` after `alignment.jl`; module docstring |
| `src/Indicators/CONTEXT.md` | +11 -10 | gali in the runs line, labels, files and exports tables, limits; 80 lines (the `check_frame(SALI)` row was merged into the `alignment_indices` row, the deferred line lost GALI_k) |
| `bench/gali.jl` | +28 new | step and `gali_values` on Y3 and Y4 N = 100 (K = 4, both types), log-form step, `check_frame`, a 1000 step Y3 orbit |

No test file, no other source file changed (`git diff -- test` empty).

## Design

- State `GALIState{T, D, K, M}` (`M = K - 1`): `SVector{K}` sums and compensation terms, `SVector{M}` current values, `value`, `n`, `step`, `code::Int8`, `column`. Isbits. Per column `add(CompensatedSum, log(alpha_i))`, the same operations as the sali reducer, so sums and `W` are bitwise those of sali (AC8b, AC8c pass with equality).
- Failure: first column with a zero or non-finite scaled norm; state and `W` unchanged, absorbing (as sali). `reducer_error` returns kind "deviation matrix", threshold 0 for a zero norm and Inf otherwise.
- Evaluation: `_gram_inplace!`, a modified Gram-Schmidt of the unit columns in place, each column orthogonalised twice against the earlier ones; `G_k = prod_{i = 2..k} abs(R_ii)` for every k at once (spec B2 iv). Column j sees only columns 1..j-1, so G_2, G_3 of a K = 4 run are bitwise those of K = 2, 3 runs (AC8a, reported 63 of 63 bitwise by the plan; passes). A zero column stays zero (divide by 1), so every later order is exactly 0.
- Log form (K = d): `L = -(((s_1 + ... + s_d) + (c_1 + ... + c_d)))`, plain left-to-right sums of the sums and of the compensation terms; at d = 2 this is the sali expression `-((s1 + s2) + (c1 + c2))` to the bit (AC8c). 0 at n = 0, `nothing` for K < d.
- Floor `gali_floor(T, k) = 2 k effective_eps(T)`; `reducer_final.floor` is the `SVector{K-1}` of the floors of the orders 2..K.
- Work matrix: `MMatrix` (stack, zero allocation) for isbits element types; a heap `Matrix` for BigFloat (CPU only). The scaled column norm `_norm_col` is `safe_norm` written as loops with the same operation order (checked bitwise on 800 random vectors of d = 2, 4, 9, 200, entries over 6 decades).
- Why loops and `@inbounds` (three failed attempts first): the unrolled `SMatrix`/`SVector` algebra at d = 200 compiled 8 s in Float64 and 901 s in Double64 (an 800-element unrolled Double64 expression); tuple-of-columns recursion allocated 208 B per step (a closure over a reassigned variable boxed it); and `@inbounds` is needed because a bounds-error path that holds the `MMatrix` defeats the stack allocation (6.4 kB per step at d = 200). A ternary inside the inner scaling loop did the same, hence `s = iszero(r) ? one(T) : r` outside it. These are Julia 1.13 compiler behaviours, not documented anywhere; the comments in the file say why.
- No string interpolation in anything the stream kernel reaches: `reducer_step`, `_gram_inplace!`, `_norm_col`, `reducer_final` hold none; the `ArgumentError` of `reducer_init` is a constant text. `gali_values` (K < 2), `gali_value` and `check_frame` interpolate in host-only error messages.

## Performance (`bench/gali.jl`, `verify` environment, `-t2`, BenchmarkTools medians, indicative, machine not idle)

| Entry | Float64 | Double64 | Allocations |
|-------|---------|----------|-------------|
| Y3 (d = 4, K = 4) `reducer_step` | 0.2 us | 2.4 us | 0 |
| Y3 `gali_values` | 0.2 us | 2.1 us | 0 |
| Y3 `reducer_step`, K = d = 4 (log form) | 0.2 us | n/a | 0 |
| Y4 N = 100 (d = 200, K = 4) `reducer_step` | 5.5 us | 64.1 us | 0 |
| Y4 N = 100 `gali_values` | 5.1 us | 64.2 us | 0 |
| `check_frame` (K = 4) | 0.2 us | n/a | 0 |
| Y3 `propagate_orbit`, 1000 steps, K = 4 | 269 us | n/a | 0 |

At d = 200, K = 4: Float64 5.3 us per step by a plain timing loop (Double64 69.6 us), against 21.8 us (Float64) and 204 us (Double64 not measured to completion earlier: its first call alone was 901 s) for the first unrolled draft. K = d = 20: 7.1 us per step. The Float64 step of Y4 N = 100 costs about as much as the banded frame action (1.2 us on 4 vectors, coupled-maps notes) times 4.

First-call compile: `reducer_step` alone, d = 4, 10, 20, 200: 0.07 to 0.14 s each (loops, not unrolled). `propagate_orbit` on Y4 N = 100 with K = 4: 40.7 s (Float64) and 42.1 s (Double64) for the first call, of which GALI adds 0.13 s: the same call with no reducer takes 43.0 s, SALI (K = 2) 4.1 s. The 40 s is the fixed-size chain state (D105) multiplied by a 200 x 4 frame (`ChainTangent * SMatrix`), not this module; coupled-maps recorded 3.1 s for a narrower frame. Record for stage 04 (criterion 20 time ratio) and for the p3 study (it pays once per process and per K).

Not measured here: GPU compile of the reducer (criterion 21, stage 04). Risk: the `MMatrix` buffers (6.4 kB at d = 200, K = 4) inside a CUDA kernel need the stack allocation to succeed in the GPU compiler as well; if not, the Y4 N = 100 chain falls back to the CPU with a recorded reason, as the spec allows (Compute). Y3 (d = 4) is 128 bytes of buffer.

## Audit

| Check | Result |
|-------|--------|
| GREEN | gali 2704 of 2704; whole suite 1 408 257 pass, 0 fail, 0 error, 1 broken (pre-existing) |
| Quality | Aqua and JET clean (quality 12 of 12); layering 2 of 2. Coverage by `--code-coverage=user` on the gali suite alone (stage 02 test environment, driver `ga02/myred.jl`): `gali.jl` 90 of 90 executable lines (100 percent). 456 `.cov` files written under `src/`, `test/` and `bench/` were all deleted (none tracked; `git status` shows only the four intended paths) |
| Size | file 189 lines; longest function `_gram_inplace!` (22 lines); `reducer_step` 12; every other function under 12 |
| Tests untouched | `git diff --stat -- test` empty |
| Numerics | no `@fastmath`; per-column log sums use `CompensatedSum` and `add`; the log form sums the sums and the compensation terms separately, then adds them (bitwise sali at d = 2) |
| Em dashes | none in the new or edited files |

Realised ratios from the gali run (printed by the tests, indicative): AC1 worst error / (2 B_meas) 0.20 (Float64) and 0.0097 (Double64); AC13 orbit D l_1 = 0.00765, l_1 - l_2 = 0.00676 at n = 3000; AC14 12 of 18 starts in band (1 excluded as lower-dimensional); AC15a median floor ratio Double64 / Float64 1.66 (k = 3), 1.83 (k = 4); AC16e 41 of 100 starts reach bitwise equal columns; AC18b worst ratio 0.0029; AC6 m_floor 17 (Float64), sharper bound violated 0 of 160. The log-form ratio against 32 eps sum Psi and the AC9 error ratios are in `scratchpad/ga03/green4.out`, not in this file (stage 04 reports them at full size).

## Readings and deviations

- Evaluation route: modified Gram-Schmidt twice (spec P2 allows "Householder or MGS2"); no Householder. Accuracy per AC1: 0.20 of the bound in Float64.
- Log form for K = d > 2 is a plain left-to-right sum of the K sums plus a plain sum of the K compensation terms, not a Neumaier total of the sums. The spec text says "negative of the left-to-right compensated total"; bitwise equality with sali at d = 2 fixes the form there (a Neumaier total would carry the rounding of `s1 + s2` and differ). For d > 2 the extra error is a few eps of the total, inside B_L = 32 eps sum Psi (AC9, AC11, AC20 pass); a Neumaier form for d > 2 would be a separate code path and is not done.
- `_norm_col` duplicates the arithmetic of `safe_norm` (loops, same operation order) instead of calling it, for compile time at d = 200; the arithmetic is verified bitwise equal. If `safe_norm` changes, `_norm_col` and the sali sums-equality tests (AC8b) will show it.
- `gali_value(W, k)` copies `W[:, 1:k]` to a `Matrix` (dynamic `k`); the pure functions on `Matrix` input convert to `SMatrix` (type-unstable, host use only).
- The frame check also fails for a NaN value (the comparison `g0 >= threshold` is false).
- `check_frame` accepts any `AbstractMatrix` with at least 2 columns; it does not check `K <= d`.

## Process deviations

- The first shell command that wrote the file failed to parse (a heredoc quote), so the first draft was written with the editor tool. No effect on content.
- Four versions of the evaluation were written and measured (see Design); only the last is in the tree.
- The gali suite was run in the stage 02 test environment (`-t 2`) as the instruction gives, and again inside `Pkg.test` (`-t 4`); coverage was measured on the first only, not inside `Pkg.test`.
- Unrelated count: `experiments` is 9269 testsets-assertions here against 9267 in the coupled-maps notes; nothing in this module touches Experiments; not investigated.

## Known limits and not done here

- Criteria 19 (ChaosTools k = 3, 4), 20 timing, 21 (GPU), the 1 against N threads comparison, the full-size runs (`GA_FULL_RUN`) and the realised constants of P5 are stage 04 (`verify/gali.jl`).
- Not edited (outside this stage): `src/precompile.jl`, `docs/src/api/*.md`, `CHANGELOG.md`, `_config/catalogs/indicators.md` (rows I3, I4 and the Status column, stage 05), `_config/architecture.md`.
- The Experiments layer must call `check_frame(GALI(), W0; orbit)` for a supplied frame and put `GALI()` last in the chain (D71, D107e); nothing in this stage wires GALI into `run_indicator_record`.
- First-call compile of the Y4 N = 100, K = 4 chain frame product is 40 s (Design); the p3 study should pay it once per process.

## Conflicts between spec, docs and code

- Spec Algorithm item 3 says "uses no library call the GPU kernel cannot compile"; the `MMatrix` stack buffer is a StaticArrays type, not a library factorisation, but whether the CUDA compiler keeps it on the stack is untested (stage 04, criterion 21).
- Spec B9 and the API table give `gali_floor(T, k) = 2 k effective_eps(T)`; the code returns `2 * k * effective_eps(T)` (same value, same operation order as the test constant).
- Test plan "API" says `reducer_final` is `(n, gali, L, sum, comp, floor)`; implemented as written. The plan notes that `check_frame` is not called by `propagate`; the code follows that, spec Algorithm item 1 still reads as if the run rejects the frame.
- `src/Indicators/CONTEXT.md` listed GALI_k as deferred to M5 (D70); updated. `_config/catalogs` still says I3 is planned (stage 05).

## Rework D109 b to d (2026-11-05)

Only `src/Indicators/gali.jl` changed (193 to 217 lines; longest function unchanged in order of 20 lines). Tests first, see the test plan.

- D109b: `gali_values(::SMatrix)` throws `ArgumentError("GALI needs d >= k")` (constant text) before any work; the `AbstractMatrix` method and `gali_value` (which checks `size(W, 1) >= size(W, 2)` itself) interpolate the sizes; `check_frame` gets it through `gali_values`.
- D109c: the loops moved into `_gram_core!(Q, D, K, G)` and `_norm_col(M, j, D::Int)` with runtime sizes; the `Val` methods call them (same code, same operation order, so the reducer hot path and the zero-allocation result are unchanged). `G` is an `SVector` (`Base.setindex`) on the static path and a `Vector` (in place) on the host path (`_put`). `gali_values(::AbstractMatrix)` uses a heap `Matrix` buffer and returns a `Vector{T}` (was an `SVector` from a dynamic-size `SMatrix`); inferred concretely, bitwise equal to the `SMatrix` path.
- D109d: `float(T) === T || return gali_values(float.(W))` in both methods; `check_frame` uses `gali_floor(float(T), K)` (so `effective_eps(float(T))`). `reducer_init(GALI(), Int, Val(3), Val(2))` does not throw today and still does not (the state is built with `zeros(SVector{K, Int})`; the reducer is only reached through `propagate`, whose element type is a float); not changed.

First call at d = 200, K = 4 (fresh process, `t1.jl` in scratch, includes compile): Float64 0.518 s before, 0.123 s after; Double64 0.968 s before, 0.118 s after; `check_frame` first call 0.012 s before, 0.006 s after (after the `gali_values` call has compiled). Later calls: Float64 4 us, Double64 47 to 49 us (unchanged within noise from the 30 us to 0.9 ms of the review; BigFloat not re-timed).

Result: gali suite 2827 of 2827 (see the test plan); whole suite `Pkg.test()` with 4 threads: 1 408 380 pass, 0 fail, 0 error, 1 existing broken, 32 min 22 s in the testset (Aqua, JET, layering pass); log `scratchpad/ga05/full-d109.out`. Coverage of `gali.jl` not measured in this rework.
