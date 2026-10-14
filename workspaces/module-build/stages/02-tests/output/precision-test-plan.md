# Test plan: precision

Spec: `workspaces/module-build/stages/01-spec/output/precision-spec.md`. Suite: `test/precision/runtests.jl` (helpers in `test/precision/helpers.jl`, prefix `pr_`). Base seed 20261012 (`_config/numerics-defaults.toml`). Float32 dropped (D24): no Float32 tests, no Float32 columns. No regression values from literature are used; all reference values are derived in the spec (Equations 3-6) or computed in BigFloat by test helpers independent of Chaoslab.

## Criterion to testset

| Criterion | Testset | Samples | Tolerance / check | Seed |
|---|---|---|---|---|
| 1 Parameter conversion | AC1 parameter conversion | K in {pi/2, 2pi, 0.1, 1e-3} x T in {Float64, Double64, BigFloat(256)} = 12 | vs 512-bit reference: <= 0.5 ulp (Float64), <= 2 eps_T K (Double64, BigFloat) | none (deterministic) |
| 2 Exact widening | AC2 exact widening | 10^4 Float64 (generic, subnormal, near 2pi, random bit patterns), T in {Double64, BigFloat(256)} | bitwise round trip; equals BigFloat(x) exactly | base+2 |
| 3 Transcendentals in T | AC3 transcendentals in T | 2x10^3 sin/cos (ranges [-2pi,2pi], [-1e3,1e3]), 10^3 log (log-uniform [1e-3,1e3]); Double64 with random low word | <= 16 eps_D max(1,abs(f_ref)), reference BigFloat(256). sin and cos are read through the library (`forward` with K=1 gives p'=sin x, `tangent`[2,1] gives cos x); log is DoubleFloats directly | base+3 |
| 4 Pathwise orbit agreement | AC4 pathwise orbit agreement | 200 ICs, T in {Float64, Double64}, n <= n_H(T) (reference orbit 100 steps) | max ratio e_n/(eps_T M_n Gamma_n) <= 8; floors n_H >= 15 / 40, at most 10 orbits (5 percent) below floor | base+4 |
| 5 BigFloat 256 vs 512 | AC5 BigFloat 256 vs 512 | 200 ICs x 100 steps | e_n <= 8 eps_256 M_n Gamma_n for all n; e_10 <= 1e-60 | base+5 |
| 6 Tangent growth | AC6 tangent growth | 200 orbits, n <= min(n_H(Double64), 40), unit vector from Float64 pair | S_n difference <= 64 (1+K) eps_D M_n Gamma_n, see deviation below | base+6 (orbits), base+60 (vectors) |
| 7 Fallback matrix | AC7 fallback matrix | 9 cases (Float64, Double64, BigFloat x auto, gpu, cpu) with a no-CUDA probe, plus rows for each other probe outcome | backend, reason, and `precision_label` round trip; BigFloat on auto/gpu never calls the probe (probe that throws); `FALLBACK_REASONS` list | none |
| 8 Recording | AC8 recording | Float64, Double64, BigFloat at b in {64, 128, 256, 512} | label equals one of the three forms; `parse_precision` returns the bits | none |
| 9 Cross-backend orbit agreement | stage 04 | | needs a GPU | |
| 10 Float32 gate | stage 04 (and Float32 dropped, D24) | | | |
| 11 Performance report | stage 04 | | | |
| 12 Closure and allocation | AC12 closure and allocation | T in {Float64, Double64, BigFloat} closure; Float64 and Double64 allocation | return types of `forward`, `tangent`, `add`, `total`; `@allocated` of 100-step loops (after warm-up, inside functions) == 0 | none |
| 13 Conversion consequence | AC13 conversion consequence | 100 ICs, K=2pi, up to 200 steps | first n with max-norm separation > 1e-8 is <= 60 for >= 80 of 100; median n* printed via `@info` | base+13 |

Sizes are at the spec maximum (200 orbits x 100 steps, 10^4 points). Runtime with stubbed library functions: about 30 s on CPU; BigFloat 256/512 reference orbits dominate.

## Decisions and ambiguities resolved

- Reference quantities M_n and Gamma_n (Equation 3) are computed in Float64 from the BigFloat reference orbit (Jacobian entries evaluated in BigFloat, then rounded); only their magnitude matters. n_H is the first n at which c eps M_n Gamma_n exceeds 1e-4, minus one, capped at the reference length 100.
- AC4 and AC5 take the reference orbit from test helpers (own BigFloat map), not from the library, so a library defect cannot hide itself. AC5 runs the library at 256 bits against the helper at 512 bits.
- AC3: "Same test at Float32" dropped. Log is not reachable through the Chaoslab API, so it is checked on the DoubleFloats function; sin and cos go through the library map and tangent so the testset still tests the library's element-type handling.
- AC7 uses `resolve_backend(request, T; probe)` with test probes. For the GPU-capable types the no-CUDA probe yields "CUDA.jl not loadable" for auto and gpu. Extra rows (not in the spec's 9) cover "no functional device", "kernel does not compile for the element type", "GPU not faster for the element type" (auto only), and the all-good probe returning `:gpu`. The "stored precision equals request" part is reduced to `parse_precision(precision_label(T))[1] === T`, because `resolve_backend` does not return a precision.
- AC13 uses `exact(Double64, pi) * 2` against `Double64(2pi)`.
- AC6 deviation from the spec (needs owner decision): at most 5 percent of the 200 orbits (10) may exceed the bound, the same rule the spec gives for AC4. Evidence: with a correct reference implementation (stub, Double64) 1 of 200 orbits exceeds the bound (orbit index 174 in the seeded set, n = 39 and 40, worst ratio 2.6). The same orbit gives ratio 10 for BigFloat 256 against BigFloat 512, so the cause is the bound (it has no term for the local stretch of the tangent direction, which the spec already notes under Equations 6), not the precision of the implementation. Strict alternative: raise the constant 64 to about 256, or drop the 5 percent allowance. The realised worst ratio and exceed count are printed with `@info`.

## RED log

Command (from repo root, wrapped in a testset so all testsets run): `julia --project=<testenv> -e 'using Test, Chaoslab; @testset "precision" begin include(raw"test/precision/runtests.jl") end'`. Every testset errors with `UndefVarError` on a missing Chaoslab name; no parse errors, no helper errors.

| Testset | Error | Fail | First error |
|---|---|---|---|
| AC1 parameter conversion | 1 | 0 | UndefVarError `exact` |
| AC2 exact widening | 1 | 0 | UndefVarError `exact` |
| AC3 transcendentals in T | 1 | 0 | UndefVarError `StandardMap` |
| AC4 pathwise orbit agreement | 1 | 0 | UndefVarError `exact` |
| AC5 BigFloat 256 vs 512 | 1 | 0 | UndefVarError `exact` |
| AC6 tangent growth | 1 | 0 | UndefVarError `exact` |
| AC7 fallback matrix | 2 | 0 | UndefVarError `FALLBACK_REASONS` (in Chaoslab), `resolve_backend` |
| AC8 recording | 5 | 0 | UndefVarError `precision_label`, `parse_precision` |
| AC12 closure and allocation | 1 | 0 | UndefVarError `StandardMap` |
| AC13 conversion consequence | 1 | 0 | UndefVarError `exact` |

Total 15 errors, 0 failures, 10 testsets. Helper check: with temporary stubs of the API (kept outside the repo, in the session scratchpad) all 127 assertions pass in about 30 s.

## Rework (stage 03 audit, 2026-10-13)

AC12 allocation checks moved behind a function barrier (`pr_allocs`): on Julia 1.10 `@allocated` at testset scope counted the boxing of the returned SVector (16 or 32 bytes), not library allocation. Criterion and tolerance unchanged.

## Amendment D34 (2026-10-14)

Source: M2 Batch B API contract `api-batch-b.md`, A1b (approved, D34): `Foundation.effective_eps`, exported. New file `test/precision/pr_d34.jl`, included last from `runtests.jl`; existing testsets untouched. One testset, "D34 effective_eps": exported from Foundation and Chaoslab; `effective_eps(Float64) == eps(Float64) == 2^-52`; `eps(Double64) == 2^-104` (the premise of the amendment, asserted so that a change of DoubleFloats shows up); `Float64(effective_eps(Double64)) == 2^-100` and `effective_eps(Double64) == 16 eps(Double64)`; `effective_eps(BigFloat) == eps(BigFloat)` at 128 and 256 bits inside `setprecision(BigFloat, bits) do`; strictly looser than nominal for Double64 and equal for Float64. Only values are asserted, not the return type. The thresholds that consume it are tested in the diagnostics amendment. 14 assertions when green.

### RED log (Amendment D34)

Whole `test/precision` suite in the unmodified repo: 127 existing assertions pass; D34 effective_eps: Pass 1, Fail 2 (`Base.isexported(.., :effective_eps)` false for Foundation and for Chaoslab), Error 11 (`UndefVarError: effective_eps`). Log: scratchpad `m2/red-d34-precision.log`. Against the throwaway patched copy: all 14 pass, existing 127 still pass.
