# Test plan: tangent-core

Spec: `workspaces/module-build/stages/01-spec/output/tangent-core-spec.md`. API: the Batch A contract (Foundation items NumericalHealthError, CompensatedSum, FALLBACK_REASONS, BackendChoice, resolve_backend; Tangent items initial_deviations, reducer protocol, propagate, propagate_orbit, chunk_size, bytes_per_orbit). Float32 is dropped (D24): no Float32 tests.

Files in `test/tangent-core/`: `runtests.jl` (includes the rest), `tc_helpers.jl` (test maps `tc_Lin`, `tc_Std`, `tc_Fault`, `tc_Boom`; reducers `tc_LogNorm` and `tc_Reader`; Fibonacci, BigFloat, KS and allocation helpers), `tc_ac01_04.jl`, `tc_ac05_08.jl`, `tc_ac09_11.jl`. Base seed 20261012 (`numerics-defaults.toml`), per-criterion offset = criterion number (AC10 uses +10, the contract checks inside AC8 use +8).

## Criterion to testset

| Criterion | Testset | Sample sizes | Tolerances |
|-----------|---------|--------------|------------|
| AC1 Y5 exactness | "AC1 Y5 exactness" | n = 1..35, one orbit, k = 2, tau = 1; also checks the reader sees count n, last n, no gap | 0 (bitwise): W_n = A^n (Fibonacci entries), x_n = A^n x0 |
| AC2 Y5 renormalization | "AC2 Y5 renormalization and exponent" | tau in {1, 5, 10}, n = 1e4, 100 orbits, k = 1; plus chain-order check (10 orbits, tau = 5) | (a) 8 n eps(Float64) against BigFloat(256) ln norm of A^n w0; (b) 2e-3 on S/n; reader sees count N/tau, last N, bad 0 (I2); LogNorm first then reader: max entry <= 1 + 4 eps; reader first: > 1.5 |
| AC3 K = 0 standard map | "AC3 K = 0 standard map" | n in {1, 10, 1000, 1e6}, k = 2; k = 1 from e_2, n = 1e4 | 0 (bitwise) for W_n and x_n; 8 n eps for S_n vs 0.5 ln(n^2 + 1); S_n/n <= ln n / n + 1e-12 |
| AC4 BigFloat product | "AC4 BigFloat tangent product" | T in {Float64, Double64} x {K = 0.5, n = 100; K = 2 pi, n = 20} x 20 seeded points | ||W_n - W_ref||_F <= 8 n eps(T) prod||M_i||_F ||W_0||_F, W_ref at BigFloat(256) along the realized orbit points; final x equals the iterated `forward` bitwise |
| AC5 zero allocations | "AC5 zero allocations" | n = 1e4, k = 2, chains (LogNorm, Reader) at tau 1 and 4; empty chain at n = 100 | exactly 0 bytes of `@allocated propagate_orbit` after warm-up, measured inside a function |
| AC6 thread and chunk invariance | "AC6 thread and chunk invariance" | 4096 seeded orbits, n = 1e3, k = 2, tau = 1, chunk in {64, 1000, default}, plus default W0 | 0: `propagate` equals the serial map of `propagate_orbit`, bitwise (`isequal`), on x, W and reducer finals (accumulated logs) |
| AC7 seeded initialization | "AC7 seeded initialization" | orbit ids 1..1e4 for (d,k) in {(2,1), (2,2), (4,4), (8,8)}, Float64; 200 orbits in Double64; 1e4 for (b) and (c) | (a) 8 k eps(T); (b) bitwise equal across repeats, reversed order, `Threads.@threads`, chunk positions (via `propagate` with the identity map, chunks 64 and 1000), 1e4 distinct matrices, seed+1 differs; (c) KS <= 1.63/sqrt(1e4) |
| AC8 numerical health | "AC8 numerical health" | 6 orbits, one flagged, chunk = 2; overflow map diag(e^3, e^-3) | (a) Inf at step 7, tau = 5: orbit, step 10, "state", :cpu, no result; (b) NaN in W at step 3: step 5, "deviation matrix"; (c) tau = 200 no error, tau = 300 and tau = 100 (N = 300) step 300; `threshold == Inf`, value non-finite; exact steps; (d) contract: N in {7 (tau 5), 0, -5} raises `ArgumentError` before any work (the map throws if called) |
| AC9 compensated summation | "AC9 compensated summation" | 1e8 Float64 terms, 1e7 Double64 terms, 4-term Neumaier case | (a) relative error vs rational <= eps(Float64); (b) plain sum error > 100x the compensated error + 1e-10 (negative control, last assertion; spec measured 1.9e-9); (c) <= eps(Double64) against the exact rational (hi + lo); `[1, 1e100, 1, -1e100]` totals exactly 2 |
| AC10 backend and provenance | "AC10 backend resolution and provenance" | custom probes; 200-orbit runs for provenance | exact strings; FALLBACK_REASONS is the seven strings in order (eight after Amendment D33); (a) default probe on CI; (b) `:cpu` reason ""; (c) `:gpu` without device; (d) BigFloat never calls the probe; (e) reasons 1, 2, 3, 6, 7 forced by probe or `cpu_only_reducer`, check order element type > cpu_only_reducer > loadable > functional > compiles > gpu_faster; (f) provenance fields, `backend` never auto, `chunk_size` equals the chunk given |
| AC11 chunk size | "AC11 chunk size" | free = 2^30, f = 0.5 | 0: 5592405, 2236962, halved 2796202 and 1118481, plus reducer_bytes, f = 0.25, Double64 (176 B per orbit), `bytes_per_orbit` 96, 240 |
| AC12 CPU vs GPU | stage 04 | | local, report in 04-verify |
| AC13 performance record | stage 04 | | local, report in 04-verify |
| AC14 single kernel | stage 04 | | inspection |
| AC10(e) reasons 4 "device out of memory" and 5 "NumericalHealthError raised only on the GPU" | stage 04 | | per-chunk reasons that need a device; cannot be forced on CPU |
| Oracle block (ChaosTools, full-BigFloat orbit, 1 vs N threads in separate processes) | stage 04 | | separate verification environment |

## Spec ambiguities resolved

- AC2 gives no initial point; the Y5 orbit x_n = A^n x0 overflows Float64 near n = 740, so the tests use x0 = 0 (the tangent of a linear map does not depend on x).
- AC4 "the same realized orbit points": the tests iterate `forward` in T and require the kernel's final x to equal it bitwise; the reference product uses those points.
- AC6 "1 thread vs N threads": per the brief, `propagate` is compared with a serial `map` over `propagate_orbit` at the process's thread count; no process spawning.
- AC8 "step 7" and "step 3" are off-by-one tolerant: the fault appears within the interval before the first handoff that must report it (10 and 5).
- `gpu_faster = false` is tested only with `:auto` (the spec says `auto` picks the GPU only if it is faster; whether an explicit `:gpu` ignores it is left open).
- AC7 (a) is also checked in Double64 (200 orbits): W0 must be orthonormal to eps(T), so the QR has to run in T, not in Float64 then widened.
- Provenance field types are not fixed by the API; the tests only require `string(backend)` in ("cpu", "gpu"), `precision == "Float64"`, `device == "cpu"` on CPU, `threads >= 1`, `cpu` a non-empty string.
- The orbit health error with several faulty orbits is not specified; the tests flag exactly one orbit.

Runtime: about 20 s in the RED run (mostly compilation), about 70 s against a throwaway serial stub with generic QR (AC7 32 s of it). Sample sizes follow the spec and were not shrunk; a static-array QR in `initial_deviations` should bring AC7 well below that. AC9 runs 1e8 plain plus 1e8 compensated additions (about 2 s).

## RED log

Run (the extra outer `@testset` makes a failing testset not abort the file; `test/runtests.jl` already nests the suite):

`julia --project=<testenv> -e 'using Test, Chaoslab; @testset "tangent-core" begin include("test/tangent-core/runtests.jl") end'`

| Testset | Result | First error |
|---------|--------|-------------|
| AC1 Y5 exactness | 1 Error | UndefVarError `propagate` |
| AC2 Y5 renormalization and exponent | 1 Error | UndefVarError `initial_deviations` |
| AC3 K = 0 standard map | 1 Error | UndefVarError `propagate` |
| AC4 BigFloat tangent product | 1 Error | UndefVarError `initial_deviations` |
| AC5 zero allocations | 3 Error | UndefVarError `propagate_orbit` |
| AC6 thread and chunk invariance | 1 Error | UndefVarError `initial_deviations` |
| AC7 seeded initialization | 1 Error | UndefVarError `initial_deviations` |
| AC8 numerical health | 1 Error | UndefVarError `propagate` (aborts the testset at the first control run; the contract `@test_throws` checks are not reached) |
| AC9 compensated summation | 2 Error | UndefVarError `CompensatedSum` |
| AC10 backend resolution and provenance | 2 Error | UndefVarError `FALLBACK_REASONS` (then `resolve_backend`) |
| AC11 chunk size | 10 Error | UndefVarError `bytes_per_orbit` (then `chunk_size`) |

Total: 24 Error, 0 Pass, 0 Fail. No parse errors and no errors inside helpers: the type and method definitions in `tc_helpers.jl` are guarded with `@static if isdefined(Chaoslab, ...)`, so the file loads before the module exists and the testsets fail on the missing names.

Helper validation (not committed): a throwaway serial stub of the API, kept in the session scratchpad, makes all 423 assertions pass. A mutant that evaluates the tangent at the new point instead of the old one fails 60 of 160 assertions in AC4 (AC1 to AC3 cannot see it, since their maps have constant tangents).

## Amendment D28 (2026-10-13)

Test file `test/tangent-core/tc_state_length.jl` (included last from `runtests.jl`). One testset, "W9 state longer than the tangent dimension". It defines a local test map `tc_Long` (state length 4, dimension 2: standard map on (x, p), a step counter, a running sum of x; optional Inf in component 3 or 4 at a chosen step) and a reducer `tc_StateLen`, so the suite does not depend on the Systems wrapped map. Seed `TC_SEED + 9`.

| Criterion | Assertion group | Sample size | Tolerance |
|---|---|---|---|
| W9 (a) | `propagate_orbit` vs manual loop, state of 4 components, W is 2 x k | k in {1, 2}, tau in {1, 3}, 30 steps | bitwise |
| W9 (b) | `propagate` vs serial map of `propagate_orbit` (x, W, reducer finals) | 257 orbits, chunk 16 and default, k in {1, 2}, tau = 3, reducers (log norm, reader) | bitwise |
| W9 (c) | Inf in component 3 or 4 | poison at step 7, tau = 5, 5 orbits, chunk 2 | error "state", step 10, orbit 1 |
| W9 (d) | zero allocation of `propagate_orbit` | 1000 steps with two reducers, 100 steps without | 0 bytes |
| W9 (e) | reducer sees the state | `tc_StateLen` | length 4 |
| W9 (f) | `bytes_per_orbit(Val(2), Val(2), T, 4)` and `chunk_size(...; s = 4)` | Float64, Double64 | 112, 208, floor(2^29 / 112); default call unchanged (96, 5592405) |
| W9 (g) | k = 3 with dimension 2 | one call | ArgumentError before any work |
| W9 (h) | existing suites unchanged | stage 03 run of all of `test/tangent-core` | 423 pass before, 423 + 42 after |

### RED log (Amendment D28)

Run: `tc_helpers.jl` plus `tc_state_length.jl` only, before any source change. W9: 2 Pass, 5 Error, 0 Fail. The 2 passes are the default-call assertions that hold today (`bytes_per_orbit` 96 and `chunk_size` 5592405). The errors are not missing names but missing capability of existing names: MethodError `bytes_per_orbit(::Val{2}, ::Val{2}, ::Type{Float64}, ::Int64)` (three assertions), MethodError `chunk_size(...; s)` unsupported keyword (one assertion), and MethodError `_orbit_body(::tc_Long{Float64}, ::SVector{4, Float64}, ::SMatrix{2, 1, Float64, 2}, ...)` ending the testset at the per-orbit body loop (assertion groups (b) to (g) are not reached).

## Amendment D29 (2026-10-13)

Test file `test/tangent-core/tc_device.jl` (included last from `runtests.jl`), one testset "D29 device hook and fallbacks". No GPU in CI: `tc_FakeDevice <: Chaoslab.Foundation.AbstractDevice` has `device_backend = CPU()`, a probe that reports functional and compiling, and a free memory of 2^28 bytes; the test sets and restores `Foundation.DEVICE[]` in `try/finally`. Faults are injected on the host side of a device chunk: a `W0` vector that throws on its n-th access (out-of-memory and compile exceptions classified by test methods of `is_oom` and `is_compile_error` on test-only exception types), a standard map whose first `forward` call anywhere returns NaN (a health error only the first run of a chunk raises), and a reducer with a `BigFloat` state.

| Criterion | Assertion group | Sample size | Tolerance |
|---|---|---|---|
| D29 hooks | defaults: `NoDevice`, `device_backend`, `workgroup`, `free_memory`, `is_oom`, `is_compile_error`, probe record, hooks not exported | one call each | exact |
| D29a | `:gpu` and `:auto` on the fake device vs `:cpu`; `:auto` with `gpu_faster = false`; chunk size from `free_memory` | 40 orbits, 100 steps, tau 4, chunk 16 | bitwise; exact reason strings |
| D29b reason 6 | reducer chain with a `BigFloat` state | 40 orbits | cpu, "reducer is CPU-only", results equal the CPU run |
| D29b reasons 4 and 3 | fault in the first, middle and last chunk | 40 orbits, chunks of 16, 16, 8 | cpu, device `cpu`, exact reason, results bitwise equal |
| D29b reason 5 | NaN only in the first run; an Inf the CPU raises too; an unclassified exception | 40 and 6 orbits | reason string; `NumericalHealthError` with backend `cpu`, orbit 3, step 10; the exception propagates |

Reason 5 cannot be forced on real hardware (a chunk whose health check fails on the device but not on the CPU needs a device-only arithmetic difference); the fake device forces it with a fault that fires once. Run on the CPU backend only; the GPU measurements are `verify/gpu.jl` (D29c).

## Amendment D33 (2026-10-13)

`tc_device.jl` grows by four testsets and the list assertions change. `tc_with_device` and the aliases `tc_F`, `tc_Tg` moved to `tc_helpers.jl`; the AC10 tests that used the default device (`tc_ac09_11.jl`: the default probe and the provenance loop) and the D29 defaults now set `Foundation.DEVICE[]` to `NoDevice()` explicitly and restore it, so no test depends on whether CUDA is loaded. `FALLBACK_REASONS` is asserted with eight strings in `tc_ac09_11.jl` and `test/precision/runtests.jl`.

| Criterion | Assertion group | Sample size | Tolerance |
|---|---|---|---|
| D33a | fake device with a recorded ns per step: chunk capped at 1100; explicit smaller chunk kept; N beyond one wave goes to the CPU with the eighth reason; the default of 100 ns (N = 1e4 gives 1000 orbits, below one wave) | 40 orbits | exact chunk size and reason string; results equal the CPU run |
| D33b | W0 vector throwing a lost-device exception in the middle chunk | 40 orbits, chunks of 16 | eighth reason, `gpu_chunks` 1, `DEVICE[]` is `NoDevice()` after the run |
| D33c | device whose `free_memory` throws out-of-memory; device whose backend `allocate` throws it | 40 orbits | reason 4, `gpu_chunks` 0, results equal the CPU run |
| D33d | provenance after a switch with a small fake free memory; clean device run; CPU run | 40 orbits | `gpu_chunk_size`, `chunk_size` greater than it after the switch, `gpu_chunks` |

The record reader of the extension cannot run in CI; it was exercised locally (absent file, unparsable file, wrong shape, wrong value types) in a scratch script, all giving an empty entry.

## Amendment D34 (2026-10-14)

Source: M2 Batch B API contract `api-batch-b.md`, "Amendments to Batch A code", A1a and A2 (approved, D34). New file `test/tangent-core/tc_d34.jl`, included last from `runtests.jl` (one line appended there); the existing testsets are untouched. Two testsets, `tc_` prefix, test reducers `tc_ProvA` and `tc_ProvB` (read-only, differing only in their provenance; the provenance methods are defined under `@static if isdefined(Chaoslab, :reducer_provenance)`, as `tc_helpers.jl` does for the reducer protocol). Seeds: `TC_SEED` + 3400 to 3410 (above the Batch A offsets).

| Criterion | Testset | Sample | Tolerance |
|---|---|---|---|
| A1a six-argument constructor | "D34a NumericalHealthError column" | `NumericalHealthError("det", 1.5, 0.5, 3, 2, :cpu)` | `column == 0`, other fields as given |
| A1a seven-argument constructor, field order | same | `NumericalHealthError(what, value, threshold, step, orbit, backend, column)` with column 3 | `fieldnames == (:what, :value, :threshold, :step, :orbit, :backend, :column)`, `column::Int` |
| A1a `showerror` | same | both constructors | still names the kind, step and orbit |
| A1a Batch A raisers keep column 0 | same | `propagate` with `tc_Fault(:inf_state)` and `(:nan_W)`, `check_residual`, a `ResidualMonitor` det error through `propagate_orbit` | `column == 0`, kind, orbit, step as before |
| A2 default and export | "D34b reducer_provenance" | `tc_LogNorm`, `tc_Reader`, a reducer with a method | `reducer_provenance(r) === NamedTuple()` by default; exported from Tangent and Chaoslab |
| A2 merge in chain order | same | 4 orbits, `tc_ProvA`, `tc_ProvB` in both orders, `chunk = 16` | later reducer wins a shared key (`shared`); keys added exactly (`alpha`, `beta`, `shared`); all Batch A keys equal to the run without reducers |
| A2 reducers without provenance | same | `(tc_LogNorm, tc_Reader)` and `(tc_LogNorm, tc_ProvA, tc_Reader)` | provenance identical to the empty chain, resp. plus `alpha` and `shared` only |

Counts: 2 testsets, 28 assertions when green. No existing assertion compares a `NumericalHealthError` field count or a provenance key set, so nothing in Batch A had to be edited.

### RED log (Amendment D34)

Run in the test environment (Chaoslab dev'd at the unmodified repo), whole `test/tangent-core` suite inside one testset: 537 existing assertions pass, 3 Fail and 8 Error in the new file, all of them missing capability of the new fields and names, none a helper bug. D34a: Pass 2, Error 2 (`FieldError: type NumericalHealthError has no field column`; `MethodError` for the seven-argument constructor, which ends the testset, so the remaining D34a assertions are not reached). D34b: Pass 4, Fail 3 (`Base.isexported(.., :reducer_provenance)` false, provenance key sets without `alpha` and `shared`), Error 6 (`UndefVarError: reducer_provenance`, `FieldError` on `alpha`). Log: scratchpad `m2/red-d34-tangent-core.log`.

Helper check: the file was also run against a throwaway patched copy of `src/` (scratchpad `m2/lb_patch`, never committed) implementing A1 to A3 as written in the contract: all 28 D34 assertions pass and the 537 existing assertions still pass.

## Amendment D38 (2026-10-15)

Runtime: the package suite carries a reduced set of types, shapes, variants and sample sizes (the reductions are marked "D38" in the test files); everything removed runs at stage 04 in `verify/` with identical assertions (decision-log D38). Where a Double64 part listed above is reduced or moved, the D38 note in the test file names it.
