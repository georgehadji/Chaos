# Implementation notes: clv-ginelli (catalog row L12; spec D122)

Code: `src/Lyapunov/clv_ginelli.jl` (304 lines) and `src/Lyapunov/clv_backward.jl` (139 lines), include and exports in `src/Lyapunov/Lyapunov.jl`, `src/Lyapunov/CONTEXT.md` (78 lines), `docs/src/api/lyapunov.md` (one phrase; the autodocs block picks up the docstrings), `bench/clv-ginelli.jl`. Spec: `../01-spec/output/clv-ginelli-spec.md`. Test plan: `../02-tests/output/clv-ginelli-test-plan.md` (its API proposal and U1 to U15 are the operative reading). Tests are untouched.

## Status

| Run | Result |
|-----|--------|
| New suite `test/clv-ginelli/` (driver `red_clv.jl`, scratch environment `env05`, Julia 1.13.2, 1 thread) | 6251 pass, 0 fail, 0 error after Amendment D123 (three consecutive runs; 32 s quiet, 1 min 46 s under load); before D123 6250 pass, 1 fail (AC10c `a == b`) |
| Whole package suite `Pkg.test()` | 1434634 pass, 2 fail, 0 error, 3 broken (the pre-existing ones), 34 min 57 s; both failures are AC10c `a == b` (Float64 and Double64) of this suite (before D123; not rerun after it) |
| Aqua, JET (`test/quality.jl`), layering (`test/layering.jl`) | inside the same run: quality 12 pass (Aqua and JET clean), layering 2 pass |
| Coverage, `--code-coverage=user`, new suite alone | `clv_backward.jl` 83 of 83 executable lines; `clv_ginelli.jl` 147 of 153 (96 percent); the six missed lines are the `n_b = 0` start of the backward pass (`fill!` and the first `_deliver!`) and the exception unwrapping of a failing task in `propagate_clv`; both were checked by hand in a scratch script (below). The `.cov` files were deleted |
| Not run | `verify/` scripts, `CLV_FULL_RUN` sizes, device runs, `--check-bounds=yes` (stage 04) |

## Design

- Two files. `clv_ginelli.jl`: `ClvGinelli{M, A}` (map, tau, h, orthonormalizer, storage, block), `ClvResult{T, X}`, `ClvEnsemble`, the forward handoff `_clv_handoff`, the two run drivers `_run_store!` and `_run_ckpt!`, `propagate_clv_orbit`, `propagate_clv`, `clv_initial_coefficients`, `clv_bytes`. `clv_backward.jl`: packed storage of R, `_bstep!` (back-substitution and normalisation), `_deliver!`/`_store!` (one sample), `_angles!`, `_bwd_range!`. Longest function: `_store!` 34 lines (`propagate_clv` 27, `_run_ckpt!` 26).
- The forward handoff is the delivered one: `W = fused_product(tangent(map, z), W); z = forward(map, z)` per step (tangent at the old point, one call of each per step, U11), then the state and frame health checks, then `orthonormalize(alg, W)` of the delivered orthonormalizers, then the window conditioning rule `G_w eps < 1/2` for CGS, MGS and RGS (not for the Householder family), exactly as `QRSpectrum`. No new factorization code. `lambda_fwd` is the compensated sum of `ln R_ii` over the window divided by `n_s tau` and by `h`: bitwise `QRSpectrum(alg).lambda ./ h` (U9, tested).
- Storage: R is stored packed upper-triangular (`k(k+1)/2` numbers per handoff) in one `Matrix{T}`; the frames n = 0 .. n_s in one `Array{T,3}` only when `:V` or `:G` is requested and there is no sink. `bytes` is the formula `clv_bytes` of the run (frames term only when frames are held), so the formula is the code's own layout. In checkpointed mode `bytes` is the checkpointed formula, an upper bound when no frames are requested (the block frame buffer is then not allocated).
- Checkpointed mode (V17): first pass through the blocks 0 .. n_c - 2 stores (z, G) at each block start and accumulates the forward-exponent sum; then block n_c - 1 down to 0: rerun the block forward storing R (and the frames), sweep backward. The last block is rerun first and adds its intervals to the window sum, so the sum is in forward order and `lambda_fwd` equals the store-all value bitwise. Forward handoffs counted n_f + (N - h_last) + N, tested.
- Backward step: own scalar loops on preallocated `Matrix{T}` (column j solves the top-left j x j block, so column j does not depend on k), the norm is Foundation's `safe_norm` on a static `SVector` (no allocation), column 1 is set to `e_1` explicitly (U13), `Lambda = 1 / norm`, local exponent `-log(norm) / (tau h)`. A zero or non-finite norm throws "backward breakdown" with `step = n_f + n + 1` (the global handoff of the interval solved) and `column`.
- Angles: principal angles between the spans of two column sets of `C_n` (equal to those of the `V_n` sets because `G` is orthonormal), `qr` and `svdvals` in the element type, `acos` in the element type (U14), ascending, `min(m1, m2)` of them. Not part of the allocation-free path.
- `clv_initial_coefficients`: `StableRNG(_mix(_mix(seed) + orbit) xor label)` with the splitmix finalizer of the tangent stream (`Tangent._mix`, imported, not duplicated; `hash` of a tuple is not stable across Julia versions), column j drawn after columns 1 .. j - 1 (independent of k), normalised in the element type, entry (1, 1) set to 1.
- `propagate_clv`: orbit id = index; chunks of `chunk` orbits (default from a quarter of the free host memory by `clv_bytes`), `ntasks` tasks (default `Threads.nthreads()`) per chunk through `Threads.@spawn`; every orbit writes its own slot, so the result does not depend on chunk, ntasks or position (tested 1 against 4 tasks, chunk 64 and 1000). A failing task is unwrapped to the original `NumericalHealthError` (first orbit that fails, as `propagate_cqr`). `:gpu` and `:auto` resolve through `resolve_backend(...; cpu_only_reducer = true)`: backend `:cpu`, reason "reducer is CPU-only".
- Layering: Lyapunov now also imports `..Tangent._mix` and `..Foundation.fused_product` (both lower layers, not exported) and `StableRNGs` (an existing dependency).

## Deviations from the spec and the test plan

| Item | Deviation |
|------|-----------|
| Exports | beside the six names of the plan, `ClvEnsemble` (the type of the `propagate_clv` result, which the plan describes only as "an object with `results` and `provenance`") |
| Result of `propagate_clv_orbit`, `ClvResult` | exactly the fields of the plan |
| Provenance of `propagate_clv` | adds `device`, `precision`, `seed`, `k`, `bytes_per_orbit` to the four keys of the plan |
| `sink` arguments | `C` and `lambda` are the run's working matrix and vector, valid during the call (copy to keep); the test copies |
| Orthogonality check | none beyond the delivered factorization status and the conditioning rule: the spec names no orthogonality guard for this module |
| GPU | none (P9); a `:gpu` request falls back with the reason of a CPU-only reducer |
| BLAS threads | not set to 1 (the spec Compute paragraph mentions it as policy): the forward and backward loops call no BLAS; `LapackQR` and `qr` of the angle sets are the only users |
| No bound or test was changed | |

## Test changed at stage 02 (Amendment D123, owner: "Amend AC10c")

AC10c (`test/clv-ginelli/clv_ac06_10.jl:202`), `@test a == b`, with `a`, `b` the `@allocated` bytes of a whole run with n_f = 100 and n_f = 300: fails at random (observed 20965 against 20981 Float64, 37020 against 37068 Double64; about 1 to 2 failures in 20 repetitions of the testset, more under load: 2 of the 2 repetitions in the coverage run). Cause: Julia's allocation counter is not deterministic for arrays larger than the pool limit. A bare `zeros(Float64, 2, 2, 101)` (3232 bytes) measures 3351 to 3399 bytes over eight consecutive calls in a fresh Julia (`@allocated`, no package code), so the result arrays of any run (`C`, `cosines`, the packed R) add a random difference of some tens of bytes to `b - a`. The same noise makes the second assertion (`0 <= b2 - a2 <= 48 + 256`) fail with some probability (the expected increase is only 48 bytes). The module itself allocates nothing per handoff: measured with `@allocated` on `_transient` (1000 handoffs, Float64 and Double64, Householder and RGS) 192 and 304 bytes, equal for 10 and 1000 handoffs (the dynamic dispatch of the call), and on the backward sweep `_bwd_range!` over 400 handoffs 0 bytes. Change made under D123: `abs(b - a) <= 256` and `-256 <= b2 - a2 <= 48 + 256` (identical runs measured up to 144 bytes apart in checkpointed Double64 at stage 03).

## Calibration margins (realized, from the `@info` lines of the run)

| Test | Bound | Realized ratio |
|------|-------|----------------|
| AC1 exponents / vectors | 1 | 0.039 / 0.195 |
| AC2 window mean against the delivered spectrum | 1 | 0.989 (as in the plan) |
| AC3a worst relative slope deviation | 0.2 | 0.0064 |
| AC3b worst distance | 1e-10 | 1.4e-12 |
| AC3c error / bound | 1 | 0.247 |
| AC5 (a) / 512 eps | 1 | 0.078 |
| AC6b / 16 eps | 1 | 0.625 |
| AC7 worst sine and cosine | 1e-6 | 6.1e-14 |
| AC11 error / B_v, Float64 F3 orbit C | 1 | 0.141 (others below 2.5e-4) |

These equal the stage 02 scratch values to the digits shown where the plan lists them, because the arithmetic is the same.

## Performance (`bench/clv-ginelli.jl`)

Benchmark environment, 1 thread, medians, orbit C (KahanLi8, h = 1/20, tau = 20), n_f = n_b = n_s = 100, k = 4 (300 handoffs):

| Entry | Time | Allocations |
|-------|------|-------------|
| store-all, no outputs | 11.1 ms | 423 (86.6 kB, the result arrays and packed R) |
| store-all, `:V` | 10.8 ms | 428 (112.7 kB) |
| checkpointed (block 25), `:V` | 16.6 ms | 431 (90.3 kB) |
| store-all, RGS | 11.3 ms | 423 |
| store-all, angles (1:2, 3:4) | 11.3 ms | 3662 (405 kB) |
| store-all, k = 1 | 10.4 ms | 401 |
| delivered discrete RGS spectrum, same 300 handoffs (forward reference) | 10.3 ms | 0 |
| standard map, 1000 handoffs, k = 2 | 0.19 ms | 86 |
| 256 standard-map orbits, 300 handoffs each, 1 task | 16.2 ms | 21791 |

The forward pass dominates (store-all is 1.08 times the discrete spectrum; the backward sweep adds little); the checkpointed run costs 1.5 times the store-all run (about the doubled forward work of V17 on the share that is reran).

## Known limits

- CPU only. The frames and R of all N handoffs are held in `:store` mode (288576 bytes for d = k = 4, N = 2000); use `storage = :checkpoint` (24528 bytes at `block = 100`) for long runs.
- `bytes` is a formula (V17), not a measurement; result arrays handed to the caller are not counted.
- The angle routine allocates (`qr`, `svdvals`) and is meant for angle requests only.
- Ensemble results are all kept (`ClvEnsemble.results`); `chunk` bounds the number of orbits in flight, not the size of the result.
- `n_b = 0`: the last column of `lambda` and `R` is NaN by convention (U2); the starting `C0` is delivered as sample `n_s` (tested by hand, not by the suite).
- Float32 not run (D24); BigFloat runs (AC11b, AC12b) allocate per handoff, as expected.

## Audit

| Check | Result |
|-------|--------|
| GREEN | new suite 6251 of 6251 after D123 (AC10c tolerance 256 bytes); every earlier module passes (whole suite above, run before D123) |
| Quality | inside the same run: quality 12 pass (Aqua and JET clean), layering 2 passITY |
| Size | longest function `_store!` 34 lines; `clv_ginelli.jl` 304 lines, `clv_backward.jl` 139, `Lyapunov.jl` 36; `src/Lyapunov/CONTEXT.md` 78 lines |
| Tests untouched | `git diff` of `test/clv-ginelli/` since 02: only the two AC10c assertions, changed under Amendment D123 (owner: "Amend AC10c"); nothing else |
| Numerics | no `@fastmath`; the forward-exponent sums are `CompensatedSum` accumulators (one per column); the local exponents are single logarithms, not accumulated sums |

## Documentation and code inconsistencies

- The spec writes the forward exponents per `tau h` of the window; the test plan U9 and the code give them as the delivered per-step value divided by `h`. The code follows the plan.
- The spec says angles "smallest principal angle only"; the plan (U2) and the code return all principal angles ascending.
- The spec's `n* = 25 / gamma` guidance is superseded by `40 / (g_min Dt)` (U3); the code has no default transient.
- The approved plan lists the field `lambda_fwd` for "the window n = 0 .. n_s - 1" while samples are n = 0 .. n_s; the code accumulates over the intervals n = 0 .. n_s - 1 as the plan says.
