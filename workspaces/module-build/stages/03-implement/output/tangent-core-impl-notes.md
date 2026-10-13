# Implementation notes: tangent-core

## Foundation part

Code: `src/Foundation/Foundation.jl` (`NumericalHealthError`, `CompensatedSum` with `add`/`total`, `FALLBACK_REASONS`, `BackendChoice`, `resolve_backend`, default probe).

### Status

`test/tangent-core/`: AC9 (compensated summation) 4/4 pass; AC10 backend resolution 25/26 pass, the remaining error is section (f), which calls `propagate` (Tangent module, not yet present). AC11 and others need Tangent.

### Deviations from the spec

- `resolve_backend` checks `:cpu` first (returns reason "" and device "cpu" without calling the probe), then BigFloat, then `cpu_only_reducer`, then the probe in the order loadable, functional, compiles, gpu_faster. An unknown request throws `ArgumentError`.
- `gpu_faster = false` only forces CPU for `:auto`; an explicit `:gpu` request is honoured (the API contract names the "GPU not faster" reason for `auto` only; the test only checks `:auto`).
- Default probe: `loadable = Base.find_package("CUDA") !== nothing`, `functional = false`, never loads CUDA. So without CUDA installed the reason is "CUDA.jl not loadable"; with CUDA installed but the Foundation probe unable to test the device, "no functional device". A package extension replaces this later.
- `NumericalHealthError` has a `Base.showerror` method; fields as in the API contract.
- `add(acc, x)` converts `x` to the accumulator type (`T(x)`).

### Performance

`add`: Float64 7.2 ns, Double64 307 ns per call, 0 allocations.

### Known limits

- Reasons 4 and 5 of `FALLBACK_REASONS` ("device out of memory", "NumericalHealthError raised only on the GPU") are only listed; they are per-chunk and set by Tangent's `propagate`.

## Tangent part

Code: `src/Tangent/` (`Tangent.jl` protocol functions and exports, `initial.jl`, `orbit.jl`, `propagate.jl`); bench `bench/tangent-core.jl`.

### Status

`test/tangent-core/` (run in the `testenv` environment, `-t 2`): 423 pass, 0 fail, 0 error (AC1 to AC11, including AC10(f)).

### Design

- `_orbit_body` is the single per-orbit loop (tangent at the old point, then advance x; nested `tau` loop, so no modulo per step). It returns `(x, W, states, status, step)` and never throws; the kernel and `propagate_orbit` both call it. A bad status (1 state, 2 deviation matrix) stops the orbit at once, so `step` is the failing handoff; a last check after the final handoff covers a reducer that returns a non-finite `W` (step `N`).
- The health check runs before the reducers at each handoff, so a reducer never sees a non-finite `W`.
- `propagate` launches one `_orbit_kernel!` per chunk (CPU backend, workgroup 16 orbits). Per chunk: host fills `W0` (seeded `initial_deviations` or the given `W0`), kernel runs, host scans the status array, lowest orbit id first, throws `NumericalHealthError` (`value` is the first non-finite entry as Float64, `threshold` Inf), otherwise copies the chunk into the result. Reducer summaries (`reducer_final`) are computed on the host.
- Reducer chain is a `Tuple` applied by tuple recursion (`_apply`), so `propagate_orbit` is allocation free for isbits `T` and isbits reducer states (AC5: 0 bytes).
- `initial_deviations`: key = splitmix64(splitmix64(seed) + orbit) feeding `StableRNG`; `randn` d x d in Float64; QR with `sign(diag(R))` fix; first k columns; for `T != Float64` two Gram-Schmidt passes in `T` (needed for orthonormality to eps(Double64)).
- Provenance fields as in the contract; `chunk_size` is the chunk size chosen (given or from `Sys.free_memory()`), not clipped to the number of orbits.

### Deviations from the spec

- Health check order: state before deviation matrix when both are bad.
- Final check after the last step is on the reducer-processed `x`, `W`, reported at step `N`.
- `backend = :gpu`/`:auto` resolving to `:gpu` (only possible with a future CUDA extension) raises an error; this stage has no GPU code. The per-chunk reasons "device out of memory" and "NumericalHealthError raised only on the GPU" are not implemented (stage 04).
- `cpu_only_reducer` is always `false` in `propagate`: there is no CPU-only marker for reducers yet.

### Performance (Julia 1.13.1, 2 threads, measured at top level so only indicative)

- `propagate_orbit`, Float64 `StandardMap`, k = 2, no reducers: about 330 ns per step. Zero allocations inside functions (AC5); a top-level `@allocated` shows 64 bytes from untyped globals.
- Double64 `StandardMap` k = 2: about 17 microseconds per step (Double64 `sin` and `cos` dominate).
- `propagate`, 20000 orbits x 300 steps, Float64: about 130 ns per orbit-step with the default chunk, 250 ns with chunk 1000 (the per-chunk host work, one QR per orbit for the seeded `W0`, is serial and becomes visible).

### Known limits

- Host-side `W0` generation (QR per orbit) is serial; for very short runs it dominates. A threaded fill is possible without changing results (the stream is keyed by orbit only).
- `reducer_init` state types must equal the types `reducer_update` returns; the per-chunk state buffer is typed from `reducer_init`.
- The run takes about 2.5 minutes mostly from AC3 (10^6 steps) and AC9 (10^8 compensated adds).

## Review fixes (stage 04)

- `add(::CompensatedSum{T}, x)` no longer narrows: accepts `x::T`, any `Integer`, and floats with `promote_type(T, typeof(x)) === T`; otherwise `ArgumentError`. The check folds at compile time.
- `propagate`: a `W0` entry with element type wider than T throws `ArgumentError` (checked per orbit as chunks are filled, before any result is returned).
- Provenance gains `headroom` (the `chunk_size` fraction, constant `CHUNK_HEADROOM = 0.5`); existing fields unchanged.
