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

## Amendment D28 (2026-10-13)

Code: `src/Tangent/orbit.jl` (`_orbit_body` binds state length `s` and tangent dimension `d` separately; arithmetic unchanged), `src/Tangent/propagate.jl` (`propagate` reads `d = dimension(sys)`, checks `s >= d` and `1 <= k <= d`, passes `Val(d)` into `_propagate`; state buffers hold `SVector{s,T}`, deviation buffers `SMatrix{d,k,T}`; `bytes_per_orbit(Val(d), Val(k), T, s = d)` and `chunk_size(...; s = d)`).

### Result
`test/tangent-core/runtests.jl` with the full Chaoslab: 465 pass, 0 fail, 0 error (423 before, 42 new in W9). Machine loaded: the run took 4 min 15 s. Bit-identity for s = d: no arithmetic changed, and every existing assertion (including the bitwise thread and chunk invariance of AC6) still passes.

### Deviations
- `propagate` now needs `dimension(sys)`; every Batch A system and every test map defines it. A map without it fails with a MethodError at the k check.
- `propagate_orbit` is untyped in the state length already; it accepts s > d with no check. A mismatch between `size(W0, 1)` and the map's tangent size shows up as a DimensionMismatch in the first step.
- Health: `_status` already scanned every component of the state, so counts in components 3 and 4 are covered without a code change; the test W9 (c) pins this.

### Performance
No measurable change for s = d (same method bodies after specialization). For s = 4, d = 2, `propagate_orbit` costs the map's own time (see the standard-map notes: about 260 ns per step wrapped, Float64, k = 2, loaded machine) and allocates 0 bytes.

### Known limits
- The memory formula adds `sizeof(T) (s - d)` per orbit (32 bytes in Float64 for the wrapped system); reducer state is separate.
- A reducer sees the full s-component state; reducers written for d components that index `x` by position keep working, since the first d components are the phase-space point by convention.

## Amendment D29 (2026-10-13)

Code: `src/Foundation/Foundation.jl` (`AbstractDevice`, `NoDevice`, `DEVICE`, `probe(::NoDevice, T)` replacing `default_probe`, `free_memory`, `is_oom`, `is_compile_error`, none exported; `resolve_backend(...; probe = nothing)` calls `probe(DEVICE[], T)` when no zero-argument probe is passed), `src/Tangent/propagate.jl` (`device_backend`, `workgroup`, `_alloc_bufs`, `_run_chunk!` on any backend, `_device_chunk!` with the single `try`; the `error("GPU propagation needs ...")` line is gone), `src/Tangent/Tangent.jl` (imports the hooks), `ext/ChaoslabCUDAExt.jl` (about 60 lines), `Project.toml` (`[weakdeps] CUDA`, `[extensions]`, compat `CUDA = "6"`, `TOML` in `[deps]`), both root manifests re-resolved. `test/layering.jl` walks `src/` only, so `ext/` needs no change there.

### Result
Full suite with the new `tc_device.jl` (42 assertions): 1617 pass, 0 fail on Julia 1.13.1 (19 min 54 s) and on Julia 1.10.12 (122 min 12 s), Aqua and JET included. `verify/gpu.jl` prints 19 lines (11 PASS, 1 FAIL, 6 FALLBACK, 1 INFO), see `gpu-verify-log.txt`.

### Deviations from the brief
- Compat is `CUDA = "6"`, not `"5"`: the depot has 6.4.2 (CUDA.jl 6 splits into CUDACore and CUDATools; `CUDA.KernelError` and `CUDA.InvalidIRError` are public there). Version 5 was not tried.
- `CUDA.free_memory()` is the free-memory call; `CUDA.available_memory` does not exist in 6.4.2.
- The default chunk is sized from free device memory (about 13 million orbits at 96 bytes), but a launch lasts as long as the number of steps per orbit says, whatever the orbit count: 1024 orbits x 1e5 steps is one launch of 1.5 to 2.9 s. The brief's remedy (fewer orbits per chunk) does not shorten such a launch. This machine has TdrDelay 8 s; with the default 2 s a run of that length would need a smaller n or a longer delay. The kernel was not changed.
- `gpu_faster` for Float64 is true (21.4 vs 91.8 ns/orbit-step on 4 threads), not false as the brief expected, so with CUDA loaded `:auto` picks the GPU for Float64 on this machine.
- Reason 5 is forced in `tc_device.jl` by a fake device and a fault that fires once; on real hardware it cannot be forced.
- `verify/Manifest.toml` was regenerated with CUDA 6.4.2 pinned: Atomix 1.6.0 to 1.5.0, ColorTypes 0.12.3 to 0.12.4 and OrderedCollections 2.0.2 to 2.0.3 changed as a side effect, 38 packages were added. Earlier verify logs were produced with the old manifest.
- AC12a "regular class" is defined by agreement of the CPU Float64 S/n with a Double64 CPU run within 1e-6 (1004 of 1024 orbits); on it the CPU-GPU maximum is 1.162e-06, a FAIL against the proposed 1e-6 by 16 percent.

### Known limits
- Double64 does not compile for the GPU (Quadmath `__subtf3` lazy initialisation, `Vector{Base.Padding}` allocation); the probe caches `compiles = false` per element type and the run falls back with reason 3.
- The `compiles` probe launches one orbit of `StandardMap` with no reducer; a reducer that does not compile is caught per chunk (reason 3 after the fact).
- The record file is read at every probe (a small TOML parse), not cached, so a new measurement takes effect without a restart.

## Amendment D33 (2026-10-13)

Code: `src/Foundation/Foundation.jl` (eighth `FALLBACK_REASONS` string, hooks `is_device_lost` and `gpu_step_ns`, neither exported), `src/Tangent/propagate.jl` (`LAUNCH_BUDGET_S`, `LAUNCH_WAVE`, `DEFAULT_GPU_NS`, `_budget_orbits`, `_fallback_reason`; `_propagate` sets up the device inside the same classification as a chunk, runs a `while` loop so the chunk size can change at a fallback, and records `gpu_chunks` and `gpu_chunk_size`; the docstring states that bitwise equality holds for CPU backends only), `ext/ChaoslabCUDAExt.jl` (`is_device_lost` for the CuError codes 4, 100, 201, 214, 700, 702, 709, 710, 714, 715, 716, 719, 999; `record_entry` reads the record inside `try`; `gpu_step_ns`; the `COMPILES` cache under a `ReentrantLock`), `verify/gpu.jl` (record written through a temporary file and `mv`; AC13 and the record first; AC12a symmetric), `verify/tangent-core.jl` (the reasons-list check has eight strings).

### Deviations and notes
- Reason for a lost device: no new string; "launch exceeds the display watchdog" is used for code 702 and, for lack of a separate string, for the other context losses. Proposed separate string if the owner wants one: "device lost".
- `DEFAULT_GPU_NS = 100` is a judgement (about 5 times the measured 21.4 ns); with it an explicit `:gpu` run without a record refuses N above 9765 with one wave. `verify/gpu.jl` measures AC13 first so a fresh machine records its own value before the long-N lines.
- The device chunk is capped by the budget even when the caller passes a larger `chunk`.
- A real 702 cannot be provoked on purpose; the lost-device path is tested with a fake exception, and `is_device_lost(CuError(ERROR_LAUNCH_TIMEOUT))` was checked in a scratch script.
