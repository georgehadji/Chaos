# Tangent

Propagation of orbit ensembles with deviation vectors: one KernelAbstractions kernel for CPU and GPU, a chain of reducers called every `tau` steps. Imports Foundation and Systems. Built by the tangent-core run (catalog L7 as a parameter, section D2).

Spec: `workspaces/module-build/stages/01-spec/output/tangent-core-spec.md` (T-; amendments D28, D29, D33, D30). Reports: `stages/04-verify/output/tangent-core-verification.md`.

## Files

| File | Holds |
|------|-------|
| `Tangent.jl` | module, exports, the reducer protocol (`reducer_init`, `reducer_update`, `reducer_final`, `reducer_step`, `reducer_failed`, `reducer_error`) |
| `initial.jl` | `initial_deviations`: seeded Haar frames keyed by (seed, orbit) |
| `orbit.jl` | `_orbit_body` (the one per-step loop), `propagate_orbit`, health status codes |
| `propagate.jl` | `propagate`, chunking, device fallback, launch budget, provenance, `chunk_size`, `bytes_per_orbit` |

## Exports

| Name | What it does |
|------|--------------|
| `propagate(sys, x0, N; k, tau, seed, reducers, backend, chunk, W0)` | run all orbits in chunks; returns `PropagationResult` |
| `PropagationResult` | final `x`, `W`, per-orbit reducer `finals`, run `provenance` |
| `propagate_orbit` | the same body for one orbit, serial; throws `NumericalHealthError` |
| `initial_deviations(seed, orbit, Val(d), Val(k), T)` | orthonormal `d x k` frame for one orbit |
| `reducer_init`, `reducer_update`, `reducer_final` | reducer protocol: isbits state, update at each handoff (may replace `W`), host summary |
| `reducer_step`, `reducer_failed`, `reducer_error` | non-throwing failure path for kernels: record in state, host raises with the orbit id (D30 patch for diagnostics) |
| `chunk_size`, `bytes_per_orbit` | memory formula and chunk from free bytes |

## Invariants a caller can rely on

| Invariant | Guaranteed by |
|-----------|---------------|
| Tangent at step `n` is taken at the old point; orbit and `W` advance together | T-I1, T-AC1 to T-AC4 |
| Handoffs occur exactly at multiples of `tau`; `N` must be a positive multiple of `tau` (`ArgumentError` before any work) | T-I2, T-AC8 |
| Replacement transparency: log factors plus final norm equal the unreplaced evolution | T-I3, T-AC2, T-AC3 |
| Orbit result is independent of thread count, chunk size and ensemble position (bitwise, CPU backends) | T-I5, T-AC6 |
| CPU versus GPU agree on `S/n` within the AC12 tolerances, never bitwise | T-AC12, risk K16 |
| Non-finite `x` or `W` is reported at the first handoff with orbit, step, backend and kind; a failed chunk stores nothing; the lowest orbit id wins | T Algorithm "Numerical health", T-AC8 |
| Column `j` of the seeded frame does not depend on `k`; frames are orthonormal to `8 k eps` | T Algorithm "Initial deviation vectors", T-AC7 |
| No wrapping of `x` inside the stream; a state longer than `dimension(sys)` (winding counts) is health-checked in full | T-I7, D28 amendment, W9 |
| The hot loop allocates nothing after warm-up | T-I8, T-AC5 |
| A device launch stays under the 1 s budget, or the run goes to the CPU with reason 8; a lost device or out-of-memory falls back with its reason; provenance names the backend that finished | D33 amendment, D33a to D33d |
| Exactly one per-step update implementation exists | T-AC14 |

## Tests and verification

- Package: `test/tangent-core/` (`tc_ac01_04.jl`, `tc_ac05_08.jl`, `tc_ac09_11.jl`, `tc_state_length.jl`, `tc_device.jl`); reducer failure path also in `test/diagnostics/`.
- Oracle runs: `verify/tangent-core.jl` (BigFloat products, full BigFloat orbit, timings), `verify/gpu.jl` (every GPU criterion, writes the per-machine record `data/gpu-benchmark.toml`, git-ignored).
- Notes: `stages/03-implement/output/tangent-core-impl-notes.md`; bench `bench/tangent-core.jl`.

## Deferred

| Item | Where |
|------|-------|
| Step-chunked GPU relaunch (state carried between launches), so long `N` runs on the device | named run before M2 exit, D33 |
| ChaosTools `lyapunovspectrum` comparisons (regular orbit, chaotic ensemble) | lyap-benettin stage 04, D27, D30 |
| Reducers do not receive the orbit id; the host attaches it | accepted deviation, `tangent-core-verification.md` |
| Adaptive `tau`, banded tangents | catalog L7b, L7c, M5 |
| Float32 | D24 |
| Serial host-side fill of seeded frames (visible for short runs) | `tangent-core-impl-notes.md`, "Known limits" |
