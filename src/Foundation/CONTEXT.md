# Foundation

Bottom layer: health error, compensated sums, backend resolution, precision labels. Imports no other `Chaoslab` layer (`_config/architecture.md`, "Layering"). Built by the tangent-core run (health, sums, backend) and the precision run (labels, `exact`).

Specs: `workspaces/module-build/stages/01-spec/output/tangent-core-spec.md` (T-) and `precision-spec.md` (P-). Reports: `stages/04-verify/output/`.

## Files

| File | Holds |
|------|-------|
| `Foundation.jl` | the whole module: error type, sums, device hooks, `resolve_backend`, precision functions |

## Exports

| Name | What it does |
|------|--------------|
| `NumericalHealthError` | exception with `what`, `value`, `threshold`, `step`, `orbit`, `backend` |
| `CompensatedSum`, `add`, `total` | Neumaier accumulator (`s`, `c`), update, value `s + c`; `add` rejects a narrowing operand with `ArgumentError` |
| `FALLBACK_REASONS` | the eight reason strings, in the order of `architecture.md` "Compute" |
| `BackendChoice`, `resolve_backend` | resolve `:auto`, `:gpu`, `:cpu` for an element type to `backend`, `device`, `reason` |
| `precision_label`, `parse_precision` | `"Float64"`, `"Double64"`, `"BigFloat(bits)"` and back; bits are read when the label is built |
| `exact` | build a `T` from an integer, float, irrational, rational or decimal string without a Float64 detour |
| `Double64` | re-export of the extended type (D20) |

Not exported, for extensions and tests: `AbstractDevice`, `NoDevice`, `DEVICE`, `probe`, `free_memory`, `is_oom`, `is_compile_error`, `is_device_lost`, `gpu_step_ns`. `ext/CONTEXT.md` says which of them the CUDA extension fills.

## Invariants a caller can rely on

| Invariant | Guaranteed by |
|-----------|---------------|
| The error record has orbit, step, backend and which quantity failed | T Algorithm "Numerical health", T-AC8 |
| A compensated sum has no error growing with the term count; the plain sum fails the same test | T Algorithm "Accumulated logarithms", T-AC9 |
| `resolve_backend` never returns `auto`; `:cpu` gives an empty reason; BigFloat on a GPU request gives reason 3; check order is element type, CPU-only reducer, loadable, functional, compiles, GPU faster | T Compute, T-AC10 |
| Each reason is recorded as the exact string of `FALLBACK_REASONS` | T-AC10(e), D25, D33 |
| `parse_precision(precision_label(T))` returns `T`; bad strings throw `ArgumentError` | P-AC8 |
| `exact` and the Float64-to-`T`-to-Float64 trip lose nothing; no step goes through Float64 | P Invariant 1, P-AC1, P-AC2 |
| Results of `add` and `exact` stay in `T`; both allocate nothing after warm-up for isbits `T` | P Invariant 2, P-AC12 |
| A fallback changes the backend, never the requested precision | P Invariant 4, P-AC7 |
| Float32 has no method (`precision_label(Float32)` is a `MethodError`) | D24 |

## Tests and verification

- Package: `test/tangent-core/tc_ac09_11.jl` (sums, backend, provenance), `test/tangent-core/tc_device.jl` (fake device, per-chunk reasons), `test/precision/runtests.jl`.
- Oracle runs: `verify/tangent-core.jl`, `verify/precision.jl`, `verify/gpu.jl`; reports `stages/04-verify/output/tangent-core-verification.md` and `precision-verification.md`.
- Implementation notes: `stages/03-implement/output/tangent-core-impl-notes.md`, `precision-impl-notes.md`.

## Deferred

| Item | Where |
|------|-------|
| Float32 element type | D24, revisit at M4 |
| CPU-only reducer marker (reason 6 fires today only for a reducer state that is not isbits) | lyap-householder-qr, D25, D30 |
| A separate reason string for a lost device (reason 8 is reused) | proposal in `tangent-core-impl-notes.md`, Amendment D33 |
| Test of `exact(Double64, ::BigFloat)` | `precision-impl-notes.md`, "Known limits" |
| Double64 on the GPU: does not compile, run falls back with reason 3 | `tangent-core-verification.md`, section GPU |
