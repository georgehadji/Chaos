# ext (ChaoslabCUDAExt)

Package extension that gives `propagate` a CUDA device. Loads when `CUDA` is present (weak dependency in `Project.toml`, compat in the same file). Decision D29; launch budget and lost-device handling D33; policy `_config/architecture.md`, "Compute".

Spec: `workspaces/module-build/stages/01-spec/output/tangent-core-spec.md`, Amendments D29 and D33. Report: `stages/04-verify/output/tangent-core-verification.md`, section GPU, and `gpu-verify-log.txt`.

## Files

| File | Holds |
|------|-------|
| `ChaoslabCUDAExt.jl` | the whole extension (one file, 83 lines) |

## What it adds

The extension exports nothing. It fills the hooks that `Foundation` and `Tangent` define:

| Hook | Filled with |
|------|-------------|
| `Foundation.DEVICE[]` | `CUDADevice()` when `CUDA.functional()`, set in `__init__` |
| `Foundation.probe(::CUDADevice, T)` | `loadable`, `functional`, `device` name, `compiles` (one-orbit launch of the real kernel, cached per `T`), `gpu_faster` (from the record) |
| `Foundation.free_memory`, `is_oom`, `is_compile_error`, `is_device_lost` | free device bytes; classifiers for the fallback reasons 4, 3 and 8 |
| `Foundation.gpu_step_ns(::CUDADevice, T)` | measured ns per orbit-step from the record, or `nothing` |
| `Tangent.device_backend`, `Tangent.workgroup` | `CUDABackend()`, work-group size 256 |

The per-machine record is `data/gpu-benchmark.toml` (git-ignored), written by `verify/gpu.jl`, keyed by device name and precision label. It is read at every probe, inside `try`: any read, parse or shape error gives `gpu_faster = false` and no ns value.

## Invariants a caller can rely on

| Invariant | Guaranteed by |
|-----------|---------------|
| Without the extension, or without a functional device, every run resolves to the CPU with a recorded reason; the library never loads CUDA itself | T-AC10, D29 |
| `:auto` picks the GPU for an element type only after a recorded measured win; an explicit `:gpu` is honoured when the device compiles the kernel | D25, T Compute |
| The same kernel source runs on both backends; no second stepping loop | T-AC14 |
| An out-of-memory, a compile failure, a health error seen only on the device, or a lost device reruns the chunk on the CPU with the matching reason | D29, D33b, D33c |
| One launch is kept under the 1 s estimated budget, or the run goes to the CPU with reason 8 | D33a |
| After a lost device `DEVICE[]` is `NoDevice()` for the rest of the session | D33b |
| GPU results agree with the CPU on `S/n` within tolerance, not bitwise | T-AC12 |

## Tests and verification

- Package (no GPU in CI): `test/tangent-core/tc_device.jl` drives the device path with a fake device on the CPU backend, including the fake lost-device case.
- Local GPU run: `verify/gpu.jl` (every GPU criterion of Batch A), on the GeForce RTX 5060 Laptop GPU from D46 (`tangent-core-verification.md`, section "Amendment D46"); the GTX 1050 runs of Batch A, Batch B and D44 stay in the reports. A real launch timeout cannot be provoked on purpose.

## Deferred

| Item | Where |
|------|-------|
| Step-chunked relaunch with state carried between launches (removes the `N` limit from the watchdog budget) | named run before M2 exit, D33 |
| Double64 on the GPU (Quadmath call and heap allocation block compilation) | reason 3 recorded today, `tangent-core-verification.md` |
| Separate "device lost" reason string | proposal, `tangent-core-impl-notes.md` |
| `DEFAULT_GPU_NS` without a record is a judgement value | `tangent-core-impl-notes.md`, Amendment D33 |
| Loading from `bin/chaoslab.jl` | D7, with the entry point |
