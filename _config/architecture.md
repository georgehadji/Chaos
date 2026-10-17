# Architecture, Languages, Patterns

Layer 3, canonical. Numerics rules: `numerics-policy.md`. Code and test rules: `coding-conventions.md`.

## Decision

Core library in **Julia** (`Chaoslab.jl`). **Python** only for N-body adapters (REBOUND ecosystem, P6) and ICM mechanical scripts. The two meet through HDF5 files with a written schema; no in-process bridge.

Why Julia: the plan treats precision as a run parameter (C1, R2 precision arm, idea 5). On native Windows, NumPy `longdouble` is float64 and Numba kernels cannot run double-double or arbitrary precision, so a Python core would turn the precision studies into float64 plus slow mpmath spot checks. In Julia the same kernel runs on Float64, Double64 (DoubleFloats.jl), MultiFloats and BigFloat. Julia also gives forward-mode AD (ForwardDiff) for tangent maps on Windows (JAX is "experimental" there), and ready oracles: ChaosTools (`lyapunov`, `lyapunovspectrum`, `gali`, 0-1 test), TaylorIntegration (`lyap_taylorinteg`, jet transport, BigFloat), OrdinaryDiffEqSymplecticRK (KahanLi8, SofSpa10, McAte, CalvoSanz4 ...).

Rejected: Rust (nothing off the shelf, `f128` still nightly-only); Python+Numba core (fails the precision constraint); juliacall bridge (threading and signal caveats, startup cost per process).

Fallback if the human overrides to Python: P1/P2 survive (float64 + mpmath oracles), R2 precision arm and idea 5 shrink, AD tangents need JAX (experimental on Windows).

## Module map

| Module | Lang | Pattern | Notes |
|--------|------|---------|-------|
| `Foundation` | Julia | value objects (immutable structs), type parameter `T<:Real` for precision | `State{T}`, errors, provenance, compute backend detection and fallback (section "Compute") |
| `Systems` | Julia | multiple dispatch on `AbstractMap`, `AbstractFlow`, `AbstractSeparableFlow`; StaticArrays for dim <= 8 | each system: `step`, `tangent` (analytic Jacobian), `invariants`; maps keep lifted coordinates, wrap only for output |
| `Integrators` | Julia | table-driven strategy: coefficients as data + ONE generic `splitstep`; composition builder (triple jump, Suzuki) | coefficients copied from the cited paper, cross-checked against OrdinaryDiffEqSymplecticRK where it ships the method; hand-rolled: Forest-Ruth, PEFRL, Blanes-Moan, SABA/SBAB, ABA864, Tao, Gauss-Legendre |
| `Tangent` | Julia | pipe-and-filter: `propagate` is the single hot loop, emits snapshots to reducers; one KernelAbstractions kernel per system serves CPU and GPU | carries deviation vectors through the SAME splitting (needed anyway: ChaosTools tangent systems cannot use partitioned symplectic solvers) |
| `Lyapunov` | Julia | strategy via dispatch: `orthonormalize(alg, W)` (pure, returns Q, R and a status; D34) for CGS, MGS, RGS, Householder (LAPACK), H2, Givens, block/sparse; estimators are reducers | LAPACK through LinearAlgebra; hand-written small-dim kernels; Ginelli CLVs with stored R factors and checkpoints |
| `Indicators` | Julia | reducers `(init, update!, finalize)` | SALI, GALI_k, MEGNO, FLI/OFLI, RLI, stretching numbers, NAFF/FMA, WBA, REM, 0-1; all computed in ONE pass over the same orbit for fair cost comparison |
| `Diagnostics` | Julia | reducers | pairing residual, det residual, orthogonality loss, energy / modified-energy error |
| `TimeSeries` | Julia | plain functions | embedding, Wolf, Rosenstein, Kantz, Eckmann-Ruelle (P4) |
| `Shadow` | Julia | symbolic -> generated code: Symbolics.jl BCH Poisson brackets, `build_function` | modified Hamiltonian H_h for R2 |
| `Experiments` | Julia | command (`RunSpec` -> `run`), repository (`Store`), content-addressed cache | TOML.jl configs, HDF5.jl store, manifests, chunked resumable sweeps |
| `Figures` | Julia | scripts | CairoMakie (single language end to end) |
| `nbody` (P6) | Python | adapters | REBOUND/REBOUNDx/celmech in venv A; ASSIST pins `rebound<5`: venv B or WSL2; CPU only, parallel as one process per independent simulation |
| `scripts/` (ICM) | Python stdlib | one script, one job | `icm_status.py`, `icm_validate.py`, `bibcheck.py` (Crossref), `promote.py`, `fetch_horizons.py` |


## Layering

`Foundation <- Systems <- Integrators <- Tangent <- {Lyapunov, Indicators, Diagnostics} <- {Shadow, TimeSeries} <- Experiments <- Figures`. A test parses each file's `using`/`import` lines and fails on an upward edge. Python never imports Julia and vice versa.

## Compute

GPU first, CPU fallback (decision D17). Each run resolves `backend` from `numerics-defaults.toml` (`auto`, `gpu`, `cpu`); `auto` means GPU when CUDA.jl loads and reports a functional device, else CPU.

- One kernel, two backends: hot loops are KernelAbstractions kernels over isbits element types; the CPU backend runs the same kernel with `Threads.nthreads()` workers. The CPU run is the verified reference; a module ships GPU as default only after its cross-backend test passes.
- What runs where: ensembles of independent orbits in Float64 (grids, indicator maps, sweeps) are GPU candidates (Float32 is dropped from M2 and returns at M4 at the earliest, D24); BigFloat, single long orbits and LAPACK-heavy spectra run on CPU threads unless the module spec measures a GPU win.
- Fallback to CPU, recorded in provenance with the reason: CUDA.jl not loadable, no functional device, kernel does not compile for the element type, device out of memory, a `NumericalHealthError` raised only on the GPU, a reducer in the chain that is CPU-only (e.g. a LAPACK factorization), GPU not faster for the element type (`auto` picks the GPU for an element type only after the tangent-core benchmark shows a GPU win on this machine, risk K15; decision D25), or a launch that would exceed the display watchdog ("launch exceeds the display watchdog": estimated launch time above 1 s even for the shortest launch of the step-chunked relaunch; decisions D33, D44), or a lost device ("device lost": a launch timeout or another error after which the device context cannot be used; decision D44). Nine reasons in all, in this order.
- CUDA.jl is a weak dependency (package extension) loaded by `bin/chaoslab.jl`; the library never requires a GPU and CI has none.
- GPU record (decision D29): `gpu_faster` per device and element type is read from `data/gpu-benchmark.toml` under the package directory (git-ignored, per machine; tables `[<device name>."<precision label>"]` with `gpu_faster`, `cpu_ns_per_step`, `gpu_ns_per_step`, `threads`, `date`), written by `verify/gpu.jl`; an absent entry means false.
- Memory: chunk size derives from free device or host memory at run start, never a fixed number; ensembles are structure-of-arrays; workspaces are preallocated per chunk; the hot loop allocates nothing.
- CPU: threads over orbits with dynamic scheduling; BLAS threads set to 1 when Julia threads > 1, so the two pools do not oversubscribe cores.
- Reproducibility: random streams are keyed by orbit id, never by thread or device, so results do not depend on thread count; backends are compared at the preregistered estimator tolerance, never bitwise (chaotic orbits amplify FMA and reduction-order differences).

## Repo shape

```
Project.toml, Manifest.toml      Julia package Chaoslab (Manifest committed for reproducibility)
src/Chaoslab.jl + src/<module>/  library
test/                            Test.jl suites mirroring src/
bench/                           BenchmarkTools suites (perf budgets)
python/nbody/                    optional N-body adapters (own pyproject, GPL deps)
bin/chaoslab.jl                  Julia entry point (panel, run, analyze, figures) into Experiments/Figures
scripts/                         ICM mechanical scripts (Python stdlib)
verify/                          oracle cross-checks for module-build stage 04: own Project.toml + Manifest (ChaosTools, CairoMakie, BenchmarkTools), one script per slug (D27)
schema/run-store.md              HDF5 layout shared by Julia and Python
```

## Scalability

| Axis | Target | Mechanism |
|------|--------|-----------|
| Orbits per study | 10^6 to 10^7 initial conditions | chunked batches, `Threads.@threads` over orbits, resumable chunks keyed by content hash |
| Phase-space dimension | 2 to 200 (coupled maps, lattices) | k-vector tangent blocks, LAPACK QR; banded/sparse QR (R3) |
| Time | 10^9 map steps; Gyr N-body | checkpoint/resume; REBOUND SimulationArchive on the Python side |
| Precision | Float64 to BigFloat | type parameter, same kernel |
| Hardware | GPU or CPU on one machine -> cluster | section "Compute"; job arrays read one TOML + chunk index |
| Storage | 2000x2000 float64 grid = 32 MB | HDF5 per run, attributes for metadata; zarr only if cloud/parallel writes are needed |

## Robustness

Test, oracle, health-check and reproducibility rules live in `coding-conventions.md`, section "Tests".

## Agents

Subagents (Sonnet; Haiku for mechanical work) write specs, tests and code inside module-build stages, one module per run, reading only that stage's Inputs. The language choice is invisible to ICM: stages point to `src/`, `test/`, `scripts/` paths.
