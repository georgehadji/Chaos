# Changelog

All notable changes to `Chaoslab` and this repository. Versions follow Semantic Versioning (`_config/build-order.md`, section "Releases").

## Unreleased

### Added (M0 bootstrap)

- ICM workspace: four workspaces, 19 stage contracts, `_config/` references, run-store schema v1, setup questionnaire.
- Julia package skeleton `Chaoslab` (empty module), test entry with Aqua, JET and a layering test that also tests itself.
- ICM scripts (Python standard library): `icm_validate.py`, `icm_status.py`, `bibcheck.py`, `promote.py`, with offline unit tests.
- CI on Windows and Ubuntu, Julia LTS and current release, two threads; a separate job runs the ICM checks and script tests.
- Compute policy: GPU first with automatic CPU fallback, one kernel for both, allocation-free hot loops, thread-count invariance, backend recorded per run (decision D17).
- Native Windows build of REBOUNDx 5.1.0 and celmech 1.5.8 (`python/nbody/windows-build/`): patches, build script, smoke tests; closes risk K2 (decision D18).

### Added (M1 literature)

- Four literature runs (lyapunov-algorithms, standard-map-stickiness, chaos-indicators, symplectic-integrators): bibliography 180 to 353 verified entries, 40 located benchmark values, sticky-orbit candidates for the P1 panel (decision D19).
- `bibcheck.py`: DataCite fallback, old-style arXiv IDs, duplicate records under a second DOI form, key suffixes per `conventions.md`, proceedings years, em dashes as `---`. `promote.py`: blocks duplicate titles; flips only keys the topic adds.
- `litmerge.py` (merge per-source search tables, identifier list) and `verify_report.py` (stage 03 outputs, duplicate and chapter drops, checkpoint overrides), named in the literature-review stage 02 and 03 contracts, with tests.

### Added (M2 Batch A)

- `Foundation` (tangent-core, precision): numerical health error, Neumaier compensated sums, backend resolution with eight recorded fallback reasons, precision labels and `exact` conversion; Double64 (DoubleFloats) is the extended type, BigFloat only for short CPU oracle runs (D20).
- `Systems` (standard-map): lifted `StandardMap`, `WrappedStandardMap` with exact integer winding counts (D28), `LinearSymplecticMap` and block builders (catalog Y5 verified; Y1 stays planned until lyap-benettin verifies the large-K known answer), symplectic and determinant residuals.
- `Tangent` (tangent-core): one KernelAbstractions kernel for orbit ensembles with k deviation vectors, reducer chain at every `tau` steps (catalog L7 and D2 verified), seeded Haar initial frames, memory-based chunking, provenance record, non-throwing reducer failure path.
- `Diagnostics` (diagnostics): det, symplecticity, orthogonality, sum and pairing residuals, thresholds as multipliers of `eps(T)`, and the `ResidualMonitor` reducer (catalog section F verified).
- GPU path (`ext/ChaoslabCUDAExt.jl`, D29, D33): CUDA as a weak dependency, probe per element type, per-machine throughput record, launches kept under a 1 s estimated budget, lost-device and out-of-memory fallback to the CPU. Float64 verified on a GeForce GTX 1050; Double64 does not compile for the device and falls back with a recorded reason.
- Verification environment `verify/` (own Project and Manifest, D27): one script per module that prints the measured value behind every acceptance criterion, plus `verify/gpu.jl`; reports, logs and figures under `workspaces/module-build/stages/04-verify/output/`.
- `bench/` timing scripts per module; `Manifest.toml` (Julia LTS) and `Manifest-v1.13.toml` committed (D13).
- Float32 is dropped from M2 (D24); two fallback reasons added (D25, D33).

### Changed

- The `Manifest.toml` and `bench/` items of "Deferred" below are done (D13, D16).
- Tangent-core spec Oracle corrected (D30): ChaosTools `lyapunov` is the two-trajectory method and `lyapunovspectrum` the tangent oracle; identical initial vectors are possible. Amendment D30 appended to the tangent-core and diagnostics specs.

### Deferred, with reason

- `bin/chaoslab.jl`, `python/nbody/`: no caller before M3 (D16).
- Coverage reporting: not yet set up; module tests exist, the report is not wired into CI.
- Step-chunked GPU relaunch (state carried between launches): a named run before M2 exit (D33).
- ChaosTools exponent comparisons: lyap-benettin stage 04 (D27).
- `CITATION.cff` and Zenodo: after `setup` fills the author data, before the first tag `v0.1.0` (M3).
