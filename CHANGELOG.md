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

### Added (M2 Batch B)

- `Lyapunov` (lyap-benettin, lyap-gram-schmidt, lyap-householder-qr, lyap-svd-stabilized, lyap-two-trajectory): Benettin `lambda_1(n)` with compensated sums; `QRSpectrum` with CGS, MGS, RGS, Householder, Givens and augmented QR in the kernel and `LapackQR` on the CPU library path; `SVDSpectrum` singular-value exponents for d = k = 2; two-trajectory `propagate_shadow` (lifted and wrapped). Catalog L1 to L6, L9 (d = 2) and Y1 verified (D43).
- Verification scripts per slug, `verify/threads.jl` (1 vs 4 threads, bitwise) and the Batch B sections of `verify/gpu.jl` (CPU vs GPU on every kernel-path estimator, K15 timing record: GPU 6.3 times the 4-thread CPU).
- Spec amendments D35, D41 (SVD c check against the ensemble 1-ulp sensitivity), D42 (two-trajectory GPU check in round-off units), D43 (CGS/MGS bitwise claim withdrawn; ChaosTools initial-frame statements corrected).

### Changed

- The `Manifest.toml` and `bench/` items of "Deferred" below are done (D13, D16).
- Package test suites run in under 30 minutes on 2 threads: slow suites keep one representative case per criterion and run the full cross-product in `verify/` (D38); a PrecompileTools workload caches the hot paths at package precompile (D39).
- The CUDA extension raises the device stack limit to at least 4 KiB at load, so `SVDSpectrum` runs on the GPU instead of failing with an illegal-address error at the 1 KiB default (D40).
- Tangent-core spec Oracle corrected (D30): ChaosTools `lyapunov` is the two-trajectory method and `lyapunovspectrum` the tangent oracle; identical initial vectors are possible. Amendment D30 appended to the tangent-core and diagnostics specs.

### Changed (D46)

- GPU stage 04 moves to the GeForce RTX 5060 Laptop GPU (compute capability 12.0): full `verify/gpu.jl` run and the D44 rerun pass there (D46); D44j accepts one launch per chunk when the budget fits n = 1e5, the D44c rows at tau = 10 run at K = 0.5.
- `propagate` and `propagate_shadow`: the provenance fields `launch_steps` and `launches_per_chunk` are 0 when the first device chunk fails before any launch completed (review finding 1 of D44, D46).
- lyap-two-trajectory criterion 13 runs the Chaoslab side in wrapped mode and passes at 1e-4; catalog L2 is verified without open criteria (D46).

### M2 exit (D48)

- All nine M2 modules verified (stages 01 to 05); catalog I1 verified on the lyap-benettin evidence; risks K19 to K21 added, K15 closed, K6 moved to M4 and M5.

### Deferred, with reason

- `bin/chaoslab.jl`, `python/nbody/`: no caller before M3 (D16).
- Coverage reporting: not yet set up; module tests exist, the report is not wired into CI.
- Step-chunked GPU relaunch (state carried between launches): a named run before M2 exit (D33).
- lyap-two-trajectory criterion 13 (external two-trajectory oracle): sent back to stage 01, open before M2 exit (D43).
- Fallback reason for a lost device (error 700 is reported as "launch exceeds the display watchdog"): with the step-chunked relaunch (D33, D43).
- `CITATION.cff` and Zenodo: after `setup` fills the author data, before the first tag `v0.1.0` (M3).
