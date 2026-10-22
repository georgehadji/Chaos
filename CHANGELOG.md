# Changelog

All notable changes to `Chaoslab` and this repository. Versions follow Semantic Versioning (`_config/build-order.md`, section "Releases").

## Unreleased

### Added

- `Chaoslab.Indicators` layer with the `SALI` reducer (Smaller Alignment Index and its log form for `k = 2`), `alignment_indices`, `sali_pair`, `sali_from_log`, `sali_floor` and `check_frame`. Verified at module-build stage 04 against closed forms, the Lyapunov L1 estimator on the identical orbit, ChaosTools GALI_2 and the local GPU (catalog row I2 `verified`; decisions D66, D67, D68).
- `FLI` reducer (Fast Lyapunov Indicator, log form, catalog row I5) in `Chaoslab.Indicators`: compensated log sums of both column lengths read before normalisation, running maxima of the basis (primary) and the first column, maximisers, `replace = false` position for a chain before `SALI` (one-replacer rule, D71, D72), `check_frame(FLI(), W0)` and `fli_unlogged`. OFLI (row I6) is deferred to M5. Verified at module-build stage 04 against closed forms, a BigFloat reference, the L1 and SALI sums on the identical orbit (bitwise), ChaosTools GALI_2 (one-sided) and the local GPU in the D68 form (catalog row I5 `verified`).
- `MEGNO` reducer (Mean Exponential Growth factor of Nearby Orbits, catalog row I8) in `Chaoslab.Indicators`: per column of the seeded frame `Y` (the (1,-1) form), its mean `Ybar` over the handoffs and the (2,0) mean `Ybar20`, from five compensated sums per column (O(1) state per orbit, weights formed in the element type), the log sums bitwise equal to those of `SALI`, `FLI` and `Benettin`, `replace = false` position for a chain before `SALI` (one-replacer rule, D71, D72) and `check_frame(MEGNO(), W0)`. Spec decisions D73, D74. Verified at module-build stage 04 against closed forms, a BigFloat reference, Mestre2011 Eq. (A2) against the FLI and L1 sums on the identical orbit, the ChaosTools tangent system and the local GPU in the D68 form (catalog row I8 `verified`).
- `RLI` reducer (Relative Lyapunov Indicator, catalog row I9) in `Chaoslab.Indicators`, run on `PairMap` of `Chaoslab.Systems` (the product of a map with itself, tangent applied block by block so each block is bitwise the single-orbit run) with a block frame carrying the same seeded `xi_0` for both orbits: per handoff the difference `DL` of the two finite-time Lyapunov indicators and its mean `RLI` over the handoffs, from compensated log sums; `pair_start` (separation along `xi_0`, minimal image, realized `Dx_r`, health error when the separation is lost) and `pair_distance`; `replace = false` position as for `FLI` and `MEGNO`; `check_frame(RLI(), W0)`. Spec decisions D75, D76.

### Changed (catalog and benchmarks)

- Catalog row I2 and the benchmarks rows for SALI: for 2D maps SALI tends to zero for regular orbits too, by a power law, so a fixed threshold gives a horizon-dependent label.

## 0.1.2 (2026-10-21)

First release archived on Zenodo (D64). The computing code is unchanged since 0.1.0.

### Changed (metadata)

- Author name in `CITATION.cff` and `Project.toml`.

### Changed (figures)

- `ChaoslabCairoMakieExt`: one method legend per figure (F1, F5), the Delta_V band named (F1), the pooled W panel titled and the bootstrap band named (F2), class colours and K-index line styles with a legend (F3), the two all-orbit bars set apart from the per-test bars (F4), the threshold drawn per distinct value and named as the median over the orbits (F5), a precision legend (F6), n axes ticked at powers of ten. The figure data and the tables file are unchanged. F1 and F2 are also drawn in two halves (`F1a`, `F1b`, `F2a`, `F2b`: K index 1 to 3, and 4 to 6 with the pooled W panel) with a 21-unit font, about 8 pt at two-column journal width; the full F1 and F2 stay.

## 0.1.1 (2026-10-21)

### Added (D63)

- User documentation with Documenter.jl (`docs/`): a home page with installation and a run example that executes at build, and the API page of every layer from the docstrings; `checkdocs = :exports` fails the build on an undocumented export.
- `README.md`: installation, example, documentation build, citing.
- CI job "Documentation build" (Ubuntu, current Julia).

### Changed

- Docstrings: `QR_PIVOT` and `QR_PASSES` documented; two cross-references in docstrings made resolvable.

## 0.1.0 (2026-10-21)

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

### Added (M3 experiments, D49)

- `Experiments` module (panel scope): screening grid with L5 at Float64, periodic orbits up to the preregistered period by 256-bit Newton with a conditioning test (D51), the preregistered class predicates, seeded draws keyed by cell and orbit id, controls, nested subsets, and a panel file whose manifest checks the preregistration hash. CPU only; byte-identical at any thread count on one machine and Julia version (D52). Package dependency `SHA` added.
- `verify/experiments.jl`: ChaosTools periodic-orbit oracle, independent 1024-bit re-evaluation, thread-count runs, full study inputs, portraits.
- Study `p1-le-benchmark`: preregistration frozen with amendments A1 (T predicate, panel block, caps, hash rule), A2 (superseded) and A3 (conditioning test); panel of 883 orbits built from a clean commit.
- `Experiments` runner and run store (D53, D54, D55, D56, D57): run plan of 15569 work items in 1363 chunks with content hashes, chunked resumable executor (atomic HDF5 file per chunk under schema v2, JSON manifest, health checks, stop and pause rules, cost caps), pilot and CPU rerun check, `bin/chaoslab.jl` (panel, run, pilot), device path for the eight Float64 kernel methods with a per-method recording-chain record (the GPU is 2.9 to 4.1 times slower on the RTX 5060, so `auto` resolves Float64 to the CPU; risk K19 closed), run blocks of successive amendments merged in file order, per-arm stop exemption with per-orbit `breakdown_tau` flags, control scope `exact_start`, per-method control bound. Package dependencies `HDF5`, `JSON`, `Dates`.
- Study `p1-le-benchmark`: amendments A4 (run block), A5 (L2 in wrapped coordinates, tau sweep exempt, control scope, measured budget lines) and A6 (control bound 5e-8 for L9), all before any study data file.
- Manifest writer (D59): the temporary write, the replace and every numbered-backup copy are retried on busy, access-denied and not-permitted errors; a backup that stays held is skipped and listed in the manifest (`backup_failures`) instead of ending the campaign.
- `bin/chaoslab.jl run|pilot --prior-record <file>` (D58): the older GPU record is written to the manifest as `prior_record` (preregistration A4 K19).
- `verify/experiments-runner.jl`, `verify/experiments-runner-reader.jl`: reduced campaign with kills and resume, independent HDF5 reader in a separate process, `-t 1/4/23` runs, pilot against the budget lines, real-horizon control check, device records.

### Added (M3 analysis, D60)

- `Analysis` (new module slug `analysis`, between `Experiments` and `Figures`): reads a finished study run and computes the preregistered tables T1 to T7, the numbers behind the figures F1 to F7, the labels of D36 and D37 and the decisions H1 to H5 by the frozen rules; thresholds come from a parameter table (`read_params`), every input is checked and hashed before and after, refusals name their cause, the tables file is byte-identical at any thread count.
- Figures: package extension `ChaoslabCairoMakieExt` (CairoMakie as a weak dependency, loaded under the verify environment) draws F1 to F7 as PDF files; `Project.toml` gained the weak dependency, the extension entry and the compat bound.
- `bin/chaoslab.jl analyze <preregistration> <manifest> <panel> <parameters> [--store] [--out] [--threads]` (exit code 0, 1 refusal, 2 figures skipped).
- Study stage 04 input: `workspaces/study/stages/04-analyze/references/p1-le-benchmark-analysis-parameters.toml`, the transcription of the frozen values of the preregistration (checked against its text by `verify/analysis.jl`).
- `verify/analysis.jl`, `verify/analysis-independent.jl`: closed forms at K = 0 on the real store, run with figures (54 s, 2.45 GiB), byte-identical output at 1 and 23 threads, independent recomputation by a script that shares no code with the module (883 of 883 orbits, all sets, H1 to H5 and 47 label-count cells equal), the command.
- Corrections of the specification found by the real store: oracle O1 does not hold for L9 (it has its own singular-value closed form); items of the tau sweep with a breakdown hold NaN rows; the Y5 control has no `K`; the ladder needs one common half point for the three precisions.
- Layering: `Analysis` added to the layer list of `test/layering.jl`; `_config/architecture.md` and `_config/build-order.md` list the module (D60).

### Fixed (CI, before v0.1.0)

- Julia floor raised from 1.10 to 1.11 (D62): `Manifest.toml` is now the Julia 1.11 manifest and the CI matrix tests 1.11 and the current release. The old 1.10 manifest lacked the dependencies added in M3 (JSON, HDF5), so the LTS jobs failed at build.
- `write_chunk`: the file descriptor is converted to `RawFD` before `fsync`/`_commit` (Julia 1.11 returns an `Int`).
- Experiments tests: the campaign helper uses a zero-cost clock when a test injects none, so the first-chunk cost cap cannot pause arms on a slow CI host; the cap itself stays tested with injected clocks. The CLI run test sets `--cap-override` for the same reason; the outside-repository path test uses an absolute path of the host platform (Ubuntu).

### M3 exit (D61)

- Modules `experiments` (D49) and `analysis` (D60) verified (stages 01 to 05); study p1-le-benchmark accepted (stages 01 to 05, results document); risk K20 closed, K13 and K18 updated.
- `CITATION.cff` added; first tag `v0.1.0`, archived by the GitHub-Zenodo integration.

### Deferred, with reason

- `bin/chaoslab.jl`, `python/nbody/`: no caller before M3 (D16).
- Coverage reporting: not yet set up; module tests exist, the report is not wired into CI.
- Step-chunked GPU relaunch (state carried between launches): a named run before M2 exit (D33).
- lyap-two-trajectory criterion 13 (external two-trajectory oracle): sent back to stage 01, open before M2 exit (D43).
- Fallback reason for a lost device (error 700 is reported as "launch exceeds the display watchdog"): with the step-chunked relaunch (D33, D43).
