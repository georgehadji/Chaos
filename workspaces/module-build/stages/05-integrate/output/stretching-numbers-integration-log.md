# Integration log: stretching-numbers

Milestone M4 (D65), on `main`. Stage 05 of `workspaces/module-build`.

## Catalog rows flipped

| File | Row or section | Change |
|------|----------------|--------|
| `_config/catalogs/indicators.md` | I10 (stretching numbers, spectra, spectral distance) | planned -> verified; main-issue column states that cells, interval and vector pair set the level and there is no threshold, that in 2D D falls as 1/Z in both classes (regular to chaotic median ratio 9 to 14 at equal Z), and the per-orbit workspace (D78); sources extended with the keys read (Maffione2011a, Darriba2012) and the three not read |
| `_config/catalogs/lyapunov.md` | L13 (finite-time LE distributions, local exponents) | planned -> verified with a qualifier: distribution, mean and second moment of one orbit's stretching numbers; no FTLE field (L14, an Experiments grid product, spec P8) |

Evidence: `stages/04-verify/output/stretching-numbers-verification.md`, accepted at the stage 04 checkpoint (2026-10-27) with Amendments D78 and D79: criteria 1 to 16 pass (package suite folder 3476 checks with the review test; ChaosTools oracle counts equal in all 4120 comparisons; GPU Float64 12 PASS, Double64 3 FALLBACK in the D68 form); independent review with no CRITICAL or HIGH finding, its MEDIUM fixed (handoff limit in the reducer).

Open: none for the module. Recorded limits: edge ties follow the rounding of the cell width, provenance cell width in the range type, std resolution about sqrt(eps), `stretching_range` for the standard maps and 2x2 linear maps only. Sources not read (spec P11): Voglis1994, Voglis1999, Contopoulos2002, PikovskyPoliti2016, Haller2001, Shadden2005; the owner may supply Voglis1999 and Voglis1994. Follow-up for Experiments (obligation in `src/Experiments/CONTEXT.md`): allocate the per-orbit workspace in the device and record paths, run `check_run`, retain counts per spec P9. 4D weak-chaos test (Maffione2011a Sect. 6) waits for `coupled-maps` (M5, spec P12). The full package suite was last run in full at stage 03 (1,396,554 checks); after it only the reducer failure branch and one test file changed, and the folder was rerun (3476).

## Commits

- 0ba578b module-build stretching-numbers stage 01: spec approved (D78): histogram workspace P1(c), cumulative D, L13 covered, L14 to Experiments
- 1ddc3e2 module-build stretching-numbers stage 02: tests (42 testsets, RED) and test plan; spec Amendment D79
- 0062e93 module-build stretching-numbers stage 03: StretchingNumbers reducer, Tangent reducer workspace (additive), bench, docs, impl notes
- 4c1cd9f module-build stretching-numbers stage 04: verification (criteria 1 to 16, CPU and GPU), review R1 fixed (handoff limit in the reducer), AC6 negative-control guard fixed
- 4b3e1d5 module-build stretching-numbers stage 05: catalog I10 and L13 verified; spec notes, Indicators limits, Experiments obligation (workspace, count retention), changelog, integration log

## Files changed by the module

Short hashes of the stretching-numbers commits that touched each file. Stage paths are relative to `workspaces/module-build/`.

- `CHANGELOG.md` (0062e93, 4b3e1d5)
- `_config/catalogs/indicators.md` (4b3e1d5)
- `_config/catalogs/lyapunov.md` (4b3e1d5)
- `_config/decision-log.md` (0ba578b, 1ddc3e2)
- `bench/stretching.jl` (0062e93)
- `docs/src/api/indicators.md` (0062e93)
- `docs/src/api/tangent.md` (0062e93)
- `ext/ChaoslabCUDAExt.jl` (0062e93)
- `src/Experiments/CONTEXT.md` (4b3e1d5)
- `src/Indicators/CONTEXT.md` (0062e93, 4b3e1d5)
- `src/Indicators/Indicators.jl` (0062e93)
- `src/Indicators/stretching.jl` (0062e93, 4c1cd9f)
- `src/Tangent/CONTEXT.md` (0062e93)
- `src/Tangent/Tangent.jl` (0062e93)
- `src/Tangent/orbit.jl` (0062e93)
- `src/Tangent/propagate.jl` (0062e93)
- `test/stretching/runtests.jl` (1ddc3e2, 4c1cd9f)
- `test/stretching/sn_ac01_03.jl` (1ddc3e2)
- `test/stretching/sn_ac04_07.jl` (1ddc3e2, 4c1cd9f)
- `test/stretching/sn_ac08_10.jl` (1ddc3e2)
- `test/stretching/sn_ac11_13.jl` (1ddc3e2)
- `test/stretching/sn_api.jl` (1ddc3e2)
- `test/stretching/sn_helpers.jl` (1ddc3e2)
- `test/stretching/sn_limits.jl` (4c1cd9f)
- `test/stretching/sn_records.jl` (1ddc3e2)
- `test/stretching/sn_selftest.jl` (1ddc3e2)
- `verify/gpu.jl` (4c1cd9f)
- `verify/stretching.jl` (4c1cd9f)
- `stages/01-spec/output/stretching-numbers-spec.md` (0ba578b, 1ddc3e2, 4b3e1d5)
- `stages/02-tests/output/stretching-numbers-test-plan.md` (1ddc3e2, 4c1cd9f)
- `stages/03-implement/output/stretching-numbers-impl-notes.md` (0062e93, 4c1cd9f)
- `stages/04-verify/output/gpu-stretching-verify-log.txt` (4c1cd9f)
- `stages/04-verify/output/stretching-numbers-chaostools.png` (4c1cd9f)
- `stages/04-verify/output/stretching-numbers-spectra.png` (4c1cd9f)
- `stages/04-verify/output/stretching-numbers-verification.md` (4c1cd9f)
- `stages/04-verify/output/stretching-numbers-verify-log.txt` (4c1cd9f)
- `stages/05-integrate/output/stretching-numbers-integration-log.md` (the commit that adds this log)
