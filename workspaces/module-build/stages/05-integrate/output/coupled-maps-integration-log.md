# Integration log: coupled-maps

Milestone M5 (first M5 module, D103), on `main`. Stage 05 of `workspaces/module-build`.

## Catalog rows flipped

| File | Row or section | Change |
|------|----------------|--------|
| `_config/catalogs/systems.md` | Y3 (coupled standard maps, Froeschle 4D map) | planned -> verified; system names the form (Skokos2010 Eq. (16)); known-answer column adds the mu = 0 shear conjugacy to Y1 and orbits C and D; sources corrected from Froeschle1972 to Skokos2010, Skokos2003a (Froeschle1972 origin, not read) |
| `_config/catalogs/systems.md` | Y4 (N coupled standard maps, chain) | planned -> verified; system names N range, boundaries, coupling and that general N is derived (D103); known-answer column adds the bitwise eps = 0 limit and the exact origin spectrum; sources corrected from Skokos2004, Skokos2007 to Maffione2011a (N = 2) |
| `_config/bibliography/benchmarks.md` | new rows after "SALI, regular orbit, 2D maps" | "SALI, 4D Froeschle map, orbits C and D" (Skokos2003a) and "Maximal Lyapunov exponent, 4D Froeschle map, chaotic orbits" (Skokos2010; the Skokos2003a text value flagged) |

Y3 and Y4 flip together (P11). The origin spectra and the K = 0 laws are derived closed forms, so they are not benchmark rows (benchmarks hold values read in a source).

Evidence: `stages/04-verify/output/coupled-maps-verification.md`, accepted at the stage 04 checkpoint (2026-11-05): criteria 1 to 21 and D104a to D105d pass (D105d after the rework at stage 04; the Double64 GPU gate is a recorded FALLBACK); folder suite 643 checks; whole package suite 1 405 552 pass, 0 fail, 1 broken (existing, weighted-birkhoff); independent review without CRITICAL or HIGH finding.

Open: none for the module. Recorded limits: the general-N chain and its boundaries are derived, no source prints N > 2 (D103); compile time of the fixed-size chain grows with N k (3.1 s at N = 100, k = 4); Double64 does not compile in the GPU kernel (FALLBACK with reason); criterion 9(c) holds for finite perturbations (spec, clarification at stage 05); `verify/coupled-maps.jl` was not rerun after the stage 04 rework (error paths, promotion and constructor only, covered by the suite).

## Commits

- 2607577 module-build coupled-maps stage 01: spec approved (D103)
- ad0dd00 module-build coupled-maps stage 02: RED tests and test plan; spec Amendment D104
- 2e7228c module-build coupled-maps stage 03: implementation (FroeschleMap, StandardMapChain and wrapped forms), bench, impl notes
- 50fba26 module-build coupled-maps Amendment D105 stage 02: RED tests for the fixed-size chain state
- 6627758 module-build coupled-maps Amendment D105 stage 03: fixed-size chain state, ChainTangent
- 0341413 module-build coupled-maps stage 04: verification, review rework (D105d constant messages, review 2 and 3)
- (this commit) module-build coupled-maps stage 05: catalog Y3 and Y4 verified, benchmarks, changelog, docs, spec clarifications, integration log

## Files changed by the module

Short hashes of the coupled-maps commits that touched each file. Stage paths are relative to `workspaces/module-build/`.

- `CHANGELOG.md` (stage 05)
- `_config/bibliography/benchmarks.md` (stage 05)
- `_config/catalogs/systems.md` (stage 05)
- `_config/decision-log.md` (2607577, ad0dd00, 50fba26)
- `bench/coupled-maps.jl` (2e7228c, 6627758)
- `docs/src/api/systems.md` (stage 05)
- `src/Systems/CONTEXT.md` (2e7228c, 6627758)
- `src/Systems/Systems.jl` (2e7228c, 6627758)
- `src/Systems/froeschle_map.jl`, `src/Systems/residuals.jl` (2e7228c)
- `src/Systems/standard_map_chain.jl` (2e7228c, 0341413)
- `src/Systems/wrapped_standard_map_chain.jl` (2e7228c, 6627758, 0341413)
- `src/Systems/chain_tangent.jl` (6627758, 0341413)
- `test/coupled-maps/` `runtests.jl` (ad0dd00, 50fba26, 0341413), `cm_helpers.jl`, `cm_ac01_04.jl`, `cm_ac05_06.jl`, `cm_ac07_12.jl`, `cm_ac13_16.jl`, `cm_ac17_19.jl` (ad0dd00), `cm_d105.jl` (50fba26), `cm_review.jl` (0341413)
- `verify/coupled-maps.jl`, `verify/gpu.jl` (0341413)
- `01-spec/output/coupled-maps-spec.md` (2607577, ad0dd00, 50fba26, stage 05)
- `02-tests/output/coupled-maps-test-plan.md` (ad0dd00, 50fba26, 0341413)
- `03-implement/output/coupled-maps-impl-notes.md` (2e7228c, 6627758)
- `04-verify/output/` `coupled-maps-verification.md`, `coupled-maps-verify-log.txt`, `gpu-coupled-maps-verify-log.txt`, `coupled-maps-orbits-cd.png` (0341413)
- `05-integrate/output/coupled-maps-integration-log.md` (this commit)

## Audit

| Check | Result |
|---|---|
| Status matches evidence | Y3 and Y4 verified; the report has no open criterion (Double64 GPU FALLBACK is a recorded result of criterion 18, not a deferred part) |
| Log | commits and files listed above; the stage 05 SHA is the commit that adds this file |
| Changelog | `CHANGELOG.md` Unreleased, Added: coupled standard maps entry |

Not changed: `src/precompile.jl` has no coupled-maps workload (the first call at N = 100 compiles in about 3 s); `_config/architecture.md` still names `step` and `invariants` and StaticArrays up to dimension 8 (code: `forward`, `inverse`, `tangent`, `inverse_tangent`; chain states up to dimension 200).
