# Integration log: gali

Milestone M5 (second M5 module, D106), on `main`. Stage 05 of `workspaces/module-build`.

## Catalog rows flipped

| File | Row or section | Change |
|------|----------------|--------|
| `_config/catalogs/indicators.md` | I3 (GALI_k) | planned -> verified (CPU-GPU orbit term D_orb for d >= 4 not formed, D108i); regular-orbit column states the torus law with the maps case (n^-2(k-N), 2D maps n^-2, spec F4, P8); main issue adds the direct-form floor 2 k eps and the log form; sources add SkokosManos2016 |
| `_config/catalogs/indicators.md` | I4 | stays planned; text corrected to "GALI near periodic orbits and fixed points", the laws of I3 tested in gali criteria 11 and 12, no separate indicator; the source `Skokos2008?` replaced by "none: no source on disk defines LGALI" (spec P7, F6) |
| `_config/bibliography/benchmarks.md` | new row after "GALI_k, regular orbit on N-dim torus" | "GALI_k, symplectic maps" (SkokosManos2016 p.19 to 21, Fig. 8, 9, Eq. (33); Skokos2007 p.38) |

The I3 flip closes the M4 clause of D70 ("No gali module in M4", k = 2 only through sali): GALI_k now runs for k = 2..d on Y3 and Y4 (spec P8, settled as proposed at D106).

Evidence: `stages/04-verify/output/gali-verification.md`, accepted at the stage 04 checkpoint (2026-11-05): criteria 1 to 21 pass at full size under Amendments D107 and D108 (21(b) chaotic class reported, D108i; Double64 GPU FALLBACK with reason); ChaosTools `gali` k = 2, 3, 4; folder suite 2724 checks; whole package suite 1 408 277 pass, 0 fail, 1 broken (existing, weighted-birkhoff); independent review without CRITICAL or HIGH finding, its MEDIUM (subnormal column) fixed tests first (D108j).

Amendment D109 (2026-11-05, after stage 05, owner instruction): the orbit-divergence term D_orb for d >= 4 is formed and gates 21(b) chaotic class (PASS on Y3, Y4 N = 2, 5; device ratio at c = 64 at most 0.024), and review findings 2 to 4 are fixed tests first (`ga_d109.jl` RED 18 fail, 5 error; GREEN 103 of 103; whole suite 1 408 380 pass, 0 fail, 1 broken). The I3 qualifier "(CPU-GPU orbit term D_orb for d >= 4 not formed, D108i)" is removed.

Open: none for the module. Not run: the GPU log form for Y4 N = 5 at K = d = 10 (the section runs K = 4 there), horizons beyond n = 8 in the chaotic class. The sali reducer accepts a subnormal column norm (not changed, D108j). Recorded limits: the full-size verify and GPU groups of the D108 rerun ran before the last signature change of `_first_bad` (behaviour identical, whole suite rerun after it); the stage 04 timings were taken on a loaded machine.

## Commits

- cb2e7ce module-build gali stage 01: spec approved (D106)
- e5281f0 module-build gali stage 02: RED tests and test plan; spec Amendment D107
- f5e4469 module-build gali stage 03: implementation, bench, impl notes
- 375cc59 module-build gali stage 04: verification, spec Amendment D108, review rework (D108j)
- 5631b00 module-build gali stage 05: catalog I3 verified, I4 text, benchmark row, changelog, docs, integration log
- (this commit) module-build gali Amendment D109: D_orb for d >= 4 (verify/gpu.jl), review findings 2 to 4 fixed tests first, I3 qualifier removed

## Files changed by the module

Short hashes of the gali commits that touched each file. Stage paths are relative to `workspaces/module-build/`.

- `CHANGELOG.md` (stage 05)
- `_config/bibliography/benchmarks.md` (stage 05)
- `_config/catalogs/indicators.md` (stage 05)
- `_config/decision-log.md` (cb2e7ce, e5281f0, 375cc59)
- `bench/gali.jl` (f5e4469)
- `docs/src/api/indicators.md` (stage 05)
- `src/Indicators/CONTEXT.md` (f5e4469, stage 05)
- `src/Indicators/Indicators.jl` (f5e4469)
- `src/Indicators/gali.jl` (f5e4469, 375cc59)
- `test/gali/` `runtests.jl` (e5281f0, 375cc59), `ga_helpers.jl`, `ga_ac09_12.jl`, `ga_ac13_16.jl` (e5281f0, 375cc59), `ga_ac01_04.jl`, `ga_ac05_08.jl`, `ga_ac17_20.jl` (e5281f0), `ga_d108.jl` (375cc59), `ga_d109.jl` (D109 commit; `runtests.jl`, `src/Indicators/gali.jl`, `verify/gpu.jl`, the spec, test plan, impl notes, verification report and GPU log, `CHANGELOG.md`, `_config/catalogs/indicators.md` and `_config/decision-log.md` also in it)
- `verify/gali.jl`, `verify/gpu.jl` (375cc59)
- `data/literature/gali/` (SkokosManos2016 full text, owner-approved download at stage 01; git-ignored)
- `01-spec/output/gali-spec.md` (cb2e7ce, e5281f0, 375cc59)
- `02-tests/output/gali-test-plan.md` (e5281f0, 375cc59)
- `03-implement/output/gali-impl-notes.md` (f5e4469)
- `04-verify/output/` `gali-verification.md`, `gali-verify-log.txt`, `gpu-gali-verify-log.txt`, `gali-laws.png` (375cc59)
- `05-integrate/output/gali-integration-log.md` (this commit)

## Audit

| Check | Result |
|---|---|
| Status matches evidence | I3 verified; the report has no open criterion under D108; the deferred part (D_orb for d >= 4) is named in the qualifier |
| Log | commits and files listed above; the stage 05 SHA is the commit that adds this file |
| Changelog | `CHANGELOG.md` Unreleased, Added: `GALI` reducer entry |

Not changed: `src/precompile.jl` has no gali workload (first call of the Y4 N = 100, K = 4 chain about 41 s, dominated by the D105 chain frame product, not by GALI); `_config/architecture.md` still says all indicators run in one pass (column reducers other than GALI are k = 2, separate passes on 4D maps, spec P4).
