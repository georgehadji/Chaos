# Integration log: clv-ginelli

Milestone M5 (catalog row L12, D122, D123, D124), on `main`. Stage 05 of `workspaces/module-build`.

## Catalog rows flipped

| File | Row or section | Change |
|------|----------------|--------|
| `_config/catalogs/lyapunov.md` | L12 | planned -> verified (CPU only: no device kernel); method cell restated (the delivered discrete QR forward, storing R, backward iteration of an upper-triangular C with column normalisation); "Gives" cell names the first k CLVs, local exponents, cosines and principal angles, the forward exponents; main-issue cell restated from the evidence (memory 288 kB per orbit at d = k = 4, n_s = n_b = 1000; checkpoint overhead 1.4 to 1.8x; transient 120 / (g_min Dt), D124; near-tangencies; angle floor about sqrt(eps)) |
| `_config/catalogs/lyapunov.md` | L12 source cell | Ginelli2007, WolfeSamelson2007, KuptsovParlitz2012 -> Ginelli2007, Ginelli2013 (Sect. 4), KuptsovParlitz2012 (arXiv versions, session scratch copies); WolfeSamelson2007 named as not on disk (spec Sources) |
| `_config/catalogs/lyapunov.md` | L12b | unchanged, planned (held, D122; the cross-check of criterion 13 is the KuptsovParlitz2012 static algorithm, not L12b) |
| `_config/bibliography/lyapunov.md` | Ginelli2007, Ginelli2013, KuptsovParlitz2012, WolfeSamelson2007, Froyland2013 | read status recorded: the three arXiv texts read (scratch copies, not in `data/literature/`), the other two not on disk (spec Sources) |
| `_config/bibliography/benchmarks.md` | new row "Lyapunov exponents of the generalized three-dimensional Henon map" | KuptsovParlitz2012 p.18 values with the stage 04 R6 reproduction (about 1e-3, not to the printed digits) |

Evidence: `stages/04-verify/output/clv-ginelli-verification.md`, accepted at the stage 04 checkpoint (owner: "Accept", 2026-11-19): 14 of 14 criteria pass after Amendment D124; oracles BigFloat(256) closed forms, identities of the symplectic pairing and the flow direction, the same-orbit discrete QR spectrum (bitwise), the KuptsovParlitz2012 static null-space algorithm in BigFloat(256) (criterion 13: worst 7.6e-11 Float64, 8.1e-28 Double64, 4.0e-46 at 256 bits against 1e-30); bitwise 1 against 8 CPU threads; independent review without CRITICAL or HIGH finding; M1, M2 fixed under D124.

Open: review findings L1 to L9 (LOW). Not applicable: a device part (no kernel in this release, spec Compute).

## Stage 05 changes

- `_config/catalogs/lyapunov.md` row L12, `_config/bibliography/lyapunov.md`, `_config/bibliography/benchmarks.md` (above).
- `CHANGELOG.md`: Unreleased, Added, the covariant Lyapunov vectors entry.
- `src/Lyapunov/CONTEXT.md`: the list of runs in the header sentence (punctuation) and a Deferred row for clv-ginelli (device kernel, L12b, adjoint CLVs; `verify/clv-ginelli.jl`, `CLV_FULL_RUN`) (79 lines).
- Local checks: `julia --project=docs docs/make.jl` exit 0, no error or warning; after D124 the package clv-ginelli suite 6258 pass, Aqua, JET and layering 14 pass, FULL1 to FULL3 44872, 17563, 151 pass; whole package suite 1434634 pass, 2 fail (AC10c, fixed by D123), 3 broken at stage 03 (not rerun after D123 and D124; CI runs it on push).

## Commits

- 0c4722c module-build clv-ginelli stage 01: spec approved (D122)
- 7d55735 module-build clv-ginelli stage 02: RED tests and test plan
- 9c0e4d4 module-build clv-ginelli stage 03: Ginelli CLVs in Lyapunov (GREEN), Amendment D123
- 4fcfe43 module-build clv-ginelli stage 04: accepted after Amendment D124
- (this commit) module-build clv-ginelli stage 05: catalog row L12 verified, changelog, integration log

## Files changed by the module

Short hashes of the clv-ginelli commits that touched each file. Stage paths are relative to `workspaces/module-build/`.

- `CHANGELOG.md` (stage 05)
- `_config/catalogs/lyapunov.md`, `_config/bibliography/lyapunov.md`, `_config/bibliography/benchmarks.md` (stage 05)
- `_config/decision-log.md` (0c4722c, 9c0e4d4, 4fcfe43)
- `bench/clv-ginelli.jl` (9c0e4d4)
- `docs/src/api/lyapunov.md` (9c0e4d4)
- `src/Lyapunov/` `clv_ginelli.jl` (9c0e4d4, 4fcfe43), `clv_backward.jl`, `Lyapunov.jl` (9c0e4d4), `CONTEXT.md` (9c0e4d4, stage 05)
- `test/clv-ginelli/` `runtests.jl` (7d55735); `clv_helpers.jl`, `clv_ac01_05.jl`, `clv_ac11_14.jl` (7d55735, 4fcfe43); `clv_ac06_10.jl` (7d55735, 9c0e4d4)
- `verify/clv-ginelli.jl` (4fcfe43)
- `01-spec/output/clv-ginelli-spec.md` (0c4722c, 4fcfe43)
- `02-tests/output/clv-ginelli-test-plan.md` (7d55735, 9c0e4d4, 4fcfe43)
- `03-implement/output/clv-ginelli-impl-notes.md` (9c0e4d4)
- `04-verify/output/` `clv-ginelli-verification.md`, `clv-ginelli-verify-log.txt`, `clv-ginelli-d124-verify-log.txt`, `clv-ginelli-convergence.png`, `clv-ginelli-angles.png` (4fcfe43)
- `05-integrate/output/clv-ginelli-integration-log.md` (this commit)

## Audit

| Check | Result |
|---|---|
| Status matches evidence | L12 verified on 14 of 14 criteria; the CPU-only qualifier names the absent device part |
| Log | commits and files listed above; the stage 05 SHA is the commit that adds this file |
| Changelog | `CHANGELOG.md` Unreleased, Added: covariant Lyapunov vectors entry |

Not changed: `_config/build-order.md` row M5 P3 (status by milestone); `src/precompile.jl` has no clv-ginelli workload; `_config/risks.md` (reviewed at the M5 exit, the next step). The arXiv copies of the three sources are in the session scratch folder only; copying them to `data/literature/clv/` is an owner decision.
