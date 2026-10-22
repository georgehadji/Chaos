# Integration log: rli

Milestone M4 (D65), on `main`. Stage 05 of `workspaces/module-build`.

## Catalog rows flipped

| File | Row or section | Change |
|------|----------------|--------|
| `_config/catalogs/indicators.md` | I9 (RLI) | planned -> verified; main-issue column states the separation law (ordered proportional to Dx, chaotic independent, 10.5 orders apart at Dx = 1e-12, n = 1e4), the pair-map construction with the same xi_0 (D75) and the label position at K = 0.5 (rli R4); SandorMaffione2016 and Maffione2011a added to the sources, Sandor2004 marked not read |

Evidence: `stages/04-verify/output/rli-verification.md`, accepted at the stage 04 checkpoint (2026-10-22) with Amendment D77: criteria 1 to 15 pass on CPU and on the local GPU in the D68 form (Float64; Double64 kernel fallback allowed); criterion 4(a) O(Dx) law reported, not gated (D77 item 1); independent review R1 (MEDIUM) and R2 to R5 (LOW) fixed in code, none open. Decision rows: D65, D66, D69, D71, D72, D75, D76, D77.

Open: none for the module. Recorded limits: wrapped state required (R8), Dx = 1e-12 default in Float64 (R6), caller-side frame check for xi_0 (`src/Indicators/CONTEXT.md`, Limits). Sandor2004 not read (spec F4, P9): the owner may supply the PDF. Follow-ups for Experiments as for fli-ofli and megno (`src/Experiments/CONTEXT.md`, "Obligations from the Indicators layer"): build the pair with `pair_start` and the block frame, run `check_frame`, record Dx, Dx_r, the direction rule and the chain.

## Commits

- bf5f0c4 module-build rli stage 05: catalog I9 verified; spec notes (R1, R4, R6, R8, cost), Indicators limits (pairing check, caller-side check_frame for xi_0)
- 69936a1 module-build rli stage 04: verification accepted (criteria 1 to 15, CPU and GPU in the D68 form); Amendment D77: pairing check, split DL, review R1 to R5 fixed, criterion 4(a) law reported
- 01b3de6 Indicators CONTEXT: spec pointers as one list (80-line limit)
- 67491cb module-build rli stage 03: PairMap in Systems, RLI reducer in Indicators (GREEN, coverage 98.7%)
- b03aab5 module-build rli stage 02: RED tests (33 testsets) and test plan; spec Amendment D76
- 8685a8d module-build rli stage 01: spec approved (D75)

## Files changed by the module

Short hashes of the rli commits that touched each file. Stage paths are relative to `workspaces/module-build/`.

- `CHANGELOG.md` (bf5f0c4, 67491cb)
- `_config/catalogs/indicators.md` (bf5f0c4)
- `_config/decision-log.md` (69936a1, b03aab5, 8685a8d)
- `bench/rli.jl` (67491cb)
- `docs/src/api/indicators.md` (67491cb)
- `docs/src/api/systems.md` (67491cb)
- `src/Indicators/CONTEXT.md` (bf5f0c4, 01b3de6, 67491cb)
- `src/Indicators/Indicators.jl` (67491cb)
- `src/Indicators/rli.jl` (69936a1, 67491cb)
- `src/Systems/CONTEXT.md` (67491cb)
- `src/Systems/Systems.jl` (67491cb)
- `src/Systems/pair.jl` (69936a1, 67491cb)
- `src/precompile.jl` (67491cb)
- `test/rli/rl_ac01_03.jl` (b03aab5)
- `test/rli/rl_ac04_07.jl` (69936a1, b03aab5)
- `test/rli/rl_ac08_10.jl` (b03aab5)
- `test/rli/rl_ac11_13.jl` (69936a1, b03aab5)
- `test/rli/rl_api.jl` (b03aab5)
- `test/rli/rl_helpers.jl` (b03aab5)
- `test/rli/rl_pair.jl` (b03aab5)
- `test/rli/rl_records.jl` (b03aab5)
- `test/rli/runtests.jl` (b03aab5)
- `verify/gpu.jl` (69936a1)
- `verify/rli.jl` (69936a1)
- `stages/01-spec/output/rli-spec.md` (bf5f0c4, 69936a1, b03aab5, 8685a8d)
- `stages/02-tests/output/rli-test-plan.md` (b03aab5)
- `stages/03-implement/output/rli-impl-notes.md` (67491cb)
- `stages/04-verify/output/gpu-rli-verify-log.txt` (69936a1)
- `stages/04-verify/output/rli-chaostools.png` (69936a1)
- `stages/04-verify/output/rli-growth-laws.png` (69936a1)
- `stages/04-verify/output/rli-verification.md` (69936a1)
- `stages/04-verify/output/rli-verify-log.txt` (69936a1)
