# Implementation notes: analysis, Amendment D99 (module-build stage 03)

Julia 1.13.2, 4 threads. Spec, tests, `schema/`, `_config/`, `data/`, `verify/` and every other module untouched. No dependency added. The committed p3 outputs are not regenerated (stage 04, O31). This file is read by stage 04.

## Result

| Run | Pass | Fail | Error | Broken | Time |
|-----|------|------|-------|--------|------|
| RED (stage 02, `an_p3t_d99.jl`) | 48 | 25 | 0 | 0 | |
| GREEN, p3 subset (driver `d99/driver.jl`, sum of testset times) | 944 | 0 | 0 | 0 | 3 min 8 s |
| GREEN, whole package suite (`Pkg.test()`, 4 threads) | 1 404 893 | 0 | 0 | 1 | 25 min 10 s |

Subset per file: fixture 41, shape 52, stats 327, params 49, log 45, refusals 5, dispatch 22, pipeline 111, d97 109, d97r 32, d99 81, p2 dispatch 33, p3 dispatch 37. `an_p3_dispatch.jl` calls `p2_reslug`, defined in `an_p2_dispatch.jl`; the driver runs the p2 file first, as `test/analysis/runtests.jl` does.

The broken result is the existing one (weighted-birkhoff). Quality 12 of 12 (Aqua, JET), layering 2 of 2, analysis 2 751, experiments 9 264.

Logs (scratchpad `d99/`): `d99-green.log`, `full.log`.

## Change (all in `src/Analysis/`, line counts of `git diff`)

| File | +/- | Change |
|------|-----|--------|
| `p3_read.jl` | 2/1 | E86: `_p3_inputs` defaults the output folder of `calibrate` to `stages/04-analyze/output`, the folder count-round and analyze read the record from; an explicit folder is used as given. |
| `p3t_tables.jl` | 12/6 | E83: `_p3t_own3` gives each row of HA item 3 the status of its own condition (point at least delta_a, lower bound above 0, an undefined bound fails); the T5 decision line states the status of item 3 and the rule sentence. E84: the T3 introduction names the omission a deviation of the analysis. E85: T2 row D has no target (`n.a.`); T8 has the line `Target failures n_f: ...; achieved failures: ...` after the `HC:` line. |
| `p3t_json.jl` | 2/2 | E84: the P120 text names the deviation of the analysis, the arm chk, N_B and the three precisions, and no study amendment. E85: `hc.target_failures` holds n_f. |

Largest changed function: `_p3t_t3` (21 lines). All files stay under 400 lines (`p3t_tables.jl` 197). The files keep their CRLF working-tree endings. No em dash, no `@fastmath`. No decision moves: `_p3t_own3` only labels rows; the decision of HA is still the one of `p3t_stats.jl` (I52).

## Readings

1. Decision line wording: "Status of item 3: fails" (or "ok"). The test reads "item 3", then at most two words, then "fails".
2. `_p3t_own3` compares exact rationals (`p3_rat` of the decimal text of delta_a and of a text point), the same comparison the decision rule makes.
3. The machine file gains one key, `hc.target_failures`, at the end of the `hc` block; the key order of the other blocks is unchanged.

## Process deviation

The implementation was written directly in this session, not by a Sonnet subagent as the stage asks.

## Not done here

Coverage of the changed functions was not measured by a tool. No benchmark change. D99g (chain on the D96 fixture), D99h (O31 before and after on the committed p3 inputs, in empty scratch folders) and D99i (p1 and p2 goldens) are stage 04. Nothing committed.
