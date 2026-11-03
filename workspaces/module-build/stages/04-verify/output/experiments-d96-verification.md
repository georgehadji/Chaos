# Verification: experiments, Amendment D96 (stage 04)

Date: 2026-11-03, Julia 1.13.2, code at c41d085. Spec: `experiments-spec.md`, "Amendment D96" (R28, I38, I39, O27, criteria D96a to D96d, P116). Test plan: `experiments-test-plan.md`, "Amendment D96". Notes: `experiments-d96-impl-notes.md`.

Verify scripts: `verify/experiments-d96.jl` (library flow, and one `count_round` step) and `verify/experiments-d96-cli.sh` (the same campaign through `bin/chaoslab.jl run --round-log`). Default project, the environment of `run` (study stage 03). Raw output: `experiments-d96-verify-log.txt`.

Fixture: the inputs of the discovery `d96` (slug p3-abstention, the arm structure of the drafted A2 block: ref, ind, ind_pa1, d64, d64_pa1, chk; specs tg, fq1, pa1, bk1, all held by closure arms only; panel seed 777, 4 K, 46 orbits, 3 rounds, N_1 80000, N_c 20480, chk horizon 1000), copied to the scratchpad folder `d96v/fixture` with a SHA-256 list (`d96v/fixture-sha256.txt`, 21 files). Nothing under `data/` was read or written.

The runner is not patched: no Experiments function is redefined. `count_round` (Analysis, D95) does not yet read the runner's cell form `k<i>:<class>:r<round>` (mismatch M1, Amendment D97), so the scripts load the three scratch Analysis patches of the discovery (`p1_cells`, `p3a_missing`, `p3b_find`, each an `@eval Chaoslab.Analysis` block) for `count_round` only. O27 therefore verifies the runner; the count side is re-verified at D97 stage 04 without patches.

## Evidence per criterion

| Criterion | Evidence | Result |
|-----------|----------|--------|
| D96a (S) | `ex_d96.jl` D96a: 17 pass (first run complete, no table entry for tg and fq, no closure arm chunk; closure run lists both with the fields of R19) | pass |
| D96b (S), I38 | `ex_d96.jl` D96b: 3 pass, goldens of the indicator table taken at 8c3bac0 (before D96) for a spec in an each_round and a closure arm, and for a block without rounds | pass |
| D96c (S) | whole package suite: 1404618 pass, 0 fail, 1 broken (existing, weighted-birkhoff); quality 12 (Aqua, JET clean), layering 2; only test change: the D96a budget line (owner decision, test plan) | pass |
| I39 | the change is inside `_indicator_manifest` only; chunk keys, hashes and plan hash do not read the table (review item 2); D86 and D94 golden chunk and plan hashes in the suite pass | pass |
| D96d (V), O27 library | rounds 1 to 3: `run_campaign` status complete (26, 22, 20 chunks), indicator table empty before closure, `count_round` continue, continue, close (15, 30, 44 orbits); closure run 282 chunks, table lists bk1, fq1, pa1, tg; resume run 0 run, 350 skipped | pass |
| D96d (V), O27 CLI | `bin/chaoslab.jl run --round-log`: exit 0 at rounds 1 to 3, closure and resume; 350 chunk files in both stores; CLI round log byte-identical to the library log (SHA-256 4c96292f40538961...) | pass |
| O24 of D95, runner side | the runner accepts the log written by `count_round` at every round and at closure, library and CLI (rows above) | pass |

## Independent review

Fresh Sonnet subagent, given the spec section and the diff, read only. No CRITICAL, HIGH or MEDIUM finding.

| # | Finding | Severity | Disposition |
|---|---------|----------|-------------|
| 1 | `findfirst` over `plan.items` returns the same item as the old `first(...)`; specs without items are absent, not empty (R28, P116, I38) | none | confirmed |
| 2 | the cache of `indicator_manifest` keys on the plan hash, built from items and chunks: the closure plan gets a new key, a table without entries cannot be reused | none | confirmed |
| 3 | the only reader of the manifest's `indicators` table is `p2_cost` (`get` with a default); resume does not read the table; the manifest is rebuilt from the plan at each write | none | confirmed |
| 4 | before closure the manifest holds `"indicators" => {}` (an empty table), not an absent key | LOW | consistent with R28 (no spec listed); D96a reads it with a default; no change |
| 5 | no other per-spec `first(...)` over plan items in `src/Experiments` | none | confirmed |

## Audit

| Check | Result |
|-------|--------|
| Evidence | every criterion has a number (table above) |
| Review closed | no open CRITICAL or HIGH finding |

## Notes for the owner

- O27 is verified on the runner. Its count side uses the scratch Analysis patches until D97; the D97 stage 04 rerun of this campaign without patches closes O24 on the Analysis side.
- No plot: the amendment has no convergence or order quantity.
