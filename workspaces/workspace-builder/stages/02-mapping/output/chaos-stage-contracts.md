# Chaos: Stage Contracts (ICM Stage 02 Mapping)

Output of `workspace-builder/stages/02-mapping`; input to Stage 03 (scaffolding). Status: DRAFT, waiting at the Stage 02 checkpoint. Paths starting with `/` are from the repo root.

## Files

| File | Holds |
|------|-------|
| `chaos-stage-contracts.md` | This index: graph, canonical sources, handoffs, status rule, audit |
| `chaos-contracts-literature-review.md` | 01-scope, 02-search, 03-verify, 04-synthesize, 05-promote |
| `chaos-contracts-module-build.md` | 01-spec, 02-tests, 03-implement, 04-verify, 05-integrate |
| `chaos-contracts-study.md` | 01-design, 02-panel, 03-run, 04-analyze, 05-report |
| `chaos-contracts-paper.md` | 01-outline, 02-draft, 03-review, 04-finalize |

## 1. Dependency graph

Arrows point from reader to what it reads (one-way). `==>` marks a gated promotion write: the final stage of a run writes into `/_config/` only after the human approves the diff.

```
                 +-----------------------------------------------+
                 | /_config/  (Layer 3, references no workspace) |
                 |  research-question  research-ideas            |
                 |  build-order  architecture  numerics-*        |
                 |  coding-conventions  value-types  publication |
                 |  catalogs/ (Status)   bibliography/           |
                 +-----------------------------------------------+
                    ^   ^ reads            ^ reads          ^ reads
   literature-review|   |                  |                |
   01 -> 02 -> 03 -> 04 -> 05 ==> bibliography/, catalog `?` keys
                        |
   module-build         |
   01 -> 02 -> 03 -> 04 -> 05 ==> catalogs/ Status = verified
          |     |
          v     v
       /test/  /src/  (Julia package Chaoslab, entry point /bin/chaoslab.jl)
                 ^
                 | runs as package          /python/nbody/ --> /data/ (HDF5)
   study         |                                               ^
   01 -> 02 -> 03 -> 04 -> 05      03-run writes /data/[study]/ -+
                      |     |      per /schema/run-store.md
                      v     v
   paper  01 -> 02 -> 03 -> 04     (reads study 04 figures, 05 results)
```

Skip edges inside a workspace (each appears in an Inputs table): lit-review 03 -> 05; module-build 01 -> 03, 01 -> 04, 03 -> 04 (code); study 01 -> 03, 01 -> 04, 03 -> 05; paper 01 -> 02 reads study outputs directly.

## 2. Cross-workspace handoffs

| From | To | What | How |
|------|----|------|-----|
| literature-review/05 | `/_config/bibliography/` | verified entries, notes, benchmarks, UNVERIFIED | `scripts/promote.py`, gated write |
| literature-review/05 | `/_config/catalogs/` | `?` keys that now resolve | gated write |
| module-build/02, 03 | `/test/`, `/src/`, `/bench/` | tests, code, budgets | files |
| module-build/05 | `/_config/catalogs/` | `Status` = `verified` | gated write |
| study/01, 03 | `/src/` | library | activated Julia package only (Pattern 14) |
| study/03 | `/data/[study]/` | run store | `bin/chaoslab.jl run`; N-body via `/python/nbody/` |
| paper/01, 02 | study 04 and 05 `output/` | figures, tables, results | files |

No workspace reads a module-build `output/` folder except module-build itself.

## 3. Canonical sources (one home per fact)

| Fact | Home | Others point there |
|------|------|--------------------|
| Question, claims C1-C3, test panel, paper sequence | `/_config/research-question.md` | study/01, paper/01, lit-review/01 |
| Idea backlog and prior art to check | `/_config/research-ideas.md` | lit-review/01, study/01 |
| Milestones and module list | `/_config/build-order.md` | catalogs `Phase` column, module-build/01 |
| Languages, module map, layering | `/_config/architecture.md` | module-build/01, 03 |
| Numerics rules | `/_config/numerics-policy.md` | module-build/01, 02 |
| Numeric default values (tau, tolerances, health thresholds, seed policy) | `/_config/numerics-defaults.toml` (read by code) | numerics-policy points here |
| Code rules (sizes, style, tests) | `/_config/coding-conventions.md` | architecture keeps a pointer |
| Methods, systems, module status | `/_config/catalogs/*.md` (`Status` column) | everything |
| Citation records | `/_config/bibliography/references.bib` | all stages |
| Key format, verification classes | `/_config/bibliography/conventions.md` | lit-review/03 |
| Benchmark values | `/_config/bibliography/benchmarks.md` | module-build/01, 02 |
| Coefficient tables | the cited paper (spec names paper and table) | code implements, spec points |
| Module equations, acceptance tolerances | `[module]-spec.md` | tests, verification |
| Study hypotheses and thresholds | `[study]-preregistration.md` | 03, 04, 05 |
| HDF5 run layout | `/schema/run-store.md` | Julia `Experiments`, `/python/nbody/` |
| Value types | `/_config/value-types.md` | study/01, paper/01 |
| Journals, citation style, author data | `/_config/publication.md` | paper/01, 04 |
| Subagent model policy | `/CLAUDE.md` | every stage names the model per step |
| ICM rules | `/_core/CONVENTIONS.md` | builder, validation |

## 4. Status rule for multi-run workspaces

The ICM default (non-empty `output/` = COMPLETE) breaks after the first run. Each workspace `CLAUDE.md` defines instead:

- `status`: list every slug found in that workspace's `output/` file names; for each slug, a stage is COMPLETE when its output artifact for that slug exists.
- `status <slug>`: the pipeline diagram for one slug.
- Root `status`: one line per workspace with slug counts, plus catalog `Status` counts (planned, spec, tests, verified).

## 5. Agent model per step

Sonnet subagents: search, fixes against fetched records, paper reading, specs, tests, code, independent review, findings, drafting. Haiku: mechanical file listing, formatting, status rendering. No AI: `scripts/*.py`, `bin/chaoslab.jl`. Named in each Process step.

## 6. Reconciliations applied to Stage 01 outputs

1. Study commands moved from `scripts/*.py` to the Julia entry point `bin/chaoslab.jl` (workflow map Sections 2, 3, 4.3; architecture repo shape). `scripts/` stays Python stdlib for ICM mechanics.
2. Build order gets a Layer 3 home, `/_config/build-order.md`; the catalog `Phase` pointer now names it.
3. Numeric defaults live in `/_config/numerics-defaults.toml`; the numerics policy states rules only.
4. Research ideas get a home, `/_config/research-ideas.md` (row added to the companion table).

To do at Stage 03 (scaffolding), not now: split the bibliography into section files plus `benchmarks.md`, `unverified.md`, `conventions.md`; add a `Status` column to every catalog; move the code rules from the architecture into `coding-conventions.md`; write `value-types.md` and `publication.md` with placeholders.

## 7. Audit (run on this draft)

| Check | Result |
|-------|--------|
| No circular references | PASS. Reads flow from workspaces to `/_config/`, `/src/`, `/schema/` and upstream study outputs only; `/_config/` references nothing. The two write-backs into `/_config/` are declared, gated promotions, not references. |
| Output consumption | PASS. Every output names a reader. Human-gate artifacts: `[topic]-verify-report.md` (also read by lit-review/05), `[study]-portraits.png`, `[paper]-review-response.md`. Final deliverables: promoted bibliography, catalog status flips, `[study]-results.md`, `[paper]/` manuscript. Logs (`[topic]-promote-log.md`, `[module]-integration-log.md`) are read by `status`. |
| Contract completeness | PASS. All 19 stages have non-empty Inputs, Process, Outputs. |
| Canonical sources | PASS after reconciliations 1-4. |
