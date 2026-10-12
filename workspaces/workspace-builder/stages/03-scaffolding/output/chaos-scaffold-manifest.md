# Chaos: Scaffold Manifest (ICM Stage 03)

The scaffold is at the repo root, not in this folder: the workspaces, `_config/` and root files must live where they run. This file records what was generated, the audit, and known gaps. Stages 04 and 05 read the scaffold at the repo root. Date: 2026-10-12.

## Generated

| Path | What | How |
|------|------|-----|
| `CLAUDE.md`, `CONTEXT.md` | Layer 0 and Layer 1: folder map, rules, triggers, routing | written |
| `LICENSE` | GPL-3.0 text (SPDX `GPL-3.0-or-later` goes into `Project.toml` at M0) | `gh api licenses/gpl-3.0` |
| `.gitignore` | `data/`, caches, venvs, coverage files | written |
| `setup/` | empty; Stage 04 writes `questionnaire.md` | `.gitkeep` |
| `_config/` | 10 files plus `catalogs/`, `bibliography/`, `setup-guides/`, each folder with CONTEXT.md | promoted from Stage 01 outputs, headings renamed to the contract names |
| `_config/bibliography/` | `references.bib` plus 8 section files and `conventions.md` | split at the Stage 01 file's sections |
| `_config/catalogs/` | 4 catalogs with a `Status` column, all `planned` | Stage 01 catalogs |
| `schema/run-store.md` | HDF5 layout v1 | written |
| `workspaces/{literature-review,module-build,study,paper}/` | CLAUDE.md, CONTEXT.md, 19 stages each with CONTEXT.md, `output/.gitkeep`, `references/.gitkeep` | stage CONTEXT files generated from the Stage 02 contract files by script (contract fidelity by construction) |
| `module-build/.../01-spec/references/spec-format.md`, `study/.../01-design/references/preregistration-format.md` | required sections that stage Inputs quote | written |

Stage 01 files now carry a one-line note naming their canonical `_config/` copy.

## Contract changes made during scaffolding

- study/01 and 03: compute limits moved from inline placeholders to `/_config/environment.md` (placeholders only in Inputs values or Layer 3 files).
- module-build/01: `spec-format.md` and `bibliography/benchmarks.md` added to Inputs; module-build/03 drops the "rules line" (now `coding-conventions.md`).
- Setup guides added to module-build/02, study/03, paper/04 Inputs.
- Canonical-sources table in `chaos-stage-contracts.md` gains environment, setup guides, spec and preregistration formats.

## Placeholders (for Stage 04)

Only genuinely unknown values: `{{AUTHOR_NAME}}`, `{{AFFILIATION}}`, `{{ORCID}}`, `{{TARGET_JOURNALS}}`, `{{CITATION_STYLE}}` (in `_config/publication.md`); `{{COMPUTE_TARGET}}`, `{{MAX_CORES}}`, `{{DATA_ROOT}}`, `{{HIGH_PRECISION_BACKEND}}` (in `_config/environment.md`). Decided values are hardcoded: repo, license, subagent models, Julia via juliaup, default precision Float64.

## Audit

| Check | Result |
|-------|--------|
| Folder structure | PASS: every stage has CONTEXT.md, `output/.gitkeep`, `references/` |
| Contract fidelity | PASS: stage CONTEXT files generated from the contract files |
| Placeholder syntax | PASS: `{{SCREAMING_SNAKE_CASE}}`, none in CLAUDE.md or root CONTEXT.md, none outside Inputs in stage files |
| .gitkeep coverage | PASS |
| CONTEXT.md size | PASS: largest stage file 45 lines |
| Reference size | PASS: every reference file under 200 lines |
| Naming | PASS: lowercase-with-hyphens (leading `_` allowed, as in `_core`) |
| Paths | PASS: every `/_config/` and `/schema/` path in an Inputs table exists |
| Sections | PASS: 23 quoted section names resolve to headings; a planted bad name was caught |
| Em dashes | PASS: none |

## Pending M0 (referenced, not yet built)

`/src/`, `/test/runtests.jl`, `/bench/`, `/bin/chaoslab.jl`, `/scripts/bibcheck.py`, `/scripts/promote.py`, `/scripts/icm_status.py`, `/scripts/icm_validate.py`, `/python/nbody/`, `/CHANGELOG.md`, `Project.toml`, `Manifest.toml`, CI. Design for `bibcheck.py` (proven while building the bibliography): stdlib only (`urllib`, `json`, a thread pool of 6), fetch `https://api.crossref.org/works/<doi>` with a `mailto` User-Agent taken from an environment variable, compare normalized title, first-author family name, author count, volume and pages, print one line per mismatch.
