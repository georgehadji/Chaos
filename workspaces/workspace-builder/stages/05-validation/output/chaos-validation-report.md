# Chaos: Validation Report (ICM Stage 05)

Scope: the scaffold at the repo root (see `../../03-scaffolding/output/chaos-scaffold-manifest.md`) and `/setup/questionnaire.md`. The vendored `_core/` and `workspaces/workspace-builder/` are out of scope. Date: 2026-10-12. Checks run by script; results below are from the final run after fixes.

## Results

| # | Check | Result | Notes |
|---|-------|--------|-------|
| 1 | Cross-reference integrity | PASS | Every `/_config/`, `/schema/` and `references/` path in a stage Inputs table exists. Paths built in M0 are listed below as known gaps, not failures. |
| 2 | No circular dependencies | PASS | 25 read edges between 19 stages; topological sort succeeds; no stage reads a later stage of its own workspace. `_config/` and `schema/` contain no path into a workspace. |
| 3 | Placeholder coverage | PASS | 9 placeholders in the repo, all covered by the questionnaire (8 asked, `{{MAX_CORES}}` derived); every question names the file holding its placeholder. |
| 4 | Conditional sections | PASS | None used. Q9 removes the whole `workspaces/paper/` folder instead. |
| 5 | Stage handoff chain | PASS | Every artifact a stage reads is listed in its producer's Outputs, including the cross-workspace reads of paper from study. |
| 6 | CONTEXT.md purity | PASS | Stage files hold only title, one sentence, and Inputs, Process, Checkpoints, Audit, Outputs. Workspace and `_config/` CONTEXT files hold routing tables. |
| 7 | Checkpoints in creative stages | PASS | All 11 creative stages have a checkpoint; every checkpoint names a valid step number. |
| 8 | Audits in creative and build stages | PASS | Every creative stage and module-build/03 has an audit, and its Process runs the audit before writing output. |
| 9 | Contract purity in spec stages | PASS | `spec-format.md` defines WHAT (equations, invariants, oracle, acceptance criteria); module-build/01's audit rejects Julia type, function or file names. |
| 10 | Line counts | PASS | Largest stage CONTEXT.md 45 lines (limit 80); every reference file under 200 lines. |
| 11 | Naming | PASS | lowercase-with-hyphens, leading `_` allowed (as in `_core`); zero-padded stage folders; every `output/` has `.gitkeep`. |
| 12 | Tool prerequisites | PASS | Guides in `_config/setup-guides/` for Julia, Python with uv and git, the REBOUND stack, LaTeX; each has Install and Verify. LaTeX is optional through Q9. The REBOUND stack is optional by decision D11 (only if P5 passes), so setup does not ask about it. |
| 13 | Quality scan | PASS | No em dashes; terms such as SALI, MEGNO, FMA are defined in the catalogs they appear in. |

## Issues found and fixed

| Issue | Fix |
|-------|-----|
| paper/02-draft (creative) had no checkpoint | Added: after the audit, present the draft and any untraced number or claim |
| module-build/04, study/03, paper/03 never said when their audit runs | Process steps now run the audit before writing output |
| `_config/environment.md` listed the stages that read it; `_config/setup-guides/latex.md` named a workspace output path | Consumer list removed; path replaced with a plain description (Pattern 3) |
| No setup guide for Python and git; REBOUND guide used `python -m venv`, which fails on the main machine (no `venv` module) | Added `python-and-git.md`; REBOUND guide now uses `uv venv` |

All fixes were made in the Stage 02 contract files first, then the stage CONTEXT files were regenerated, so contracts and stages still match.

## Accepted exception

`/setup/questionnaire.md` names `workspaces/paper/` because Q9 deletes that folder. Onboarding acts on the workspace by design, so this is not a reference into a stage.

## Known gaps (built in milestone M0)

`/src/`, `/test/runtests.jl`, `/bench/`, `/bin/chaoslab.jl`, `/scripts/bibcheck.py`, `/scripts/promote.py`, `/scripts/icm_status.py`, `/scripts/icm_validate.py`, `/python/nbody/`, `/CHANGELOG.md`, `Project.toml`, `Manifest.toml`, CI. `scripts/icm_validate.py` should absorb the checks used here.

## Open item outside the scaffold

celmech 1.5.8 and REBOUNDx 5.1.0 fail to build with MSVC on native Windows (gcc-only flags, C99 complex arithmetic, variable-length arrays). A native Windows build is in progress; the outcome goes to `_config/risks.md` (K2) and `_config/setup-guides/rebound.md`.
