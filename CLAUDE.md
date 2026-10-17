# Chaos

Research on Lyapunov exponents, chaos indicators and symplectic integration, built with ICM (Interpretable Context Methodology). Library: Julia package `Chaoslab`. License: GPL-3.0-or-later.

## Folder Map

```
Chaos/
├── CLAUDE.md, CONTEXT.md   (Layer 0, Layer 1: you are here)
├── _core/                  (ICM conventions and templates, MIT)
├── _config/                (Layer 3: shared references, canonical facts)
├── setup/                  (onboarding questionnaire)
├── workspaces/
│   ├── workspace-builder/  (ICM builder; produced this structure)
│   ├── literature-review/  (one run per topic)
│   ├── module-build/       (one run per library module)
│   ├── study/              (one run per study)
│   └── paper/              (one run per paper, optional)
├── schema/                 (HDF5 run-store contract, Julia and Python)
├── src/ test/ bench/       (Julia package; built in milestone M0 onward; bin/ planned for M3)
├── ext/                    (CUDA package extension, D29)
├── verify/                 (stage-04 oracle and GPU scripts, own environment, D27)
├── scripts/                (ICM mechanics, Python stdlib)
├── python/nbody/           (REBOUND adapters, P6 only)
└── data/                   (run outputs and full texts, git-ignored)
```

## Rules

- Subagents run on Sonnet; mechanical subagent tasks (listing, formatting, status) on Haiku. Never Opus for subagents.
- Never cite from memory. Only keys in `_config/bibliography/references.bib`; never anything in `unverified.md`.
- Commits are fine at the stage that says so; push whenever needed, never force-push (D45); every PR or submission needs explicit human approval.
- No em dashes. Names in `workspaces/`, `_config/`, `setup/`, `schema/` are lowercase-with-hyphens; code and tooling files follow their ecosystem (`Project.toml`, `src/Chaoslab.jl`, `scripts/icm_validate.py`).
- Each fact has one home (see `_config/CONTEXT.md`). Point to it; do not copy it.
- Structural decisions go into `_config/decision-log.md`; review `_config/risks.md` at every milestone exit.

## Triggers

| Keyword | Action |
|---------|--------|
| `setup` | Run `setup/questionnaire.md` |
| `status` | `python scripts/icm_status.py`: one line per workspace (slugs, stages complete) plus catalog `Status` counts |
| `status <slug>` | Pipeline diagram for one topic, module, study or paper |

## Routing

| Task | Go To |
|------|-------|
| Find or verify literature | `workspaces/literature-review/CONTEXT.md` |
| Build or change a library module | `workspaces/module-build/CONTEXT.md` |
| Run a study | `workspaces/study/CONTEXT.md` |
| Write a paper | `workspaces/paper/CONTEXT.md` |
| Change the workspace structure | `workspaces/workspace-builder/CONTEXT.md` |
| Look up a shared fact | `_config/CONTEXT.md` |
