# Paper

Turn accepted study results into one manuscript per run. Each run is named by a slug, [paper], which prefixes every output file.

## Folder Map

```
paper/
├── CLAUDE.md          (you are here)
├── CONTEXT.md         (start here for task routing)
└── stages/
    ├── 01-outline/      (Outline)
    ├── 02-draft/        (Draft)
    ├── 03-review/       (Review and revise)
    └── 04-finalize/     (Finalize)
```

Shared references live in `/_config/`; onboarding is the repo-level `/setup/`.

## Triggers

| Keyword | Action |
|---------|--------|
| `setup` | Run `/setup/questionnaire.md` (repo level) |
| `status` | For every slug found in `stages/*/output/`, show which stages have that slug's artifact |
| `status <slug>` | Show the pipeline for one slug |

`status <slug>` renders, with COMPLETE when `output/[paper]-*` exists in that stage, else PENDING:

```
[01-outline]  -->  [02-draft]  -->  [03-review]  -->  [04-finalize]
```

## Routing

| Task | Go To |
|------|-------|
| Outline | `stages/01-outline/CONTEXT.md` |
| Draft | `stages/02-draft/CONTEXT.md` |
| Review and revise | `stages/03-review/CONTEXT.md` |
| Finalize | `stages/04-finalize/CONTEXT.md` |

## What to Load

| Task | Load These | Do NOT Load |
|------|-----------|-------------|
| Outline | `stages/01-outline/CONTEXT.md` and only its Inputs rows | module-build and literature-review outputs; study data in `data/` |
| Draft | `stages/02-draft/CONTEXT.md` and only its Inputs rows | module-build and literature-review outputs; study data in `data/` |
| Review and revise | `stages/03-review/CONTEXT.md` and only its Inputs rows | module-build and literature-review outputs; study data in `data/` |
| Finalize | `stages/04-finalize/CONTEXT.md` and only its Inputs rows | module-build and literature-review outputs; study data in `data/` |

## Stage Handoffs

Each stage writes to its own `output/`, file names prefixed by the slug. The next stage reads from there. Edit an output and the next stage picks up the edit.
