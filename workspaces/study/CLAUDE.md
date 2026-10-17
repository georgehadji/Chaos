# Study

Run one preregistered numerical study per run. Each run is named by a slug, [study], which prefixes every output file.

## Folder Map

```
study/
├── CLAUDE.md          (you are here)
├── CONTEXT.md         (start here for task routing)
└── stages/
    ├── 01-design/       (Preregister)
    ├── 02-panel/        (Build the orbit panel)
    ├── 03-run/          (Run)
    ├── 04-analyze/      (Analyze)
    └── 05-report/       (Report results)
```

Shared references live in `/_config/`; onboarding is the repo-level `/setup/`.

## Triggers

| Keyword | Action |
|---------|--------|
| `setup` | Run `/setup/questionnaire.md` (repo level) |
| `status` | For every slug found in `stages/*/output/`, show which stages have that slug's artifact |
| `status <slug>` | Show the pipeline for one slug |

`status <slug>` renders, with COMPLETE when `output/[study]-*` exists in that stage, else PENDING:

```
[01-design]  -->  [02-panel]  -->  [03-run]  -->  [04-analyze]  -->  [05-report]
```

## Routing

The task table is `CONTEXT.md`, section "Task Routing" (one home; this file adds only what to load per task).

## What to Load

| Task | Load These | Do NOT Load |
|------|-----------|-------------|
| Preregister | `stages/01-design/CONTEXT.md` and only its Inputs rows | `workspaces/module-build/` outputs (Pattern 14: the library is used as a package, never via build notes) |
| Build the orbit panel | `stages/02-panel/CONTEXT.md` and only its Inputs rows | `workspaces/module-build/` outputs (Pattern 14: the library is used as a package, never via build notes) |
| Run | `stages/03-run/CONTEXT.md` and only its Inputs rows | `workspaces/module-build/` outputs (Pattern 14: the library is used as a package, never via build notes) |
| Analyze | `stages/04-analyze/CONTEXT.md` and only its Inputs rows | `workspaces/module-build/` outputs (Pattern 14: the library is used as a package, never via build notes) |
| Report results | `stages/05-report/CONTEXT.md` and only its Inputs rows | `workspaces/module-build/` outputs (Pattern 14: the library is used as a package, never via build notes) |

## Stage Handoffs

Each stage writes to its own `output/`, file names prefixed by the slug. The next stage reads from there. Edit an output and the next stage picks up the edit.
