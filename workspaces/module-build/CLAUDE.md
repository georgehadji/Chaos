# Module Build

Build one verified Julia library module per run, test first. Each run is named by a slug, [module], which prefixes every output file.

## Folder Map

```
module-build/
├── CLAUDE.md          (you are here)
├── CONTEXT.md         (start here for task routing)
└── stages/
    ├── 01-spec/         (Write the math spec)
    ├── 02-tests/        (Write failing tests)
    ├── 03-implement/    (Write the code)
    ├── 04-verify/       (Verify against oracles)
    └── 05-integrate/    (Flip status and commit)
```

Shared references live in `/_config/`; onboarding is the repo-level `/setup/`.

## Triggers

| Keyword | Action |
|---------|--------|
| `setup` | Run `/setup/questionnaire.md` (repo level) |
| `status` | For every slug found in `stages/*/output/`, show which stages have that slug's artifact |
| `status <slug>` | Show the pipeline for one slug |

`status <slug>` renders, with COMPLETE when `output/[module]-*` exists in that stage, else PENDING:

```
[01-spec]  -->  [02-tests]  -->  [03-implement]  -->  [04-verify]  -->  [05-integrate]
```

## Routing

| Task | Go To |
|------|-------|
| Write the math spec | `stages/01-spec/CONTEXT.md` |
| Write failing tests | `stages/02-tests/CONTEXT.md` |
| Write the code | `stages/03-implement/CONTEXT.md` |
| Verify against oracles | `stages/04-verify/CONTEXT.md` |
| Flip status and commit | `stages/05-integrate/CONTEXT.md` |

## What to Load

| Task | Load These | Do NOT Load |
|------|-----------|-------------|
| Write the math spec | `stages/01-spec/CONTEXT.md` and only its Inputs rows | study, paper and literature-review outputs; catalog rows of other modules |
| Write failing tests | `stages/02-tests/CONTEXT.md` and only its Inputs rows | study, paper and literature-review outputs; catalog rows of other modules |
| Write the code | `stages/03-implement/CONTEXT.md` and only its Inputs rows | study, paper and literature-review outputs; catalog rows of other modules |
| Verify against oracles | `stages/04-verify/CONTEXT.md` and only its Inputs rows | study, paper and literature-review outputs; catalog rows of other modules |
| Flip status and commit | `stages/05-integrate/CONTEXT.md` and only its Inputs rows | study, paper and literature-review outputs; catalog rows of other modules |

## Stage Handoffs

Each stage writes to its own `output/`, file names prefixed by the slug. The next stage reads from there. Edit an output and the next stage picks up the edit.
