# Literature Review

Find, verify and promote references for one topic per run. Each run is named by a slug, [topic], which prefixes every output file.

## Folder Map

```
literature-review/
├── CLAUDE.md          (you are here)
├── CONTEXT.md         (start here for task routing)
└── stages/
    ├── 01-scope/        (Plan a search)
    ├── 02-search/       (Run and screen the search)
    ├── 03-verify/       (Verify candidates)
    ├── 04-synthesize/   (Annotate papers (optional))
    └── 05-promote/      (Merge into the bibliography)
```

Shared references live in `/_config/`; onboarding is the repo-level `/setup/`.

## Triggers

| Keyword | Action |
|---------|--------|
| `setup` | Run `/setup/questionnaire.md` (repo level) |
| `status` | For every slug found in `stages/*/output/`, show which stages have that slug's artifact |
| `status <slug>` | Show the pipeline for one slug |

`status <slug>` renders, with COMPLETE when `output/[topic]-*` exists in that stage, else PENDING:

```
[01-scope]  -->  [02-search]  -->  [03-verify]  -->  [04-synthesize]  -->  [05-promote]
```

## Routing

| Task | Go To |
|------|-------|
| Plan a search | `stages/01-scope/CONTEXT.md` |
| Run and screen the search | `stages/02-search/CONTEXT.md` |
| Verify candidates | `stages/03-verify/CONTEXT.md` |
| Annotate papers (optional) | `stages/04-synthesize/CONTEXT.md` |
| Merge into the bibliography | `stages/05-promote/CONTEXT.md` |

## What to Load

| Task | Load These | Do NOT Load |
|------|-----------|-------------|
| Plan a search | `stages/01-scope/CONTEXT.md` and only its Inputs rows | other workspaces; `_config/catalogs/` except rows named in the plan |
| Run and screen the search | `stages/02-search/CONTEXT.md` and only its Inputs rows | other workspaces; `_config/catalogs/` except rows named in the plan |
| Verify candidates | `stages/03-verify/CONTEXT.md` and only its Inputs rows | other workspaces; `_config/catalogs/` except rows named in the plan |
| Annotate papers (optional) | `stages/04-synthesize/CONTEXT.md` and only its Inputs rows | other workspaces; `_config/catalogs/` except rows named in the plan |
| Merge into the bibliography | `stages/05-promote/CONTEXT.md` and only its Inputs rows | other workspaces; `_config/catalogs/` except rows named in the plan |

## Stage Handoffs

Each stage writes to its own `output/`, file names prefixed by the slug. The next stage reads from there. Edit an output and the next stage picks up the edit.
