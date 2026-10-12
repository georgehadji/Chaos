# Module Build

Build one verified Julia library module per run, test first.

## Task Routing

| Task Type | Go To | Description |
|-----------|-------|-------------|
| Write the math spec | `stages/01-spec/CONTEXT.md` | stage 01 |
| Write failing tests | `stages/02-tests/CONTEXT.md` | stage 02 |
| Write the code | `stages/03-implement/CONTEXT.md` | stage 03 |
| Verify against oracles | `stages/04-verify/CONTEXT.md` | stage 04 |
| Flip status and commit | `stages/05-integrate/CONTEXT.md` | stage 05 |

## Shared Resources

| Resource | Location | Contains |
|----------|----------|----------|
| Shared references | `/_config/CONTEXT.md` | Canonical facts, catalogs, bibliography |
| Build order | `/_config/build-order.md` | Which slugs to run, in which milestone |
