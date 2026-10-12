# Shared References (Layer 3)

Canonical facts. Workspaces read these; nothing here points to a workspace. Stage CONTEXT files name the exact section to load.

| File | Holds |
|------|-------|
| `research-question.md` | Question, claims C1-C3, ground truth, test panel, outputs retained, paper sequence |
| `research-ideas.md` | Ranked idea backlog, prior art to check |
| `build-order.md` | Milestones M0-M8 and their slugs |
| `architecture.md` | Languages, module map, layering, repo shape, scalability |
| `numerics-policy.md` | Numerics rules |
| `numerics-defaults.toml` | Default numeric values (read by code) |
| `coding-conventions.md` | Code and test rules |
| `value-types.md` | What a study or paper can deliver |
| `publication.md` | Author, journals, citation style |
| `environment.md` | Compute target, cores, data root, precision backend |
| `decision-log.md` | Structural decisions, why, when to revisit |
| `risks.md` | Risk list, reviewed at every milestone exit |
| `catalogs/CONTEXT.md` | Methods, indicators, integrators, systems, module `Status` |
| `bibliography/CONTEXT.md` | Verified references, benchmarks, UNVERIFIED list |
| `setup-guides/` | `julia.md`, `rebound.md`, `latex.md` |

Gated writes: only literature-review/05 (bibliography, catalog `?` keys) and module-build/05 (catalog `Status`) change these files, after human approval.
