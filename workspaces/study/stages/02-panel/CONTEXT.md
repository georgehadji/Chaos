# Stage 02: Panel

Build and inspect the orbit panel. Skip for grid scans.

## Inputs

| Source | File/Location | Section/Scope | Why |
|--------|--------------|---------------|-----|
| Previous stage | `../01-design/output/[study]-preregistration.md` | "Panel" section | orbit classes to build |
| Reference | `references/panel-rounds.md` | Full file | only when the panel block has the D94 keys (rounds) |

## Process

1. Run `julia --project bin/chaoslab.jl panel <preregistration> --out output/[study]-panel.toml` (no AI): picks initial conditions per class and writes the panel file with its `.provenance.toml` (and, for a round panel, `.pools.toml`) beside it; run the audit.
2. Draw the portraits from the panel with a verify script, under `--project=verify` (the command draws none): `verify/experiments.jl` section PORTRAITS (p1), `verify/p3-abstention-portraits.jl <panel> <png>` (p3).
3. **[Checkpoint]** Present the portraits.
4. If a class is wrong, the human amends `[study]-preregistration.md` in 01's output (logged amendment, allowed only before 03-run; a study whose preregistration seals the panel hash list allows none after the panel build), then this stage reruns.

## Checkpoints

| After Step | Agent Presents | Human Decides |
|------------|---------------|---------------|
| 2 | Phase portraits with the chosen orbits marked (round panel: with the `round_cells` status of the panel manifest) | Orbits match their classes, or amend |

## Audit

| Check | Pass Condition |
|-------|---------------|
| Classes | Every orbit class of the preregistration has at least one orbit in the panel |
| Seeds | Every initial condition carries its recorded seed |

## Outputs

| Artifact | Location | Format |
|----------|----------|--------|
| Panel | `output/[study]-panel.toml` | Initial conditions, class labels, seeds (read by 03) |
| Portraits | `output/[study]-portraits.png` | Phase portraits (read by the human at the checkpoint) |
