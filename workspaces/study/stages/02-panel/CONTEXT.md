# Stage 02: Panel

Build and inspect the orbit panel. Skip for grid scans.

## Inputs

| Source | File/Location | Section/Scope | Why |
|--------|--------------|---------------|-----|
| Previous stage | `../01-design/output/[study]-preregistration.md` | "Panel" section | orbit classes to build |

## Process

1. Run `julia --project bin/chaoslab.jl panel` with the preregistration (no AI): picks initial conditions per class, writes the panel file and phase portraits.
2. **[Checkpoint]** Present the portraits.
3. If a class is wrong, the human amends `[study]-preregistration.md` in 01's output (logged amendment, allowed only before 03-run), then this stage reruns.

## Checkpoints

| After Step | Agent Presents | Human Decides |
|------------|---------------|---------------|
| 1 | Phase portraits with the chosen orbits marked | Orbits match their classes, or amend |

## Outputs

| Artifact | Location | Format |
|----------|----------|--------|
| Panel | `output/[study]-panel.toml` | Initial conditions, class labels, seeds (read by 03) |
| Portraits | `output/[study]-portraits.png` | Phase portraits (read by the human at the checkpoint) |
