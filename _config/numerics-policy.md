# Numerics Policy

Layer 3, canonical. Rules only; values live in `numerics-defaults.toml`.

## Rules

Default values (tau, tolerances, health-check thresholds, seed policy) live in `numerics-defaults.toml`, which `Experiments` reads at run time; this file never restates a number.

Lifted coordinates internally, wrapping only for display and FMA windows; renormalization interval `tau` is a recorded run parameter; deviation-vector initialization random-orthonormal with recorded seed; thresholds pre-registered per study, never hard-coded; precision recorded per run; two-trajectory perturbation size swept, not fixed.
