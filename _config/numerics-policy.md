# Numerics Policy

Layer 3, canonical. Rules only; values live in `numerics-defaults.toml`.

## Rules

Default values (tau, tolerances, health-check thresholds, seed policy) live in `numerics-defaults.toml`, which `Experiments` reads at run time; this file never restates a number.

The map is defined on lifted coordinates; long runs store the state wrapped (both angle-type coordinates reduced each step with exact reduction modulo 2 pi) plus integer winding counts per coordinate, so the lifted orbit stays recoverable without losing resolution (decision D28), reset at each checkpoint with the integer offsets stored (decision D31); the raw lifted mode is kept for the wrapped-vs-lifted sensitivity study; renormalization interval `tau` is a recorded run parameter; deviation-vector initialization random-orthonormal with recorded seed; thresholds pre-registered per study, never hard-coded; precision recorded per run; two-trajectory perturbation size swept, not fixed.
