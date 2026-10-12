# Spec Format

Every `[module]-spec.md` has these sections, in this order. WHAT, not HOW: no Julia type, function or file names (ICM Pattern 10).

## Scope
Catalog rows covered (IDs), prerequisites (verified modules).

## Equations
Each equation with its source key and equation number.

## Algorithm
Math or pseudocode. Coefficient tables named by paper and table number, not copied.

## Invariants
Properties that must hold at every step (e.g. det = 1, symplecticity, sum of exponents).

## Compute
Which work is data-parallel (independent orbits, grid points) and so a GPU candidate; what stays on CPU and why; memory per orbit; the CPU-vs-GPU agreement tolerance. Policy: `/_config/architecture.md`, section "Compute".

## Oracle
Independent check: package and settings, or analytic result.

## Acceptance criteria
Numbered, each with a measurable pass condition and a tolerance.

## Sources
Verified keys used, with locations.
