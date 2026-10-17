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

## Amendments
A change after the spec checkpoint is appended, never edited into the base sections: one section `## Amendment Dnn (date)` per decision, holding what changed and, when criteria change, a subsection `### Acceptance criteria (amendment)` with the new or replaced criteria by number; small owner corrections to an amendment are appended to it as "Corrections ..." paragraphs. Stage 04 reads the base "Acceptance criteria" and "Oracle" plus every amendment's criteria and corrections, not the whole file. Sub-criteria of an amendment are named by its decision ID and a letter (D44a, D44b); the letter IDs live in the spec only, not in the decision log.
