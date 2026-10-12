# Coding Conventions

Layer 3, canonical. Applies to `src/`, `test/`, `bench/`, `bin/`, `scripts/`, `python/nbody/`.

## Code

No abstract type with one subtype; no registry beyond a `Dict{String,Any}` mapping TOML names to constructors; functions under 50 lines; files under 400.

No `@fastmath`. Accumulated logarithms use compensated (Kahan) summation. Precision is a type parameter, never a code fork.

## Tests

- Registration: a module's suite is `test/<module>/runtests.jl`; `test/runtests.jl` runs every such file automatically, after the quality (Aqua, JET) and layering checks.

- Invariant tests at random seeded points: det DT = 1, J^T Omega J = Omega for maps and one-step integrator maps; Jacobian vs complex-step / ForwardDiff.
- Integrator order from log-log slope; time reversibility; bounded energy error; modified energy conserved to higher order.
- Asymptotic-equivalence test (claim C1 as a test): all QR variants agree on lambda_1 for strongly chaotic orbits.
- Oracles: ChaosTools (lyapunov, gali), TaylorIntegration BigFloat spectra, REBOUND MEGNO and frequency analysis, nafflib/PyNAFF via HDF5.
- Cited regression values only from the verified bibliography (never from memory).
- Runtime health checks raise `NumericalHealthError`: NaN/Inf, orthogonality loss above tolerance, pairing residual above tolerance.
- Reproducibility: StableRNGs seeds, committed Manifest, provenance block in every store (git SHA, Julia and package versions, CPU, threads, precision, config hash); no `@fastmath`; compensated (Kahan) summation for accumulated logs.
- Quality gates: Aqua.jl, JET.jl, coverage >= 80%, CI on Windows + Ubuntu (Julia LTS + current). Property tests use plain Test.jl loops over seeded random points (Supposition.jl not adopted: last release 2025-11).
