# Chaoslab

Chaoslab is a Julia library for Lyapunov exponents, chaos indicators and their numerical checks on symplectic maps. It computes finite-time Lyapunov exponents with several algorithms on the same orbits, in Float64, Double64 or BigFloat, on the CPU or a CUDA GPU, and records the residuals that show whether a run can be trusted.

Version 0.1 covers the standard map and linear symplectic maps, the Benettin, QR-family, singular-value and two-trajectory estimators, run diagnostics, and the campaign and analysis layers used by the P1 benchmark study.

## Installation

Chaoslab needs Julia 1.11 or later. It is not in the General registry; add it from the repository:

```julia
using Pkg
Pkg.add(url = "https://github.com/georgehadji/Chaos", rev = "v0.1.2")
```

The GPU path is a package extension: it loads when `CUDA` is installed and imported. The figures of the analysis layer are an extension on `CairoMakie`.

## A first exponent

The largest Lyapunov exponent of the standard map at K = 6 from one orbit, with the Benettin reducer. `exact` builds the parameter in the element type of the state, so the same code runs in Double64 or BigFloat.

```@example quickstart
using Chaoslab, StaticArrays

sys = StandardMap(exact(Float64, 6))
res = propagate(sys, [SVector(0.1, 0.2)], 100_000; k = 1, tau = 10, seed = 1, reducers = (Benettin(),))
res.finals[1][1].lambda[1]
```

The full spectrum with a QR reducer (modified Gram-Schmidt) and two deviation vectors. For an area-preserving map the two exponents have opposite signs and sum to zero up to rounding.

```@example quickstart
res = propagate(sys, [SVector(0.1, 0.2)], 100_000; k = 2, tau = 1, seed = 1, reducers = (QRSpectrum(MGS()),))
res.finals[1][1].lambda
```

`res.provenance` records the backend that ran, the chunk size and any fallback reason. A non-finite state or a failed orthonormalisation throws [`NumericalHealthError`](@ref) with the step and the orbit.

## Layers

| Layer | Holds |
|-------|-------|
| `Foundation` | health errors, compensated sums, backend resolution, precision |
| `Systems` | maps, tangent maps, symplectic residuals |
| `Tangent` | the propagation kernel and the reducer protocol |
| `Diagnostics` | reducers that check a run |
| `Lyapunov` | the exponent estimators |
| `Indicators` | chaos indicators as reducers (SALI) |
| `Experiments` | preregistered campaigns and the HDF5 run store |
| `Analysis` | analysis of a finished study run |

The [API](api.md) page lists every exported name with its docstring.

## Citing

Cite the tagged release you used. The repository holds `CITATION.cff`, and each release is archived on Zenodo with a DOI. Chaoslab is licensed under GPL-3.0-or-later.
