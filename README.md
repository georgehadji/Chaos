# Chaoslab

Julia library for Lyapunov exponents, chaos indicators and their numerical checks on symplectic maps. Several estimators run on the same orbits in Float64, Double64 or BigFloat, on the CPU or a CUDA GPU, and every run records the residuals that show whether it can be trusted.

## Install

Julia 1.11 or later:

```julia
using Pkg
Pkg.add(url = "https://github.com/georgehadji/Chaos", rev = "v0.1.2")
```

## Example

```julia
using Chaoslab, StaticArrays

sys = StandardMap(exact(Float64, 6))
res = propagate(sys, [SVector(0.1, 0.2)], 100_000; k = 2, tau = 1, seed = 1, reducers = (QRSpectrum(MGS()),))
res.finals[1][1].lambda    # [1.1348588476487835, -1.1348588476487835]
```

## Documentation

The user guide and the API reference build locally:

```sh
julia --project=docs -e "using Pkg; Pkg.instantiate()"
julia --project=docs docs/make.jl     # opens as docs/build/index.html
```

## Citing and license

Cite the tagged release you used: see `CITATION.cff`; each release is archived on Zenodo with a DOI. Changes per release: `CHANGELOG.md`. License: GPL-3.0-or-later (`LICENSE`).

The repository also holds the research workspaces that produced the library and its studies (`CLAUDE.md` is the entry point).
