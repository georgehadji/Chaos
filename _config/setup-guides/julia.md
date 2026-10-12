# Setup: Julia

Julia is the language of the core library `Chaoslab`. Needed by module-build and study.

## Install

1. Install juliaup (the Julia version manager) from https://julialang.org/install/ . On Windows it is offered through the Microsoft Store and winget.
2. Open a new terminal. Run `juliaup status`; the default channel should be the current stable release.

## Verify

- `julia --version` prints the version.
- From the repo root: `julia --project -e "using Pkg; Pkg.instantiate(); Pkg.test()"` installs the pinned packages from `Manifest.toml` and runs the tests.

## How this repo uses it

- Library code in `src/`, tests in `test/`, benchmarks in `bench/`, entry point `bin/chaoslab.jl`.
- Record the Julia version in every run manifest (done automatically by `Experiments`).
