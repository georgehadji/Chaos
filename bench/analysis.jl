using BenchmarkTools, Chaoslab, StableRNGs

const SUITE = BenchmarkGroup()
const AN = Chaoslab.Analysis

# The bootstrap at the preregistered size order: 1e4 resamples of a set of 200 orbits (H2, H3), and a pointwise band over 201 recorded steps for a cell of 40.
const VALS = rand(StableRNG(20261012), 200)
SUITE["bootstrap_interval 200 orbits 1e4 resamples"] = @benchmarkable AN.bootstrap_interval($VALS, "bench"; seed = 20261012, resamples = 10_000, level = 0.95)
const MAT = rand(StableRNG(20261013), 40, 201)
SUITE["bootstrap_band 40 orbits 201 steps 1e3 resamples"] = @benchmarkable AN.bootstrap_band($MAT, "bench"; seed = 20261012, resamples = 1000, level = 0.95)

# Per-orbit work on a recorded grid of 201 steps: convergence time over the eight methods of V and the spread.
const STEPS = unique(round.(Int, 10 .^ range(0, 7; length = 201)))
const LAM = Dict{String, Vector{Float64}}(m => 0.3 .* (1 .+ 0.02 ./ STEPS) for m in ("L1", "L3", "L4", "L4b", "L5", "L5b", "L6", "L9"))
SUITE["orbit_convergence 8 methods"] = @benchmarkable AN.orbit_convergence($STEPS, $LAM, 0.3, 0.05)
SUITE["relative_spread 8 methods"] = @benchmarkable AN.relative_spread($LAM, 0.3)
