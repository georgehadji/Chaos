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

# Study p2-indicators (Amendment D89): the bootstrap of H1 and H2 at the sizes of the real study (W2 214 orbits, S2 448 orbits, 1e4 resamples, seven indicators), the
# agreement of the 28 pairs of I_8 on 214 orbits and the AUC of two classes of 300 and 150 scores.
const P2P = (budgets = [200 * 2^k for k in 0:9], b_primary = 12800, bootstrap_resamples = 10_000, bootstrap_level = 0.95, bootstrap_seed = 20261012)
const P2RNG = StableRNG(20261014)
const P2VS = (:chaotic, :regular, :undecided)
p2recs(n) = [(id = o, V = Dict(X => Symbol[P2VS[rand(P2RNG, 1:3)] for _ in 1:10] for X in AN.P2_I7), T = Dict(X => 200 * 2^rand(P2RNG, 0:7) for X in AN.P2_I7)) for o in 1:n]
const P2W, P2S = p2recs(214), p2recs(448)
SUITE["p2_bootstrap H1 214 orbits 1e4 resamples"] = @benchmarkable AN.p2_bootstrap(:H1, ($P2W,), "p2|bench|H1", $P2P)
SUITE["p2_bootstrap H2 448 orbits 1e4 resamples"] = @benchmarkable AN.p2_bootstrap(:H2, ($P2S,), "p2|bench|H2", $P2P)
const P2V8 = Dict(X => Symbol[P2VS[rand(P2RNG, 1:3)] for _ in 1:214] for X in AN.P2_I8)
SUITE["p2_unanimity 28 pairs 214 orbits"] = @benchmarkable AN.p2_unanimity($P2V8, $(AN.P2_I8))
const P2C, P2R = rand(P2RNG, 300), rand(P2RNG, 150)
SUITE["p2_auc 300 against 150 scores"] = @benchmarkable AN.p2_auc($P2C, $P2R; min_set = 20)
