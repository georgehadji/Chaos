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

# Study p3-abstention (Amendment D93): the hot paths of `calibrate` at the sizes of the real study (about 3000 orbits, 2088 elements, 200 optimism resamples, 20000 runs per path point).
# The parameter table is the one of the frozen preregistration for the budgets, the costs and the sizes; the thresholds are those of a fixture (the speed does not depend on them).
const P3Q = (budgets = [200 * 2^k for k in 0:9], b_primary = 12800, N_c = 102400, indicators_i7 = ["LE", "SALI", "FLI", "MEGNO", "E1", "E2", "WBA"], voter_bits = ["SALI", "FLI", "MEGNO", "E1", "E2", "WBA"],
    w_range = [0, 1, 2], v_range = [0, 1], min_voters = 3, min_set = 20, cost_l1 = 207.6, cost_tg = 2658.0, cost_pa1 = 1810.5, cost_bk1 = 3000.6, cost_fq1 = 16850.4,
    spec_of = Dict("LE" => "l1", "SALI" => "tg", "FLI" => "tg", "MEGNO" => "tg", "E1" => "fq1", "E2" => "pa1", "WBA" => "bk1"), e_cal = 0.2, a_cal = 0.3, s_cal = 0.15, c_min = 0.5, delta_tol = 0.05,
    alpha_hb = 0.05, rho_hb = 0.5, n_t_min = 100, n_w_min = 300, p_hb = 0.8, b_opt = 1, rule_rate_floor = 0.01, rule_rate_factor = 2.0, test_seed = 20261030, r_max = 5.0, rho_path_start = 1.0,
    rho_path_second = 1.25, y_t = 57, y_w = 247, ceil_step = 5, shared_default = 0.5)
const P3RNG = StableRNG(20261015)
const P3VS = (:chaotic, :regular, :undecided)
# 3000 orbits with random verdicts at the ten budgets, a regular share, and the sets of the real study in size (S2 448, W2 214 of which T_s 30, gate 80 in the regular set)
function p3rows(n)
    rows = NamedTuple[]
    for id in 1:n
        rho = rand(P3RNG) < 0.6 ? :C : :R
        V = Dict(X => Symbol[rand(P3RNG) < 0.2 ? :undecided : (rho === :C) == (rand(P3RNG) < 0.9) ? :chaotic : :regular for _ in 1:10] for X in AN.P2_I7)
        push!(rows, (id = id, source = id <= 883 ? :panel : :grid, kindex = 1 + id % 9, in_cal = true, rho = rho, S2 = rho === :C && id % 7 == 0, Ts = rho === :C && id % 60 == 0, W2 = rho === :C && id % 9 == 0,
                     gate = rho === :R && id % 40 == 0, k0 = false, V = V))
    end
    return rows
end
const P3ROWS = p3rows(3000)
SUITE["p3_select 2088 elements 3000 orbits (C3 to C6)"] = @benchmarkable AN.p3_select($P3ROWS, $P3Q)
SUITE["p3_optimism one resample (C3 to C6 on a resample of 3000 orbits)"] = @benchmarkable AN.p3_optimism($P3ROWS, $P3Q; threads = 1)
const P3MODEL = (rule = Dict(:T => (p = 0.1, d = 0.9), :RW => (p = 0.1, d = 0.9)), comp = Dict(X => (T = (p = 0.4, s = 0.5), RW = (p = 0.4, s = 0.5)) for X in AN.P2_I7),
                 ts = [:SALI, :FLI, :MEGNO, :E1, :E2, :WBA], w2 = AN.P2_I7)
SUITE["p3_simulate_hb one path point (285 + 950 orbits) 100 runs"] = @benchmarkable AN.p3_simulate_hb($P3MODEL, 285, 950, $P3Q; point = 16, runs = 100, threads = 1)
const P3QHC = merge(P3Q, (auc_p = 0.55, auc_shrink = 0.5, nonfail_per_fail = 20, nf_start = 20, nf_second = 25, delta_c = 0.1, z_one_sided = 1.645, d_expected_total = 1048, d_worst_total = 1309, p_hc = 0.8, alpha_ci = 0.05))
SUITE["p3_simulate_hc n_f = 65 50 runs"] = @benchmarkable AN.p3_simulate_hc($P3QHC, 0.75, 0.55; nf = 65, runs = 50, threads = 1)
