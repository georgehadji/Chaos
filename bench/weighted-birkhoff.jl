using BenchmarkTools, Chaoslab, StaticArrays

# Budget (stage 03, orientation only, the spec gates no time): the reader step at (K, H) = (4, 1) under 0.5 microsecond in Float64 with 0 allocations,
# one orbit of 1000 steps at (4, 1) within 3 times the L1 stream alone, at (4, 11) within 8 times; the index and the pure function in the microseconds.
const SUITE = BenchmarkGroup()

sys = WrappedStandardMap(0.5)
x0 = torus_state(SVector(1.0, 2.2))
W1 = SMatrix{2, 1, Float64}(1, 0)
xw = SVector(0.3, 2.2, 1.0, 0.0)
obs1 = (CosObservable(1),)
obs4 = (CosObservable(1), SinObservable(1), CosObservable(2), LiftedObservable(2, 4))
grid11 = Tuple(2^k for k in 4:14)   # eleven horizons, N = 16 .. 16384 (2 max(window) handoffs are not run here: the steps are the cost)

# the pure function and the index
series(n, K, ::Type{T} = Float64) where {T} = T[sin(T(0.7) * j + k) for j in 0:(n - 1), k in 1:K]
s1k, s10k = series(2000, 4), series(20_000, 4)
SUITE["weighted_birkhoff_K4_N1000"] = @benchmarkable weighted_birkhoff($s1k, 1000)
SUITE["weighted_birkhoff_K4_N10000"] = @benchmarkable weighted_birkhoff($s10k, 10_000)
SUITE["weighted_birkhoff_K4_N1000_double64"] = @benchmarkable weighted_birkhoff($(series(2000, 4, Double64)), 1000)
SUITE["birkhoff_index_float64"] = @benchmarkable birkhoff_index($(0.3), $(0.3 + 1e-9), $(1.0))

# the reducer per handoff, (K, H) = (1, 1), (4, 1), (4, 11)
for (name, obs, win) in (("K1_H1", obs1, 100), ("K4_H1", obs4, 100), ("K4_H11", obs4, grid11))
    r = WeightedBirkhoff(; window = win, observables = obs)
    st = reducer_step(r, reducer_init(r, Float64, Val(4), Val(1)), 1, xw, W1)[1]
    SUITE["step_$(name)_float64"] = @benchmarkable reducer_step($r, s, 2, x, $W1) setup = (s = $st; x = $xw) evals = 1000
    SUITE["final_$(name)_float64"] = @benchmarkable reducer_final($r, s) setup = (s = $st) evals = 1000
end
r4d = WeightedBirkhoff(; window = 100, observables = obs4)
SUITE["step_K4_H1_double64"] = @benchmarkable reducer_step($r4d, s, 2, x, $(SMatrix{2, 1, Double64}(1, 0))) setup =
    (s = $(reducer_step(r4d, reducer_init(r4d, Double64, Val(4), Val(1)), 1, SVector{4, Double64}(0.3, 2.2, 1.0, 0.0), SMatrix{2, 1, Double64}(1, 0))[1]);
     x = SVector{4, Double64}(0.31, 2.2, 1.0, 0.0)) evals = 1000

# cost against the L1 single column and E2 on the same orbit (spec R5); window 400 gives 800 of the 1000 handoffs
SUITE["orbit_1000_steps_L1"] = @benchmarkable propagate_orbit($sys, $x0, $W1, 1000, 1, ($(Benettin()),))
SUITE["orbit_1000_steps_e2_lifted"] = @benchmarkable propagate_orbit($sys, $x0, $W1, 1000, 1, ($(WeightedPhaseAdvance(; window = 400, observable = LiftedAngleObservable(1, 3))),))
for (name, obs, win) in (("K1_H1", obs1, 400), ("K4_H1", obs4, 400), ("K4_H11", obs4, grid11))
    SUITE["orbit_1000_steps_wb_$name"] = @benchmarkable propagate_orbit($sys, $x0, $W1, 1000, 1, ($(WeightedBirkhoff(; window = win, observables = obs)),))
end
SUITE["orbit_1000_steps_wb_megno_chain"] = @benchmarkable propagate_orbit($sys, $x0, $(SMatrix{2, 2, Float64}(1, 0, 0, 1)), 1000, 1,
    ($(WeightedBirkhoff(; window = 400, observables = obs4)), $(MEGNO())))
SUITE["propagate_1024_orbits_1000_steps_wb_K4_H1"] = @benchmarkable propagate($sys, $([torus_state(SVector(0.1 * i, 2.2)) for i in 1:1024]), 1000;
    k = 1, tau = 1, seed = 1, reducers = ($(WeightedBirkhoff(; window = 400, observables = obs4)),), backend = :cpu)
