using BenchmarkTools, Chaoslab, StaticArrays

const SUITE = BenchmarkGroup()

sys = StandardMap(6.0)
x = SVector(1.0, 2.0)
W = SMatrix{2, 2, Float64}(3, 4, 0, 5)
r = SALI()
st = reducer_init(r, Float64, Val(2), Val(2))
W1 = SMatrix{2, 2, Float64}(1, 0, 0, 1)

SUITE["sali_pair_float64"] = @benchmarkable sali_pair($(Ref(SVector(3.0, 4.0)))[], $(Ref(SVector(0.0, 5.0)))[])
SUITE["sali_from_log_float64"] = @benchmarkable sali_from_log($(Ref(-12.0))[])
SUITE["check_frame_float64"] = @benchmarkable check_frame($r, $(Ref(W))[])
SUITE["step_float64"] = @benchmarkable reducer_step($r, $st, 1, $x, $(Ref(W))[])
SUITE["step_double64"] = @benchmarkable reducer_step($r, $(reducer_init(r, Double64, Val(2), Val(2))), 1,
    $(SVector{2, Double64}(1, 2)), $(Ref(SMatrix{2, 2, Double64}(3, 4, 0, 5)))[])
SUITE["orbit_1000_steps"] = @benchmarkable propagate_orbit($sys, $x, $W1, 1000, 1, ($r,))
SUITE["orbit_1000_steps_benettin"] = @benchmarkable propagate_orbit($sys, $x, $(SMatrix{2, 1, Float64}(1, 0)), 1000, 1, ($(Benettin()),))
SUITE["propagate_1024_orbits_1000_steps"] = @benchmarkable propagate($sys, $([SVector(0.1 * i, 2.0) for i in 1:1024]), 1000;
    k = 2, tau = 1, seed = 1, reducers = ($r,), backend = :cpu)
