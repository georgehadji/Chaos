using BenchmarkTools, Chaoslab, StaticArrays

const SUITE = BenchmarkGroup()

sys = StandardMap(6.0)
x = SVector(1.0, 2.0)
W = SMatrix{2, 2, Float64}(3, 4, 0, 5)
r = MEGNO()
st = reducer_init(r, Float64, Val(2), Val(2))
W1 = SMatrix{2, 2, Float64}(1, 0, 0, 1)

SUITE["check_frame_float64"] = @benchmarkable check_frame($r, $(Ref(W1))[])
SUITE["step_float64"] = @benchmarkable reducer_step($r, $st, 1, $x, $(Ref(W))[])
SUITE["step_float64_reader"] = @benchmarkable reducer_step($(MEGNO(replace = false)), $st, 1, $x, $(Ref(W))[])
SUITE["step_double64"] = @benchmarkable reducer_step($r, $(reducer_init(r, Double64, Val(2), Val(2))), 1,
    $(SVector{2, Double64}(1, 2)), $(Ref(SMatrix{2, 2, Double64}(3, 4, 0, 5)))[])
SUITE["orbit_1000_steps"] = @benchmarkable propagate_orbit($sys, $x, $W1, 1000, 1, ($r,))
SUITE["orbit_1000_steps_sali"] = @benchmarkable propagate_orbit($sys, $x, $W1, 1000, 1, ($(SALI()),))
SUITE["orbit_1000_steps_chain"] = @benchmarkable propagate_orbit($sys, $x, $W1, 1000, 1, ($(MEGNO(replace = false)), $(SALI())))
SUITE["propagate_1024_orbits_1000_steps"] = @benchmarkable propagate($sys, $([SVector(0.1 * i, 2.0) for i in 1:1024]), 1000;
    k = 2, tau = 1, seed = 1, reducers = ($r,), backend = :cpu)
