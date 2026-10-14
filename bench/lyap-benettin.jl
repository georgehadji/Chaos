using BenchmarkTools, Chaoslab, StaticArrays

const SUITE = BenchmarkGroup()

sys = StandardMap(2.0)
x = SVector(1.0, 2.0)
W = SMatrix{2, 1, Float64}(3, 4)
st = reducer_init(Benettin(), Float64, Val(2), Val(1))
W1 = SMatrix{2, 1, Float64}(1, 0)

SUITE["safe_norm_float64"] = @benchmarkable safe_norm($(SVector(3.0, 4.0)))
SUITE["step_float64"] = @benchmarkable reducer_step($(Benettin()), $st, 1, $x, $W)
SUITE["step_double64"] = @benchmarkable reducer_step($(Benettin()), $(reducer_init(Benettin(), Double64, Val(2), Val(1))), 1,
    $(SVector{2, Double64}(1, 2)), $(SMatrix{2, 1, Double64}(3, 4)))
SUITE["orbit_1000_steps"] = @benchmarkable propagate_orbit($sys, $x, $W1, 1000, 1, ($(Benettin()),))
SUITE["orbit_1000_steps_monitor"] = @benchmarkable propagate_orbit($sys, $x, $W1, 1000, 1,
    ($(Benettin()), $(ResidualMonitor(sys, 4.0, Inf))))
