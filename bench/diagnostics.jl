using BenchmarkTools, Chaoslab, StaticArrays

const SUITE = BenchmarkGroup()

sys = StandardMap(2.0)
mon = ResidualMonitor(sys, 4.0, 1.0e-6)
x = SVector(1.0, 2.0)
W2 = SMatrix{2, 2, Float64}(1, 0, 0, 1)
st = reducer_init(mon, Float64, Val(2), Val(2))
lam = SVector(0.9, -0.9)

SUITE["update_float64"] = @benchmarkable reducer_update($mon, $st, 1, $x, $W2)
SUITE["update_double64"] = @benchmarkable reducer_update($(ResidualMonitor(StandardMap(Double64(2)), 4.0, 1.0e-6)),
    $(reducer_init(mon, Double64, Val(2), Val(2))), 1, $(SVector{2, Double64}(1, 2)), $(SMatrix{2, 2, Double64}(1, 0, 0, 1)))
SUITE["exponents_float64"] = @benchmarkable monitor_exponents($st, $lam, 6.85, 1; pairing_tol = 16.0, tau = 1, M_max = 2.0,
    step = 1, orbit = 1)
SUITE["orbit_1000_steps"] = @benchmarkable propagate_orbit($sys, $x, $(SMatrix{2, 1, Float64}(1, 0)), 1000, 1, ($mon,))
