using BenchmarkTools, Chaoslab, StaticArrays

const SUITE = BenchmarkGroup()

sys = LinearSymplecticMap(SMatrix{2, 2, Float64}(2, 1, 1, 1))
W2 = SMatrix{2, 2, Float64}(2.0, 1.0, 1.0, 1.0)
x2 = SVector(0.0, 0.0)
r = SVDSpectrum()
st = reducer_init(r, Float64, Val(2), Val(2))

SUITE["log_singular_values"] = @benchmarkable log_singular_values($(SVector(12.0, -12.0)), 1.27)
SUITE["step_svd_2x2"] = @benchmarkable reducer_step($r, $st, 1, $x2, $W2)
SUITE["orbit_1000_steps_svd"] = @benchmarkable propagate_orbit($sys, $x2, $(SMatrix{2, 2, Float64}(1, 0, 0, 1)), 1000, 1, ($r,))
