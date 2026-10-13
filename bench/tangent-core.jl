using BenchmarkTools, Chaoslab, StaticArrays

# Hot functions of the tangent stream: the per-orbit body, initial deviations, and the chunked kernel launch.
const SUITE = BenchmarkGroup()
const SM = StandardMap(0.5)
const X0 = SVector(1.0, 0.5)
const W0 = SMatrix{2, 2, Float64}(1, 0, 0, 1)

SUITE["propagate_orbit 100 steps"] = @benchmarkable propagate_orbit($SM, $X0, $W0, 100, 1, ())
SUITE["initial_deviations d=4 k=4"] = @benchmarkable initial_deviations(1, 7, Val(4), Val(4), Float64)
SUITE["propagate 4096 orbits x 100 steps"] = @benchmarkable propagate($SM, xs, 100; k = 2, seed = 1) setup = (
    xs = [SVector(0.01 * j, 0.5) for j in 1:4096]
)
