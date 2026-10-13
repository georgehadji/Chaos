using BenchmarkTools, Chaoslab, StaticArrays

const SUITE = BenchmarkGroup()
const SM = StandardMap(5.0)
const Z0 = SVector(1.25, -0.5)
const WSM = WrappedStandardMap(5.0)
const S0 = SVector(1.25, 5.75, 0.0, 0.0)

SUITE["forward"] = @benchmarkable forward($SM, $Z0)
SUITE["tangent"] = @benchmarkable tangent($SM, $Z0)
SUITE["forward 1000 steps"] = @benchmarkable begin
    z = $Z0
    for _ in 1:1000
        z = forward($SM, z)
    end
    z
end
SUITE["linear symplectic forward (N=3)"] = @benchmarkable forward(m, z) setup = (
    m = LinearSymplecticMap(block_diagonal(rotation_block(Float64, 0.3), hyperbolic_block(Float64, 0.2), hyperbolic_block(Float64, 0.4)));
    z = SVector(0.1, 0.2, 0.3, 0.4, 0.5, 0.6)
)
SUITE["wrapped forward"] = @benchmarkable forward($WSM, $S0)
SUITE["wrapped forward 1000 steps"] = @benchmarkable begin
    s = $S0
    for _ in 1:1000
        s = forward($WSM, s)
    end
    s
end
