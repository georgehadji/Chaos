using BenchmarkTools, Chaoslab, StaticArrays

const SUITE = BenchmarkGroup()

sys = LinearSymplecticMap(SMatrix{2, 2, Float64}(2, 1, 1, 1))
W2 = SMatrix{2, 2, Float64}(2.0, 1.0, 1.0, 1.0)
W4 = SMatrix{4, 4, Float64}(4.0, 0.3, 0.2, 0.1, 1.0, 3.0, 0.4, 0.2, 0.5, 1.0, 5.0, 0.3, 0.2, 0.1, 1.0, 2.0)
x2 = SVector(0.0, 0.0)

for alg in (CGS(), MGS(), RGS())
    name = string(nameof(typeof(alg)))
    SUITE["orthonormalize_4x4_" * name] = @benchmarkable orthonormalize($alg, $W4)
end
r = QRSpectrum(RGS())
st = reducer_init(r, Float64, Val(2), Val(2))
SUITE["step_rgs_2x2"] = @benchmarkable reducer_step($r, $st, 1, $x2, $W2)
SUITE["orbit_1000_steps_rgs"] = @benchmarkable propagate_orbit($sys, $x2, $(SMatrix{2, 2, Float64}(1, 0, 0, 1)), 1000, 1, ($r,))
