using BenchmarkTools, Chaoslab, StaticArrays

const SUITE = BenchmarkGroup()

sys = StandardMap(2.0)
W4 = SMatrix{4, 4, Float64}(4.0, 0.3, 0.2, 0.1, 1.0, 3.0, 0.4, 0.2, 0.5, 1.0, 5.0, 0.3, 0.2, 0.1, 1.0, 2.0)
W2 = SMatrix{2, 2, Float64}(2.0, 1.0, 1.0, 1.0)
x2 = SVector(1.0, 2.0)

for alg in (HouseholderQR(), AugmentedQR(), GivensQR(), LapackQR())
    name = string(nameof(typeof(alg)))
    SUITE["orthonormalize_4x4_" * name] = @benchmarkable orthonormalize($alg, $W4)
    r = QRSpectrum(alg)
    st = reducer_init(r, Float64, Val(2), Val(2))
    SUITE["step_2x2_" * name] = @benchmarkable reducer_step($r, $st, 1, $x2, $W2)
end
r = QRSpectrum(HouseholderQR())
SUITE["orbit_1000_steps_householder"] = @benchmarkable propagate_orbit($sys, $x2, $(SMatrix{2, 2, Float64}(1, 0, 0, 1)), 1000, 1, ($r,))
