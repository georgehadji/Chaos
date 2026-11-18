using BenchmarkTools, Chaoslab, StaticArrays

const SUITE = BenchmarkGroup()

# Continuous QR (L8): the right-hand side of the penalty form, one Q-step of each order, a full run of 100 time units (2000 orbit steps of the
# order-8 table, h = 1/20) on Henon-Heiles orbit C, order 4 against order 2, and a 256 orbit ensemble. The discrete RGS spectrum (R1) is the reference.
const HH = HenonHeilesFlow()
const m8 = KahanLi(Float64, 8)
const im = IntegratorMap(m8, HH, 1 / 20)
const zC = first(energy_start(HH, 1 / 8, -0.015, 0.0))          # orbit C (p1 rebuilt)
const I4 = one(SMatrix{4, 4, Float64})
const A4 = vector_field_jacobian(HH, zC)
const cq4 = ContinuousQR(im, HH; order = 4)
const cq2 = ContinuousQR(im, HH; order = 2)

SUITE["rhs_form84"] = @benchmarkable cqr_rhs($(Ref(I4))[], $(Ref(A4))[], 1.0; form = :f84)
SUITE["growth"] = @benchmarkable cqr_growth($(Ref(I4))[], $(Ref(A4))[])
SUITE["qstep_order4_2_orbit_steps"] = @benchmarkable propagate_cqr_orbit($cq4, $zC, $I4, 2)
SUITE["qstep_order2_1_orbit_step"] = @benchmarkable propagate_cqr_orbit($cq2, $zC, $I4, 1)
SUITE["orbit_100_time_units_order4"] = @benchmarkable propagate_cqr_orbit($cq4, $zC, $I4, 2000)
SUITE["orbit_100_time_units_order2"] = @benchmarkable propagate_cqr_orbit($cq2, $zC, $I4, 2000)
SUITE["discrete_rgs_100_time_units"] = @benchmarkable propagate_orbit($im, $zC, $I4, 2000, 20, ($(QRSpectrum(RGS())),))
SUITE["propagate_256_orbits_100_time_units"] = @benchmarkable propagate_cqr($cq4, $([SVector(zC[1], zC[2] + 1.0e-4 * i, zC[3], zC[4]) for i in 1:256]), 2000;
    k = 4, seed = 1, backend = :cpu)
