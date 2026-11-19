using BenchmarkTools, Chaoslab, StaticArrays

const SUITE = BenchmarkGroup()

# Covariant Lyapunov vectors, Ginelli's method (L12): the same orbit as the continuous-QR bench (Henon-Heiles orbit C, order-8 table, h = 1/20, tau = 20),
# n_f = n_b = 100 and n_s = 100 handoffs of one time unit, store-all against checkpointed, with and without the phase-space vectors, the delivered discrete
# RGS spectrum on the same 300 handoffs as the reference of the forward cost (R1), and a 256 orbit ensemble.
const HH = HenonHeilesFlow()
const m8 = KahanLi(Float64, 8)
const im = IntegratorMap(m8, HH, 1 / 20)
const zC = first(energy_start(HH, 1 / 8, -0.015, 0.0))          # orbit C (p1 rebuilt)
const I4 = one(SMatrix{4, 4, Float64})
const clv = ClvGinelli(im; tau = 20, h = 1 / 20)
const clvc = ClvGinelli(im; tau = 20, h = 1 / 20, storage = :checkpoint, block = 25)
const clvr = ClvGinelli(im; tau = 20, h = 1 / 20, orthonormalizer = RGS())
const std = StandardMap(5.0)
const cls = ClvGinelli(std)

SUITE["orbit_store_outputs_none"] = @benchmarkable propagate_clv_orbit($clv, $zC, 100, 100, 100; k = 4, seed = 1)
SUITE["orbit_store_V"] = @benchmarkable propagate_clv_orbit($clv, $zC, 100, 100, 100; k = 4, seed = 1, outputs = (:V,))
SUITE["orbit_checkpoint_V"] = @benchmarkable propagate_clv_orbit($clvc, $zC, 100, 100, 100; k = 4, seed = 1, outputs = (:V,))
SUITE["orbit_store_RGS"] = @benchmarkable propagate_clv_orbit($clvr, $zC, 100, 100, 100; k = 4, seed = 1)
SUITE["orbit_store_angles"] = @benchmarkable propagate_clv_orbit($clv, $zC, 100, 100, 100; k = 4, seed = 1, angle_sets = (1:2, 3:4))
SUITE["orbit_store_k1"] = @benchmarkable propagate_clv_orbit($clv, $zC, 100, 100, 100; k = 1, seed = 1)
SUITE["discrete_rgs_300_handoffs"] = @benchmarkable propagate_orbit($im, $zC, $I4, 6000, 20, ($(QRSpectrum(RGS())),))
SUITE["standard_map_1000_handoffs"] = @benchmarkable propagate_clv_orbit($cls, $(SVector(1.3, 0.7)), 100, 500, 400; k = 2, seed = 1)
SUITE["propagate_256_orbits_standard_map"] = @benchmarkable propagate_clv($cls, $([SVector(1.3 + 1.0e-3 * i, 0.7) for i in 1:256]), 100, 100, 100; k = 2, seed = 1, ntasks = 1)
