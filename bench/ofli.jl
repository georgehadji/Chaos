using BenchmarkTools, Chaoslab, StaticArrays

const SUITE = BenchmarkGroup()

# OFLI (I6): one reducer update (Binary64 and Double64, replacing and reader), the s = 0 terms, the frame check, and a full run of 1000 time units
# (2 10^4 steps of the order-8 table, tau = 20) on Henon-Heiles orbit C with k = 4, alone and in a chain with SALI (k = 2), and a 1024 orbit ensemble.
const HH = HenonHeilesFlow()
const m8 = KahanLi(Float64, 8)
const im = IntegratorMap(m8, HH, 1 / 20)
const zC = first(energy_start(HH, 1 / 8, -0.015, 0.0))          # orbit C (p1 rebuilt)
const I4 = one(SMatrix{4, 4, Float64})
const W4 = SMatrix{4, 4, Float64}(1.0, 0.2, 0.1, 0.0, 0.3, 2.0, 0.0, 0.1, 0.0, 0.5, 3.0, 0.2, 0.1, 0.0, 0.4, 5.0)
const W2 = SMatrix{4, 2, Float64}(I4[:, 1:2])
const r = OFLI(HH)
const st = reducer_init(r, Float64, Val(4), Val(4))
const rD = OFLI(HH)
const zD = SVector{4, Double64}(zC)

SUITE["step_float64"] = @benchmarkable reducer_step($r, $st, 20, $(Ref(zC))[], $(Ref(W4))[])
SUITE["step_float64_reader"] = @benchmarkable reducer_step($(OFLI(HH; replace = false)), $st, 20, $(Ref(zC))[], $(Ref(W4))[])
SUITE["step_double64"] = @benchmarkable reducer_step($rD, $(reducer_init(rD, Double64, Val(4), Val(4))), 20, $(Ref(zD))[],
    $(Ref(SMatrix{4, 4, Double64}(W4)))[])
SUITE["ofli_initial"] = @benchmarkable ofli_initial($HH, $(Ref(zC))[], $(Ref(I4))[])
SUITE["check_frame"] = @benchmarkable check_frame($r, $(Ref(I4))[], $(Ref(zC))[])
SUITE["orbit_1000_time_units_k4"] = @benchmarkable propagate_orbit($im, $zC, $I4, 20_000, 20, ($r,))
SUITE["orbit_1000_time_units_chain_sali"] = @benchmarkable propagate_orbit($im, $zC, $W2, 20_000, 20, ($(OFLI(HH; replace = false)), $(SALI())))
SUITE["propagate_1024_orbits_200_time_units"] = @benchmarkable propagate($im, $([SVector(zC[1], zC[2] + 1.0e-4 * i, zC[3], zC[4]) for i in 1:1024]), 4000;
    k = 4, tau = 20, seed = 1, reducers = ($r,), backend = :cpu)
