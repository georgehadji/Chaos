using BenchmarkTools, Chaoslab, StaticArrays

const SUITE = BenchmarkGroup()

sys = WrappedStandardMap(6.0)
x0 = torus_state(SVector(1.0, 2.0))
W1 = SMatrix{2, 2, Float64}(1, 0, 0, 1)
W = SMatrix{2, 2, Float64}(3, 4, 0, 5)
range = stretching_range(sys)
r100 = StretchingNumbers(; range = range, bins = 100)
r1000 = StretchingNumbers(; range = range, bins = 1000)
st = reducer_init(r1000, Float64, Val(2), Val(2))
ws = zeros(UInt32, reducer_workspace(r1000, Float64, Val(2), Val(2)))
x8 = zero(SVector{2, Float64})
counts = rand(0:50, 1002)
counts2 = rand(0:50, 1002)
counts2[1] += sum(counts) - sum(counts2)   # equal totals

SUITE["check_frame_float64"] = @benchmarkable check_frame($r1000, $(Ref(W1))[])
SUITE["spectrum_cell_float64"] = @benchmarkable spectrum_cell($(Ref(0.3))[], $(range[1]), $(range[2]), 1000)
SUITE["spectrum_cell_double64"] = @benchmarkable spectrum_cell($(Ref(Double64(0.3)))[], $(Double64(range[1])), $(Double64(range[2])), 1000)
SUITE["spectral_distance_B1000"] = @benchmarkable spectral_distance($counts, $counts2, 0.002)
SUITE["step_float64"] = @benchmarkable reducer_step($r1000, $st, 1, $x8, $(Ref(W))[], $ws)
SUITE["step_float64_reader"] = @benchmarkable reducer_step($(StretchingNumbers(; range = range, replace = false)), $st, 1, $x8, $(Ref(W))[], $ws)
SUITE["step_double64"] = @benchmarkable reducer_step($(StretchingNumbers(; range = range)), $(reducer_init(r1000, Double64, Val(2), Val(2))), 1,
    $(zero(SVector{2, Double64})), $(Ref(SMatrix{2, 2, Double64}(3, 4, 0, 5)))[], $ws)
SUITE["final_with_counts_B1000"] = @benchmarkable reducer_final($r1000, $st, $ws)
# cost against the L1 single column, MEGNO and RLI-free chains on the same orbit (R6)
SUITE["orbit_1000_steps_L1"] = @benchmarkable propagate_orbit($sys, $x0, $(SMatrix{2, 1, Float64}(1, 0)), 1000, 1, ($(Benettin()),))
SUITE["orbit_1000_steps_megno"] = @benchmarkable propagate_orbit($sys, $x0, $W1, 1000, 1, ($(MEGNO()),))
SUITE["orbit_1000_steps_B100"] = @benchmarkable propagate_orbit($sys, $x0, $W1, 1000, 1, ($r100,); workspace = $(zeros(UInt32, 204)))
SUITE["orbit_1000_steps_B1000"] = @benchmarkable propagate_orbit($sys, $x0, $W1, 1000, 1, ($r1000,); workspace = $(zeros(UInt32, 2004)))
SUITE["orbit_1000_steps_chain"] = @benchmarkable propagate_orbit($sys, $x0, $W1, 1000, 1, ($(StretchingNumbers(; range = range, replace = false)), $(MEGNO())))
SUITE["propagate_1024_orbits_1000_steps_B1000"] = @benchmarkable propagate($sys, $([torus_state(SVector(0.1 * i, 2.0)) for i in 1:1024]), 1000;
    k = 2, tau = 1, seed = 1, reducers = ($r1000,), backend = :cpu)
SUITE["propagate_1024_orbits_1000_steps_megno"] = @benchmarkable propagate($sys, $([torus_state(SVector(0.1 * i, 2.0)) for i in 1:1024]), 1000;
    k = 2, tau = 1, seed = 1, reducers = ($(MEGNO()),), backend = :cpu)
# the RLI pair on the same map (one pair is one orbit of the product map)
pair = PairMap(StandardMap(6.0))
zp = pair_start(StandardMap(6.0), SVector(1.0, 2.0), SVector(0.6, 0.8), 1e-9).z
Wp = SMatrix{4, 2, Float64}(0.6, 0.8, 0, 0, 0, 0, 0.6, 0.8)
SUITE["orbit_1000_steps_rli_pair"] = @benchmarkable propagate_orbit($pair, $zp, $Wp, 1000, 1, ($(RLI()),))
