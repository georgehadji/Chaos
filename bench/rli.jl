using BenchmarkTools, Chaoslab, StaticArrays

const SUITE = BenchmarkGroup()

base = StandardMap(6.0)
sys = PairMap(base)
u = SVector(0.6, 0.8)
start = pair_start(base, SVector(1.0, 2.0), u, 1e-9)
z = start.z
W0 = SMatrix{4, 2, Float64}(0.6, 0.8, 0, 0, 0, 0, 0.6, 0.8)
W = SMatrix{4, 2, Float64}(3, 4, 0, 0, 0, 0, 5, 12)
r = RLI()
st = reducer_init(r, Float64, Val(4), Val(2))
x8 = zero(SVector{4, Float64})

SUITE["check_frame_float64"] = @benchmarkable check_frame($r, $(Ref(W0))[])
SUITE["pair_start_float64"] = @benchmarkable pair_start($base, $(Ref(SVector(1.0, 2.0)))[], $u, 1e-9)
SUITE["pair_tangent_action"] = @benchmarkable tangent($sys, $z) * $(Ref(W0))[]
SUITE["step_float64"] = @benchmarkable reducer_step($r, $st, 1, $x8, $(Ref(W))[])
SUITE["step_float64_reader"] = @benchmarkable reducer_step($(RLI(replace = false)), $st, 1, $x8, $(Ref(W))[])
SUITE["step_double64"] = @benchmarkable reducer_step($r, $(reducer_init(r, Double64, Val(4), Val(2))), 1,
    $(zero(SVector{4, Double64})), $(Ref(SMatrix{4, 2, Double64}(3, 4, 0, 0, 0, 0, 5, 12)))[])
SUITE["pair_orbit_1000_steps"] = @benchmarkable propagate_orbit($sys, $z, $W0, 1000, 1, ($r,))
SUITE["pair_orbit_1000_steps_chain"] = @benchmarkable propagate_orbit($sys, $z, $W0, 1000, 1, ($(RLI(replace = false)), $(MEGNO())))
SUITE["single_orbit_1000_steps"] = @benchmarkable propagate_orbit($base, $(SVector(1.0, 2.0)), $(SMatrix{2, 2, Float64}(1, 0, 0, 1)), 1000, 1, ($(MEGNO()),))
SUITE["propagate_1024_pairs_1000_steps"] = @benchmarkable propagate($sys, $([pair_start(base, SVector(0.1 * i, 2.0), u, 1e-9).z for i in 1:1024]), 1000;
    k = 2, tau = 1, seed = 1, reducers = ($r,), backend = :cpu, W0 = $(fill(W0, 1024)))
