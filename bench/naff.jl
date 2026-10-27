using BenchmarkTools, Chaoslab, StaticArrays

const SUITE = BenchmarkGroup()

sys = WrappedStandardMap(0.5)
x0 = torus_state(SVector(1.0, 2.2))
W1 = SMatrix{2, 1, Float64}(1, 0)
W2 = SMatrix{2, 2, Float64}(1, 0, 0, 1)

# two lines, the second a tenth of the first and 0.0113 away (the F4 family of the spec): the series of the window pair
series(n, ::Type{T} = Float64) where {T} = Complex{T}[cis(2T(pi) * T(0.1234567891) * j) + T(0.1) * cis(2T(pi) * T(0.1347567891) * j) for j in 0:(n - 1)]
z1k, z10k, z100k = series(2 * 1000), series(2 * 10_000), series(2 * 100_000)

# series path (E1): one window, the pair, the weights, at N = 10^3, 10^4, 10^5 samples per window
SUITE["window_frequency_N1000"] = @benchmarkable window_frequency($(view(z1k, 1:1000)))
SUITE["window_frequency_N10000"] = @benchmarkable window_frequency($(z10k[1:10_000]))
SUITE["window_frequency_N100000"] = @benchmarkable window_frequency($(z100k[1:100_000])) samples = 5 evals = 1
SUITE["frequency_map_N1000"] = @benchmarkable frequency_map($z1k, 1000)
SUITE["frequency_map_N10000"] = @benchmarkable frequency_map($z10k, 10_000)
SUITE["frequency_map_N1000_double64"] = @benchmarkable frequency_map($(series(2 * 1000, Double64)), 1000)
SUITE["window_weights_N10000"] = @benchmarkable window_weights(Float64, 10_000)
SUITE["birkhoff_weights_N10000"] = @benchmarkable birkhoff_weights(Float64, 10_000)

# E2: the pure function on the increments of the pair, and the reducer per handoff
inc = phase_increments(z10k)
SUITE["phase_increments_N10000"] = @benchmarkable phase_increments($z10k)
SUITE["weighted_phase_frequency_N9999"] = @benchmarkable weighted_phase_frequency($inc, 9_999)
r = WeightedPhaseAdvance(; window = 100, observable = PlanarObservable(1, 2))
rl = WeightedPhaseAdvance(; window = 100, observable = LiftedAngleObservable(1, 3))
st = reducer_step(r, reducer_init(r, Float64, Val(2), Val(1)), 1, SVector(0.3, 0.4), W1)[1]
xw = SVector(0.3, 2.2, 1.0, 0.0)
stl = reducer_step(rl, reducer_init(rl, Float64, Val(4), Val(1)), 1, xw, W1)[1]
SUITE["step_planar_float64"] = @benchmarkable reducer_step($r, s, 2, x, $W1) setup = (s = $st; x = SVector(0.31, 0.41)) evals = 1000
SUITE["step_lifted_float64"] = @benchmarkable reducer_step($rl, s, 2, x, $W1) setup = (s = $stl; x = $xw) evals = 1000
SUITE["step_planar_double64"] = @benchmarkable reducer_step($r, s, 2, x, $(SMatrix{2, 1, Double64}(1, 0))) setup = (s = $(reducer_step(r, reducer_init(r, Double64, Val(2), Val(1)), 1, SVector(Double64(0.3), Double64(0.4)), SMatrix{2, 1, Double64}(1, 0))[1]); x = SVector(Double64(0.31), Double64(0.41))) evals = 1000
SUITE["final_float64"] = @benchmarkable reducer_final($r, s) setup = (s = $st) evals = 1000

# cost against the L1 single column, MEGNO and the stretching numbers on the same orbit (spec R5); window 400 gives 801 of the 1000 handoffs
rw = WeightedPhaseAdvance(; window = 400, observable = LiftedAngleObservable(1, 3))
srange = stretching_range(sys)
SUITE["orbit_1000_steps_L1"] = @benchmarkable propagate_orbit($sys, $x0, $W1, 1000, 1, ($(Benettin()),))
SUITE["orbit_1000_steps_e2_lifted"] = @benchmarkable propagate_orbit($sys, $x0, $W1, 1000, 1, ($rw,))
SUITE["orbit_1000_steps_megno"] = @benchmarkable propagate_orbit($sys, $x0, $W2, 1000, 1, ($(MEGNO()),))
SUITE["orbit_1000_steps_stretching"] = @benchmarkable propagate_orbit($sys, $x0, $W2, 1000, 1, ($(StretchingNumbers(; range = srange)),); workspace = $(zeros(UInt32, 2004)))
SUITE["orbit_1000_steps_e2_megno_chain"] = @benchmarkable propagate_orbit($sys, $x0, $W2, 1000, 1, ($rw, $(MEGNO())))
SUITE["propagate_1024_orbits_1000_steps_e2"] = @benchmarkable propagate($sys, $([torus_state(SVector(0.1 * i, 2.2)) for i in 1:1024]), 1000;
    k = 1, tau = 1, seed = 1, reducers = ($rw,), backend = :cpu)
