using BenchmarkTools, Chaoslab, StaticArrays

const SUITE = BenchmarkGroup()

wrapped = WrappedStandardMap(5.0)
lifted = StandardMap(5.0)
xw = torus_state(SVector(1.0, 2.0))
xl = SVector(1.0, 2.0)
u0 = SVector(0.6, 0.8)
starts = [torus_state(SVector(0.01 * i, 0.02 * i)) for i in 1:4096]

SUITE["separation_wrapped"] = @benchmarkable separation($wrapped, $xw, $(displaced(wrapped, xw, 1e-8 * u0)))
SUITE["orbit_1000_steps_wrapped"] = @benchmarkable propagate_shadow_orbit($wrapped, $xw, $u0, 1000, 1, 1e-8)
SUITE["orbit_1000_steps_lifted"] = @benchmarkable propagate_shadow_orbit($lifted, $xl, $u0, 1000, 1, 1e-8)
SUITE["propagate_shadow_4096x1000"] = @benchmarkable propagate_shadow($wrapped, $starts, 1000; d0 = 1e-8, tau = 1, seed = 1)
