using BenchmarkTools, Chaoslab

# Hot Foundation functions used by the precision paths: compensated accumulation, exact conversion.
const SUITE = BenchmarkGroup()
SUITE["add Float64"] = @benchmarkable add(acc, 0.1) setup = (acc = CompensatedSum{Float64}())
SUITE["add Double64"] = @benchmarkable add(acc, x) setup = (acc = CompensatedSum{Double64}(); x = Double64(0.1))
SUITE["exact(Double64, pi)"] = @benchmarkable exact(Double64, pi)
SUITE["resolve_backend"] = @benchmarkable resolve_backend(:auto, Float64; probe = () -> (loadable = false, functional = false, device = "", compiles = false, gpu_faster = false))
