# Stage 04 verification of module lyap-two-trajectory: the package suite at full size (D38).
# Run from the repo root:  julia --project=verify -t 2 verify/lyap-two-trajectory.jl
# Spec: workspaces/module-build/stages/01-spec/output/lyap-two-trajectory-spec.md (criteria 1 to 12; 13 to 17 need the external oracle and a GPU, not run here).
# The package suite (test/lyap-two-trajectory/runtests.jl) keeps one representative case per criterion and shortened sweeps; with L2_FULL_RUN defined it runs the
# full Double64 grids (AC1, AC2, AC5), the full sweeps (AC6 decades and n = 1e4, AC7 200 starts, AC8 K = 2 pi and every decade, AC9 reported cases, AC11 4096 orbits
# and chunk 1000, AC12 lifted at n = 1e5 and 100 starts). Same assertions and tolerances; the summary table has one line per criterion with passes and time.
using Printf, Test
using Chaoslab, StaticArrays, DoubleFloats, StableRNGs, LinearAlgebra

const L2_FULL_RUN = true
t = @elapsed ts = @testset verbose = true "lyap-two-trajectory verify" begin
    include(joinpath(dirname(@__DIR__), "test", "lyap-two-trajectory", "runtests.jl"))
end
@printf("lyap-two-trajectory verify | %s | %.1f s\n", all(r -> !(r isa Test.Fail || r isa Test.Error), ts.results) ? "PASS" : "FAIL", t)
