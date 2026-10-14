# Stage 04 verification of module lyap-householder-qr: the full cross-product of the package-suite testsets (decision D38).
# Run from the repo root:  julia --project=verify -t 2 verify/lyap-householder-qr.jl
# Spec: workspaces/module-build/stages/01-spec/output/lyap-householder-qr-spec.md (criteria 1 to 13 and 15; 14 and 17 are not in this script).
# The package suite (test/lyap-householder-qr/) keeps one representative case per criterion because it is compile-bound (3 minute budget);
# this script includes the same files with HH_FULL_RUN = true, so every shape, variant, element type and K of the original suite runs,
# with identical assertions (helpers hh_sel_algs, hh_pick, hh_sel, hh_ac1_cases in hh_helpers.jl). Output: the Test summary table, one row per testset.
using Test
using Chaoslab, StaticArrays, DoubleFloats, StableRNGs, LinearAlgebra

const HH_FULL_RUN = true
const TESTDIR = joinpath(dirname(@__DIR__), "test", "lyap-householder-qr")

t0 = time()
# verbose summary: one row per testset (AC<n> ...) with pass/fail counts and time; a failing testset throws after the table is printed
@testset verbose = true "lyap-householder-qr full cross-product" begin
    include(joinpath(TESTDIR, "runtests.jl"))
end
println("lyap-householder-qr verify: ", round(time() - t0; digits = 1), " s")
