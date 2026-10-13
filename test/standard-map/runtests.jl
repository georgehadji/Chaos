# Suite for module standard-map (Systems: Y1 standard map, Y5 linear symplectic maps).
# Spec: workspaces/module-build/stages/01-spec/output/standard-map-spec.md; plan: stages/02-tests/output/standard-map-test-plan.md.
# Criterion 19 (ChaosTools cross-check) and the GPU and heavy parts of 14 and 18 are stage 04 and have no test here.
using Test, Chaoslab, StaticArrays, LinearAlgebra, StableRNGs, DoubleFloats, ForwardDiff

# The outer testset keeps a failing criterion from aborting the criteria after it when this file is run on its own.
@testset "standard-map" verbose = true begin
    include("sm_helpers.jl")
    include("sm_ac01_06.jl")
    include("sm_ac07_13.jl")
    include("sm_ac14.jl")
    include("sm_y5.jl")
    include("sm_ac18.jl")
    include("sm_wrapped.jl")
end
