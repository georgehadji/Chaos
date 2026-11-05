# Suite for module coupled-maps (Systems: Y3 Froeschle 4D map, Y4 chain of N coupled standard maps, wrapped variants).
# Spec: workspaces/module-build/stages/01-spec/output/coupled-maps-spec.md; plan: stages/02-tests/output/coupled-maps-test-plan.md.
# Criteria 20 (published behaviour of orbits C and D) and 21 (ChaosTools oracle), the timing part of 17, the GPU parts of 18 and the
# 1 vs 8 thread comparison are stage 04 and have no test here.
using Test, Chaoslab, StaticArrays, LinearAlgebra, StableRNGs, DoubleFloats, ForwardDiff

# The outer testset keeps a failing criterion from aborting the criteria after it when this file is run on its own.
@testset "coupled-maps" verbose = true begin
    include("cm_helpers.jl")
    include("cm_ac01_04.jl")
    include("cm_ac05_06.jl")
    include("cm_ac07_12.jl")
    include("cm_ac13_16.jl")
    include("cm_ac17_19.jl")
    include("cm_d105.jl")
end
