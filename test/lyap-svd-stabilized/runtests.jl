# lyap-svd-stabilized suite (stage 02, RED). Spec: workspaces/module-build/stages/01-spec/output/lyap-svd-stabilized-spec.md
# API: scratchpad m2/api-batch-b.md (decision D34). AC10 and AC11 are stage 04 (see output/lyap-svd-stabilized-test-plan.md).
# Shared helpers (theta_sum, references, runners) come from the lyap-gram-schmidt folder, included once.
@isdefined(GS_HELPERS_LOADED) || include(joinpath(@__DIR__, "..", "lyap-gram-schmidt", "gs_helpers.jl"))
include(joinpath(@__DIR__, "svd_helpers.jl"))
@testset "lyap-svd-stabilized" begin
    include(joinpath(@__DIR__, "svd_ac01_05.jl"))
    include(joinpath(@__DIR__, "svd_ac06_12.jl"))
end
