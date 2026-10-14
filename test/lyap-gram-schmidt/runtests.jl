# lyap-gram-schmidt suite (stage 02, RED). Spec: workspaces/module-build/stages/01-spec/output/lyap-gram-schmidt-spec.md
# API: scratchpad m2/api-batch-b.md (decision D34). AC12 and AC13 are stage 04 (see output/lyap-gram-schmidt-test-plan.md).
# The helper file is shared with test/lyap-svd-stabilized; it defines functions only, no top-level Lyapunov call.
# The outer testset keeps a failing criterion from aborting the rest when this file is included on its own.
@isdefined(GS_HELPERS_LOADED) || include(joinpath(@__DIR__, "gs_helpers.jl"))
@testset "lyap-gram-schmidt" begin
    include(joinpath(@__DIR__, "gs_api.jl"))
    include(joinpath(@__DIR__, "gs_ac01_05.jl"))
    include(joinpath(@__DIR__, "gs_ac06_09.jl"))
    include(joinpath(@__DIR__, "gs_ac10_14.jl"))
end
