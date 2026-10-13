# tangent-core suite (stage 02, RED). Spec: workspaces/module-build/stages/01-spec/output/tangent-core-spec.md
# AC12 to AC14 and the two per-chunk fallback reasons are stage 04 (see output/tangent-core-test-plan.md).
include(joinpath(@__DIR__, "tc_helpers.jl"))
include(joinpath(@__DIR__, "tc_ac01_04.jl"))
include(joinpath(@__DIR__, "tc_ac05_08.jl"))
include(joinpath(@__DIR__, "tc_ac09_11.jl"))
include(joinpath(@__DIR__, "tc_state_length.jl"))
