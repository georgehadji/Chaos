# experiments suite (stage 02, RED). Spec: workspaces/module-build/stages/01-spec/output/experiments-spec.md
# API under test: submodule Chaoslab.Experiments (names listed in workspaces/module-build/stages/02-tests/output/experiments-test-plan.md).
# Nothing here calls Experiments at load time: a missing module errors the testsets that need it and nothing else.
include(joinpath(@__DIR__, "ex_helpers.jl"))
include(joinpath(@__DIR__, "ex_inputs.jl"))
include(joinpath(@__DIR__, "ex_geometry.jl"))
include(joinpath(@__DIR__, "ex_periodic.jl"))
include(joinpath(@__DIR__, "ex_screening.jl"))
include(joinpath(@__DIR__, "ex_predicates.jl"))
include(joinpath(@__DIR__, "ex_draw.jl"))
include(joinpath(@__DIR__, "ex_subsets.jl"))
include(joinpath(@__DIR__, "ex_panel.jl"))
include(joinpath(@__DIR__, "ex_determinism.jl"))
