# analysis suite (stage 02, RED). Spec: workspaces/module-build/stages/01-spec/output/analysis-spec.md
# API under test: submodule Chaoslab.Analysis (names listed in workspaces/module-build/stages/02-tests/output/analysis-test-plan.md).
# Nothing here calls Analysis at load time: a missing module errors the testsets that need it and nothing else.
include(joinpath(@__DIR__, "an_helpers.jl"))
include(joinpath(@__DIR__, "an_fixture.jl"))      # fixture self-test: passes before the module exists
include(joinpath(@__DIR__, "an_params.jl"))
include(joinpath(@__DIR__, "an_h1.jl"))
include(joinpath(@__DIR__, "an_sets.jl"))
include(joinpath(@__DIR__, "an_bootstrap.jl"))
include(joinpath(@__DIR__, "an_h4.jl"))
include(joinpath(@__DIR__, "an_labels.jl"))
include(joinpath(@__DIR__, "an_h5.jl"))
include(joinpath(@__DIR__, "an_pipeline.jl"))
include(joinpath(@__DIR__, "an_cli.jl"))
include(joinpath(@__DIR__, "an_breakdown.jl"))     # regression of stage 04: breakdown_tau items on the real store
include(joinpath(@__DIR__, "an_review.jl"))        # regressions of the stage 04 review
