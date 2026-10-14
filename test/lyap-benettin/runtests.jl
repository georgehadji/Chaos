# lyap-benettin suite (stage 02, RED). Spec: workspaces/module-build/stages/01-spec/output/lyap-benettin-spec.md
# API contract: M2 Batch B api-batch-b.md (D34). AC12 to AC16 are stage 04 (verification environment, local GPU); see
# workspaces/module-build/stages/02-tests/output/lyap-benettin-test-plan.md. Every testset builds what it needs inside the block,
# so a missing Chaoslab name (Benettin, safe_norm, effective_eps, ...) errors that testset only.
include(joinpath(@__DIR__, "lb_helpers.jl"))
include(joinpath(@__DIR__, "lb_ac01_04.jl"))
include(joinpath(@__DIR__, "lb_ac05_07.jl"))
include(joinpath(@__DIR__, "lb_ac08_11.jl"))
include(joinpath(@__DIR__, "lb_ac17_19.jl"))
