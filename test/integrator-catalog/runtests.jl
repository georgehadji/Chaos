# integrator-catalog suite (module-build stage 02, RED). Spec: workspaces/module-build/stages/01-spec/output/integrator-catalog-spec.md (approved at the checkpoint, D111).
# API proposal: workspaces/module-build/stages/02-tests/output/integrator-catalog-test-plan.md. Every testset begins by naming the new symbols it uses, so a missing module
# errors that testset (UndefVarError); no top-level code and no helper calls the module at include time. Package sizes follow D38 (3 minute ceiling); verify/integrator-catalog.jl sets
# IC_FULL_RUN = true before including this file for the full sizes. Stage 04 only: AC3(c) (O6), AC3(d) (O8), the print comparisons AC3(b2), AC18 GPU part, AC15(a)(d) at full size.
include(joinpath(@__DIR__, "ic_helpers.jl"))
include(joinpath(@__DIR__, "ic_flows.jl"))
include(joinpath(@__DIR__, "ic_ac01_03.jl"))
include(joinpath(@__DIR__, "ic_ac04_08.jl"))
include(joinpath(@__DIR__, "ic_ac09_12.jl"))
include(joinpath(@__DIR__, "ic_ac13_16.jl"))
include(joinpath(@__DIR__, "ic_ac17_23.jl"))
