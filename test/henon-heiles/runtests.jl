# henon-heiles suite (module-build stage 02, RED). Spec: workspaces/module-build/stages/01-spec/output/henon-heiles-spec.md (approved at the checkpoint, D116).
# API proposal: workspaces/module-build/stages/02-tests/output/henon-heiles-test-plan.md. Every testset begins by naming the new symbols it uses, so a missing module
# errors that testset (UndefVarError); no top-level code and no helper calls the module at include time. Package sizes follow D38 (3 minute ceiling); verify/henon-heiles.jl
# sets HH_FULL_RUN = true before including this file for the full sizes. Stage 04 only: the Taylor and ChaosTools spectra (O3, O4), the GPU part of AC16, the full-size portraits.
include(joinpath(@__DIR__, "hh_helpers.jl"))
include(joinpath(@__DIR__, "hh_ac01_05.jl"))
include(joinpath(@__DIR__, "hh_ac06_09.jl"))
include(joinpath(@__DIR__, "hh_ac10_12.jl"))
include(joinpath(@__DIR__, "hh_ac13_17.jl"))
