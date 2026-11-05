# gali suite (module-build stage 02, RED). Spec: workspaces/module-build/stages/01-spec/output/gali-spec.md (approved at the checkpoint, D106).
# API proposal (layer Chaoslab.Indicators, the reducer protocol of the tangent stream): see workspaces/module-build/stages/02-tests/output/gali-test-plan.md.
# AC19 (ChaosTools oracle), AC21 (local GPU) and the timing part of AC20 are stage 04 (verification environment); the test plan maps them.
# Every testset begins by naming the new symbols it uses (GALI, gali_values, gali_value, gali_floor, check_frame methods), so a missing module errors that testset
# (UndefVarError) and no top-level code calls into the module. Package-suite sizes follow D38 (3 minute ceiling); verify/gali.jl sets GA_FULL_RUN = true for the full sizes.
include(joinpath(@__DIR__, "ga_helpers.jl"))
include(joinpath(@__DIR__, "ga_ac01_04.jl"))
include(joinpath(@__DIR__, "ga_ac05_08.jl"))
include(joinpath(@__DIR__, "ga_ac09_12.jl"))
include(joinpath(@__DIR__, "ga_ac13_16.jl"))
include(joinpath(@__DIR__, "ga_ac17_20.jl"))
include(joinpath(@__DIR__, "ga_d108.jl"))
include(joinpath(@__DIR__, "ga_d109.jl"))
