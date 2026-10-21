# sali suite (module-build stage 02, RED). Spec: workspaces/module-build/stages/01-spec/output/sali-spec.md (approved at the checkpoint, D66).
# API proposal (layer Chaoslab.Indicators, the reducer protocol of the Lyapunov layer): see workspaces/module-build/stages/02-tests/output/sali-test-plan.md.
# AC16 to AC18 are stage 04 (verification environment: ChaosTools GALI_2 oracle, strong chaos at long horizon, local GPU); the test plan maps them.
# Every testset begins by naming the Indicators symbols it uses, so a missing Chaoslab.Indicators errors that testset only (UndefVarError), and no
# top-level code calls into the module. Package-suite sizes follow D38 (3 minute ceiling); verify/sali.jl sets SL_FULL_RUN = true for the full sizes.
include(joinpath(@__DIR__, "sl_helpers.jl"))
include(joinpath(@__DIR__, "sl_ac01_03.jl"))
include(joinpath(@__DIR__, "sl_ac04_06.jl"))
include(joinpath(@__DIR__, "sl_ac07_09.jl"))
include(joinpath(@__DIR__, "sl_ac10_12.jl"))
include(joinpath(@__DIR__, "sl_ac13_15.jl"))
