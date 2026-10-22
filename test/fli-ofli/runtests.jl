# fli-ofli suite (module-build stage 02, RED): the Fast Lyapunov Indicator, catalog row I5 (OFLI, row I6, is not in this module, D71).
# Spec: workspaces/module-build/stages/01-spec/output/fli-ofli-spec.md (approved at the checkpoint, D71: FLI only, chain option (b) with the one-replacer rule).
# API proposal (layer Chaoslab.Indicators, the reducer protocol of the Tangent layer): see workspaces/module-build/stages/02-tests/output/fli-ofli-test-plan.md.
# Criteria 15 (ChaosTools gali oracle line) and 16 (CPU vs GPU, D68 form) are stage 04 (verification environment, local GPU); the test plan maps them.
# Every testset begins by naming the Indicators symbols it uses, so a missing Chaoslab.Indicators errors that testset only (UndefVarError), and no top-level
# code calls into the module. Package-suite sizes follow D38 (3 minute ceiling); verify/fli-ofli.jl sets FL_FULL_RUN = true for the full sizes.
include(joinpath(@__DIR__, "fl_helpers.jl"))
include(joinpath(@__DIR__, "fl_ac01_05.jl"))
include(joinpath(@__DIR__, "fl_ac06_08.jl"))
include(joinpath(@__DIR__, "fl_ac09_11.jl"))
include(joinpath(@__DIR__, "fl_ac12_14.jl"))
include(joinpath(@__DIR__, "fl_api.jl"))
