# megno suite (module-build stage 02, RED): the Mean Exponential Growth factor of Nearby Orbits, catalog row I8.
# Spec: workspaces/module-build/stages/01-spec/output/megno-spec.md (approved at the stage 01 checkpoint, D73: own reducer with the D72 replace option; Y and Ybar in the
# (1,-1) and (2,0) forms for both columns of the seeded frame, column 1 primary, no combination of the columns).
# API proposal (layer Chaoslab.Indicators, the reducer protocol of the Tangent layer): see workspaces/module-build/stages/02-tests/output/megno-test-plan.md.
# Criteria 14 (ChaosTools one-column tangent oracle line) and 15 (CPU vs GPU, D68 form) are stage 04 (verification environment, local GPU); the test plan maps them.
# Every testset begins by naming the Indicators symbols it uses, so a missing Chaoslab.Indicators errors that testset only (UndefVarError), and no top-level
# code calls into the module. Package-suite sizes follow D38 (3 minute ceiling); verify/megno.jl sets MG_FULL_RUN = true for the full sizes.
include(joinpath(@__DIR__, "mg_helpers.jl"))
include(joinpath(@__DIR__, "mg_ac01_05.jl"))
include(joinpath(@__DIR__, "mg_ac06_08.jl"))
include(joinpath(@__DIR__, "mg_ac09_10.jl"))
include(joinpath(@__DIR__, "mg_ac11_13.jl"))
include(joinpath(@__DIR__, "mg_api.jl"))
