# ofli suite (module-build stage 02, RED): the Orthogonal Fast Lyapunov Indicator for Hamiltonian flows, catalog row I6 (OFLI2, row I7, is catalog only).
# Spec: workspaces/module-build/stages/01-spec/output/ofli-spec.md (approved at the checkpoint, D118, open points P1 to P12 as proposed).
# API proposal (layer Chaoslab.Indicators, the reducer protocol of the Tangent layer): see workspaces/module-build/stages/02-tests/output/ofli-test-plan.md.
# Criterion 17 (Taylor reference, verification environment) and the device part of criterion 18 (local GPU) are stage 04; the test plan maps them, and the package
# suite carries a Double64 stand-in for 17 and the CPU part of 18. Every testset begins by naming the new symbols it uses, so a missing module errors that testset
# only (UndefVarError), and no top-level code and no helper calls the module at include time. Package sizes follow D38 (3 minute ceiling); verify/ofli.jl sets
# OF_FULL_RUN = true before including this file for the full sizes of the spec.
include(joinpath(@__DIR__, "of_helpers.jl"))
include(joinpath(@__DIR__, "of_ac01_05.jl"))
include(joinpath(@__DIR__, "of_ac06_10.jl"))
include(joinpath(@__DIR__, "of_ac11_15.jl"))
include(joinpath(@__DIR__, "of_ac16_18.jl"))
