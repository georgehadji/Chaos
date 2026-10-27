# stretching suite (module-build stage 02, RED): the stretching numbers of one orbit, their dynamical spectra (SSN) and the spectral distance D between the spectra of two unit deviation
# vectors, catalog row I10 (and L13 as the same histogram object).
# Spec: workspaces/module-build/stages/01-spec/output/stretching-numbers-spec.md (approved at the stage 01 checkpoint, D78: a two-column reducer with cumulative spectra of B cells over a
# derived range, overflow cells, D every interval; the counts live in a per-orbit unsigned-integer workspace requested by the reducer, an additive amendment to the Tangent layer, P1 (c)).
# API proposal (Indicators: StretchingNumbers, spectral_distance, spectrum_cell, spectrum_density, stretching_range, check_run, check_frame method; Tangent: the workspace surface
# reducer_workspace, reducer_bytes, the six-argument reducer_step, the three-argument reducer_final, the `workspace` keyword of propagate_orbit): see
# workspaces/module-build/stages/02-tests/output/stretching-numbers-test-plan.md.
# Criteria 14 (ChaosTools independent tangent line), 15 (ChaosTools local growth line) and 16 (CPU vs GPU, D68 form) are stage 04 (verification environment, local GPU); the plan maps them.
# Every testset begins by naming the new symbols it uses, so a missing Chaoslab symbol errors that testset only (UndefVarError), and no top-level code calls into the new API.
# Package-suite sizes follow D38 (3 minute class); verify/stretching-numbers.jl sets SN_FULL_RUN = true for the full sizes.
include(joinpath(@__DIR__, "sn_helpers.jl"))
include(joinpath(@__DIR__, "sn_selftest.jl"))
include(joinpath(@__DIR__, "sn_ac01_03.jl"))
include(joinpath(@__DIR__, "sn_ac04_07.jl"))
include(joinpath(@__DIR__, "sn_ac08_10.jl"))
include(joinpath(@__DIR__, "sn_ac11_13.jl"))
include(joinpath(@__DIR__, "sn_api.jl"))
include(joinpath(@__DIR__, "sn_records.jl"))
