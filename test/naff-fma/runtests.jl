# naff-fma suite (module-build stage 02, RED): the windowed Fourier frequency (NAFF-type, E1) of two consecutive windows of an orbit series, the frequency map distance and its decimal logarithm (catalog I12,
# I12b, D80 with Amendment D80: the cosine window family of Laskar2003, window order as an option), and the weighted mean phase advance (E2) as an O(1) stream reducer.
# Spec: workspaces/module-build/stages/01-spec/output/naff-fma-spec.md (approved at the stage 01 checkpoint, D80: P1 option (a2) for E1 as a pure function of a series supplied by Experiments, option (c)
# for E2 as a reader reducer with the total length given in advance, no Float64 buffer in the Tangent workspace; Amendment D80 sources the iteration and the windows from Laskar2003, owner confirmed D80c).
# API proposal (Indicators: window_weights, window_frequency, frequency_map, check_frequency_map, circle_distance, frequency_distance, frequency_index, birkhoff_weight, birkhoff_weights, phase_increments,
# weighted_phase_frequency, WeightedPhaseAdvance, PlanarObservable, LiftedAngleObservable, and the check_run method): see workspaces/module-build/stages/02-tests/output/naff-fma-test-plan.md.
# Criteria 14 (independent orbit, ChaosTools, Python NAFF), 15 (CPU vs GPU) are stage 04; criterion 16 (the series feed of Experiments) is conditional on the feed being built and is not written here (plan,
# Items for the owner 2). Criterion 7(c) is a statement about the caller's series builder and is gated in nf_selftest.jl (green now).
# Every testset begins by naming the new symbols it uses, so a missing Chaoslab symbol errors that testset only (UndefVarError), and no top-level code calls into the new API.
# Package-suite sizes follow D38 (a few minutes); verify/naff-fma.jl sets NF_FULL_RUN = true for the full sizes.
include(joinpath(@__DIR__, "nf_helpers.jl"))
include(joinpath(@__DIR__, "nf_selftest.jl"))
include(joinpath(@__DIR__, "nf_ac01_03.jl"))
include(joinpath(@__DIR__, "nf_ac04_06.jl"))
include(joinpath(@__DIR__, "nf_ac07_09.jl"))
include(joinpath(@__DIR__, "nf_ac10_13.jl"))
include(joinpath(@__DIR__, "nf_d80.jl"))
include(joinpath(@__DIR__, "nf_api.jl"))
include(joinpath(@__DIR__, "nf_records.jl"))
include(joinpath(@__DIR__, "nf_d81.jl"))
