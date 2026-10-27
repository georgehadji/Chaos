# weighted-birkhoff suite (module-build stage 02, RED): the weighted Birkhoff average of K observables of the orbit point over H horizons, the pair of consecutive windows, the digits index (catalog I13,
# D80, D83), as an O(1) stream reader (WeightedBirkhoff) and a pure function of a stored series (weighted_birkhoff) with the index birkhoff_index.
# Spec: workspaces/module-build/stages/01-spec/output/weighted-birkhoff-spec.md (approved at the stage 01 checkpoint, D83: a new reader reducer for K observables and H horizons, the pure function, the
# index, the observable rules; the verified weight of NF4, the compensated sums and the run protocol are reused; weighted FLI, extrapolation and thresholds are out of scope).
# API proposal (Indicators: weighted_birkhoff, birkhoff_index, WeightedBirkhoff, CosObservable, SinObservable, LiftedObservable, observable_name, the check_run method; the verified birkhoff_weight(s),
# WeightedPhaseAdvance, FLI, SALI, MEGNO are used as they are): see workspaces/module-build/stages/02-tests/output/weighted-birkhoff-test-plan.md.
# Criteria 12 (V) and 13 (G) are encoded as far as a package suite can: 12(a), 12(b) run here (the BigFloat Poisson quadrature and the BigFloat orbit); 12(c) has nothing to test; 13 runs only when a device
# is present (verify/weighted-birkhoff.jl at stage 04 loads CUDA first) and is skipped otherwise.
# Every testset begins by naming the new symbols it uses, so a missing Chaoslab symbol errors that testset only (UndefVarError), and no top-level code calls into the new API; wb_selftest.jl uses no new
# name and passes now. Package-suite sizes follow D38 (a few minutes); verify/weighted-birkhoff.jl sets WB_FULL_RUN = true for the sizes of the spec.
include(joinpath(@__DIR__, "wb_helpers.jl"))
include(joinpath(@__DIR__, "wb_selftest.jl"))
include(joinpath(@__DIR__, "wb_ac01_03.jl"))
include(joinpath(@__DIR__, "wb_ac04_06.jl"))
include(joinpath(@__DIR__, "wb_ac07_09.jl"))
include(joinpath(@__DIR__, "wb_ac10_13.jl"))
include(joinpath(@__DIR__, "wb_api.jl"))
include(joinpath(@__DIR__, "wb_records.jl"))
