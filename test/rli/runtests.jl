# rli suite (module-build stage 02, RED): the Relative Lyapunov Indicator, catalog row I9.
# Spec: workspaces/module-build/stages/01-spec/output/rli-spec.md (approved at the stage 01 checkpoint, D75: the orbit and its shadow as one orbit of the product of the map with itself,
# a new additive construct in the Systems layer, on the existing tangent stream with k = 2 and a supplied block frame, and a small column reducer with O(1) state in the Indicators layer;
# separation along xi_0, Euclidean minimal-image norm, realized Dx_r recorded; mean over the handoffs; wrapped state by default; raw DL_k and smoothed RLI delivered).
# API proposal (Systems: PairMap, pair_start, pair_distance; Indicators: RLI and its check_frame method): see workspaces/module-build/stages/02-tests/output/rli-test-plan.md.
# Criteria 14 (ChaosTools two one-column tangent systems) and 15 (CPU vs GPU, D68 form) are stage 04 (verification environment, local GPU); the test plan maps them.
# Every testset begins by naming the new symbols it uses, so a missing Chaoslab symbol errors that testset only (UndefVarError), and no top-level code calls into the new API.
# Package-suite sizes follow D38 (3 minute ceiling); verify/rli.jl sets RL_FULL_RUN = true for the full sizes.
include(joinpath(@__DIR__, "rl_helpers.jl"))
include(joinpath(@__DIR__, "rl_ac01_03.jl"))
include(joinpath(@__DIR__, "rl_ac04_07.jl"))
include(joinpath(@__DIR__, "rl_ac08_10.jl"))
include(joinpath(@__DIR__, "rl_ac11_13.jl"))
include(joinpath(@__DIR__, "rl_pair.jl"))
include(joinpath(@__DIR__, "rl_api.jl"))
include(joinpath(@__DIR__, "rl_records.jl"))
