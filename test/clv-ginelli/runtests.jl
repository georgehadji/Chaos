# clv-ginelli suite (module-build stage 02, RED): covariant Lyapunov vectors by the Ginelli method (forward QR with stored R, backward back-substitution), catalog row L12.
# Spec: workspaces/module-build/stages/01-spec/output/clv-ginelli-spec.md (approved at the checkpoint, D122, open points P1 to P16 as proposed).
# API proposal (layer Chaoslab.Lyapunov, its own per-orbit loop, not a reducer of the tangent stream): see workspaces/module-build/stages/02-tests/output/clv-ginelli-test-plan.md.
# Criterion 13 (static null-space cross-check, 256 bits) and the recorded lines R1 to R6 are stage 04 (verify/clv-ginelli.jl); the plan maps them. Every testset begins by naming the new
# symbols it uses, so a missing module errors that testset only (UndefVarError), and no top-level code and no helper calls the module at include time. Package sizes follow D38
# (3 minute ceiling for the folder); verify/clv-ginelli.jl sets CLV_FULL_RUN = true before including this file for the full sizes of the spec.
include(joinpath(@__DIR__, "clv_helpers.jl"))
include(joinpath(@__DIR__, "clv_ac01_05.jl"))
include(joinpath(@__DIR__, "clv_ac06_10.jl"))
include(joinpath(@__DIR__, "clv_ac11_14.jl"))
