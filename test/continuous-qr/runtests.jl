# continuous-qr suite (module-build stage 02, RED): the continuous QR method for Hamiltonian flows, catalog row L8 (continuous SVD, L10, is catalog only).
# Spec: workspaces/module-build/stages/01-spec/output/continuous-qr-spec.md (approved at the checkpoint, D120, open points P1 to P15 as proposed).
# API proposal (layer Chaoslab.Lyapunov, its own per-orbit loop, not a reducer of the tangent stream): see workspaces/module-build/stages/02-tests/output/continuous-qr-test-plan.md.
# Criterion 21 (ChaosTools cross-check) and the Taylor form of 20 (verification environment), and the device part of 22 (local GPU) are stage 04; the test plan maps them,
# and the package suite carries a Double64 stand-in for 20 and the CPU part of 22. Every testset begins by naming the new symbols it uses, so a missing module errors that
# testset only (UndefVarError), and no top-level code and no helper calls the module at include time. Package sizes follow D38 (3 minute ceiling for the folder);
# verify/continuous-qr.jl sets CQ_FULL_RUN = true before including this file for the full sizes of the spec.
include(joinpath(@__DIR__, "cq_helpers.jl"))
include(joinpath(@__DIR__, "cq_ac01_05.jl"))
include(joinpath(@__DIR__, "cq_ac06_10.jl"))
include(joinpath(@__DIR__, "cq_ac11_15.jl"))
include(joinpath(@__DIR__, "cq_ac16_22.jl"))
