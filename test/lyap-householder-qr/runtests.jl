# lyap-householder-qr suite (stage 02, RED). Spec: workspaces/module-build/stages/01-spec/output/lyap-householder-qr-spec.md
# API: Batch B contract (decision D34). H1 (treppen path) and AC16 are out of M2 (D30); AC14 and AC17 are stage 04
# (see workspaces/module-build/stages/02-tests/output/lyap-householder-qr-test-plan.md). No top-level call into Chaoslab.Lyapunov:
# a missing module errors the testsets, nothing else.
include(joinpath(@__DIR__, "hh_helpers.jl"))
include(joinpath(@__DIR__, "hh_factor.jl"))
include(joinpath(@__DIR__, "hh_exact.jl"))
include(joinpath(@__DIR__, "hh_panel.jl"))
include(joinpath(@__DIR__, "hh_cross.jl"))
include(joinpath(@__DIR__, "hh_runtime.jl"))
include(joinpath(@__DIR__, "hh_backend.jl"))
