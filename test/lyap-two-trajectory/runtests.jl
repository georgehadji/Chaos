# lyap-two-trajectory suite (stage 02, RED). Spec: workspaces/module-build/stages/01-spec/output/lyap-two-trajectory-spec.md
# API: Batch B contract (decision D34): propagate_shadow, propagate_shadow_orbit, separation, displaced; own per-orbit loop (D30).
# AC13 to AC17 are stage 04 (see workspaces/module-build/stages/02-tests/output/lyap-two-trajectory-test-plan.md). No top-level call into
# Chaoslab.Lyapunov: a missing module errors the testsets, nothing else.
# Package suite size (D38, under 3 minutes): L2_FULL (l2_helpers.jl) is false here; verify/lyap-two-trajectory.jl runs the same files at full size.
include(joinpath(@__DIR__, "l2_helpers.jl"))
include(joinpath(@__DIR__, "l2_analytic.jl"))
include(joinpath(@__DIR__, "l2_sweeps.jl"))
include(joinpath(@__DIR__, "l2_runtime.jl"))
# Amendment D44 (step-chunked relaunch, item 15, D44k): output/tangent-core-d44-test-plan.md
include(joinpath(@__DIR__, "l2_d44.jl"))
