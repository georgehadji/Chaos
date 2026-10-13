"""
    Tangent

Propagation of orbits and deviation vectors, one KernelAbstractions kernel (tangent-core).
"""
module Tangent

using ..Foundation
using ..Foundation: AbstractDevice, NoDevice, DEVICE, free_memory, is_oom, is_compile_error, is_device_lost, gpu_step_ns   # device hooks, not exported
using ..Systems
using StaticArrays
using LinearAlgebra
using KernelAbstractions
using StableRNGs

export initial_deviations, reducer_init, reducer_update, reducer_final, reducer_step, reducer_failed, reducer_error
export PropagationResult, propagate, propagate_orbit, chunk_size, bytes_per_orbit

"""
    reducer_init(r, T, Val(d), Val(k)) -> state

Initial per-orbit state of reducer `r`: an immutable isbits value (GPU-compatible).
"""
function reducer_init end

"""
    reducer_update(r, state, n, x, W) -> (state, W)

Called at every handoff (steps tau, 2 tau, ...) with the step index `n`, point `x` and deviation matrix `W`.
May return a replaced `W`; later reducers in the chain see the replacement (spec E7).
"""
function reducer_update end

"""
    reducer_final(r, state) -> NamedTuple

Host-side summary of the final state of reducer `r` for one orbit.
"""
function reducer_final end

"""
    reducer_step(r, state, n, x, W) -> (state, W)

Entry point the per-orbit body calls at each handoff; defaults to [`reducer_update`](@ref). A reducer that must
not throw inside the kernel overrides it to record a failure in its state (see [`reducer_failed`](@ref)).
"""
reducer_step(r, state, n, x, W) = reducer_update(r, state, n, x, W)

"""
    reducer_failed(r, state) -> Bool

Whether the state records a health failure (default `false`). A failed orbit stops at that handoff.
"""
reducer_failed(r, state) = false

"""
    reducer_error(r, state, orbit, backend) -> NumericalHealthError

Host-side error for a state with `reducer_failed(r, state)`, carrying the orbit id and backend.
"""
function reducer_error end

include("initial.jl")
include("orbit.jl")
include("propagate.jl")

end
