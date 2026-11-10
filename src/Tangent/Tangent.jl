"""
    Tangent

Propagation of orbits and deviation vectors, one KernelAbstractions kernel (tangent-core).
"""
module Tangent

using ..Foundation
using ..Foundation: fused_product   # frame product with explicit fma (D115), not exported
using ..Foundation: AbstractDevice, NoDevice, DEVICE, free_memory, is_oom, is_compile_error, is_device_lost, gpu_step_ns   # device hooks, not exported
using ..Systems
using StaticArrays
using LinearAlgebra
using KernelAbstractions
using StableRNGs

export initial_deviations, reducer_init, reducer_update, reducer_final, reducer_step, reducer_failed, reducer_error, reducer_provenance
export reducer_workspace, reducer_bytes
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

"""
    reducer_provenance(r) -> NamedTuple

Provenance keys reducer `r` adds to the run record (default `NamedTuple()`); `propagate` merges them over the chain
in order, later reducers winning on a shared key (D34 A2).
"""
reducer_provenance(r) = NamedTuple()

"""
    reducer_workspace(r, T, Val(d), Val(k)) -> Int

Number of `UInt32` entries reducer `r` wants per orbit (default `0`). The workspace is allocated once per chunk, zeroed at
the first step of a run, carried across relaunches with the state and read back by the three-argument
[`reducer_final`](@ref). At most one reducer of a chain may request one; it is handed to every reducer of the chain
(additive amendment of the stretching-numbers run, P1 (c); the spec section is appended to the tangent-core family).
"""
reducer_workspace(r, ::Type, ::Val, ::Val) = 0

"""
    reducer_step(r, state, n, x, W, ws) -> (state, W)

Six-argument form for a reducer with a workspace: `ws` is the orbit's own `AbstractVector{UInt32}` slice of length
[`reducer_workspace`](@ref), written in place. The default forwards to the five-argument form, so a reducer without a
workspace is unchanged. [`reducer_update`](@ref) has the same six-argument form (the throwing one).
"""
reducer_step(r, state, n, x, W, ws) = reducer_step(r, state, n, x, W)

reducer_update(r, state, n, x, W, ws) = reducer_update(r, state, n, x, W)

"""
    reducer_final(r, state, ws) -> NamedTuple

Host summary of a reducer with a workspace (`ws` the orbit's slice, never kept: the summary holds copies). The default
forwards to the two-argument form.
"""
reducer_final(r, state, ws) = reducer_final(r, state)

"""
    reducer_bytes(r, T, Val(d), Val(k)) -> Int

Reducer bytes per orbit of the memory formula: `sizeof` of the initial state when it is isbits (else 0) plus 4 bytes per
workspace entry.
"""
function reducer_bytes(r, ::Type{T}, v::Val, k::Val) where {T}
    s = reducer_init(r, T, v, k)
    return (isbitstype(typeof(s)) ? sizeof(s) : 0) + 4 * reducer_workspace(r, T, v, k)
end

include("initial.jl")
include("orbit.jl")
include("propagate.jl")

end
