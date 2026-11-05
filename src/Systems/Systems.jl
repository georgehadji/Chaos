"""
    Systems

Maps with analytic tangents: standard map, linear symplectic maps (standard-map), the product of a map with itself (rli),
the Froeschle 4D map and the chain of coupled standard maps (coupled-maps).
"""
module Systems

using StaticArrays
using LinearAlgebra

export AbstractMap, StandardMap, LinearSymplecticMap, WrappedStandardMap, torus_state, lift
export rotation_block, hyperbolic_block, block_diagonal, conjugated
export dimension, forward, inverse, tangent, inverse_tangent
export PairMap, PairTangent
export FroeschleMap, WrappedFroeschleMap, StandardMapChain, WrappedStandardMapChain, ChainTangent
export forward!, inverse!, tangent!, inverse_tangent!, tangent_action!, inverse_tangent_action!
export wrap, symplectic_form, det_residual, symplectic_residual

"Supertype of the maps of this module; subtypes define `forward`, `inverse`, `tangent`, `inverse_tangent`, `dimension`."
abstract type AbstractMap end

"""
    dimension(sys) -> Int

Phase-space dimension of `sys`.
"""
function dimension end

"""
    forward(sys, z)

One lifted step of `sys` from the state `z` (an `SVector` of any element type).
"""
function forward end

"""
    inverse(sys, z)

One inverse step of `sys`: `inverse(sys, forward(sys, z)) == z` up to rounding.
"""
function inverse end

"""
    tangent(sys, z)

Jacobian of one forward step, evaluated at the pre-step point `z`.
"""
function tangent end

"""
    inverse_tangent(sys, z_image)

Jacobian of one inverse step, evaluated at the image point `z_image = forward(sys, z)`.
"""
function inverse_tangent end

"""
    forward!(out, sys, z)

In-place [`forward`](@ref) for the chain maps: fills and returns `out` (not aliasing `z`), no allocation after warm-up.
"""
function forward! end

"In-place [`inverse`](@ref) for the chain maps; same contract as [`forward!`](@ref)."
function inverse! end

"In-place [`tangent`](@ref): fills and returns the dense matrix `M`."
function tangent! end

"In-place [`inverse_tangent`](@ref): fills and returns the dense matrix `M`."
function inverse_tangent! end

"""
    tangent_action!(out, sys, z, W)

`out = M W` for the tangent `M` of one step at `z` and a `2N x k` frame `W`, using the banded structure (cost linear in `N`).
"""
function tangent_action! end

"`out = M^-1 W` for the inverse tangent at the image `z`; see [`tangent_action!`](@ref)."
function inverse_tangent_action! end

include("standard_map.jl")
include("wrapped_standard_map.jl")
include("linear_symplectic.jl")
include("pair.jl")
include("froeschle_map.jl")
include("standard_map_chain.jl")
include("chain_tangent.jl")
include("wrapped_standard_map_chain.jl")
include("residuals.jl")

end
