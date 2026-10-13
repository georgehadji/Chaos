"""
    Systems

Maps with analytic tangents: standard map, linear symplectic maps (standard-map).
"""
module Systems

using StaticArrays
using LinearAlgebra

export AbstractMap, StandardMap, LinearSymplecticMap
export rotation_block, hyperbolic_block, block_diagonal, conjugated
export dimension, forward, inverse, tangent, inverse_tangent
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

include("standard_map.jl")
include("linear_symplectic.jl")
include("residuals.jl")

end
