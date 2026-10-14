"""
    Lyapunov

Lyapunov-exponent estimators as reducers of the tangent stream: Benettin (lyap-benettin) and its shared helpers.
"""
module Lyapunov

using ..Foundation
using ..Systems
using ..Tangent
using StaticArrays
using LinearAlgebra

import ..Tangent: reducer_init, reducer_update, reducer_final, reducer_step, reducer_failed, reducer_error, reducer_provenance

export safe_norm, Benettin

include("norm.jl")
include("benettin.jl")

end
