"""
    Indicators

Chaos indicators as reducers of the tangent stream, and the pure functions behind them: the Smaller Alignment Index
SALI (sali) and the Fast Lyapunov Indicator FLI (fli). Imports Foundation, Systems and Tangent; never a Lyapunov module (a peer layer), so the scaled norm is its own.
"""
module Indicators

using ..Foundation
using ..Systems
using ..Tangent
using StaticArrays

import ..Tangent: reducer_init, reducer_update, reducer_final, reducer_step, reducer_failed, reducer_error, reducer_provenance

export SALI, alignment_indices, sali_pair, sali_from_log, sali_floor, check_frame, FLI, fli_unlogged

include("sali.jl")
include("alignment.jl")
include("fli.jl")

end
