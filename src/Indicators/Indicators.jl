"""
    Indicators

Chaos indicators as reducers of the tangent stream, and the pure functions behind them: the Smaller Alignment Index
SALI (sali), the Fast Lyapunov Indicator FLI (fli), the Mean Exponential Growth factor of Nearby Orbits MEGNO (megno) and the Relative Lyapunov Indicator RLI (rli) and the stretching-number spectra with their spectral distance (stretching). Imports Foundation, Systems and Tangent; never a Lyapunov module (a peer layer). The shared norm `safe_norm` is in Foundation (D69).
"""
module Indicators

using ..Foundation
using ..Systems
using ..Tangent
using StaticArrays

import ..Tangent: reducer_init, reducer_update, reducer_final, reducer_step, reducer_failed, reducer_error, reducer_provenance, reducer_workspace

export SALI, alignment_indices, sali_pair, sali_from_log, sali_floor, check_frame, FLI, fli_unlogged, MEGNO, RLI, pair_start, pair_distance
export StretchingNumbers, spectral_distance, spectrum_cell, spectrum_density, stretching_range, check_run

include("sali.jl")
include("alignment.jl")
include("fli.jl")
include("megno.jl")
include("rli.jl")
include("stretching.jl")

end
