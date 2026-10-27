"""
    Indicators

Chaos indicators as reducers of the tangent stream, and the pure functions behind them: the Smaller Alignment Index
SALI (sali), the Fast Lyapunov Indicator FLI (fli), the Mean Exponential Growth factor of Nearby Orbits MEGNO (megno) and the Relative Lyapunov Indicator RLI (rli), the stretching-number spectra with their spectral distance (stretching), the frequency map analysis of an orbit series (naff: the windowed Fourier frequency, the frequency map distance and the weighted phase advance reducer) and the weighted Birkhoff averages of observables of the orbit (wbirkhoff: the pure function, the pair index and the reader reducer). Imports Foundation, Systems and Tangent; never a Lyapunov module (a peer layer). The shared norm `safe_norm` is in Foundation (D69).
"""
module Indicators

using ..Foundation
using ..Systems
using ..Tangent
using StaticArrays

import ..Tangent: reducer_init, reducer_update, reducer_final, reducer_step, reducer_failed, reducer_error, reducer_provenance, reducer_workspace

export SALI, alignment_indices, sali_pair, sali_from_log, sali_floor, check_frame, FLI, fli_unlogged, MEGNO, RLI, pair_start, pair_distance
export StretchingNumbers, spectral_distance, spectrum_cell, spectrum_density, stretching_range, check_run
export window_weights, window_frequency, WindowFrequency, frequency_map, FrequencyMap, check_frequency_map, circle_distance, frequency_distance, frequency_index
export birkhoff_weight, birkhoff_weights, phase_increments, weighted_phase_frequency, WeightedPhaseAdvance, PlanarObservable, LiftedAngleObservable
export weighted_birkhoff, birkhoff_index, WeightedBirkhoff, CosObservable, SinObservable, LiftedObservable, observable_name

include("sali.jl")
include("alignment.jl")
include("fli.jl")
include("megno.jl")
include("rli.jl")
include("stretching.jl")
include("naff_fft.jl")
include("naff.jl")
include("naff_map.jl")
include("naff_reducer.jl")
include("wbirkhoff.jl")
include("wbirkhoff_reducer.jl")

end
