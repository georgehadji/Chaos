"""
    Lyapunov

Lyapunov-exponent estimators as reducers of the tangent stream: Benettin (lyap-benettin), the Gram-Schmidt QR spectrum (lyap-gram-schmidt) and their shared helpers.
"""
module Lyapunov

using ..Foundation
using ..Systems
using ..Tangent
using StaticArrays
using LinearAlgebra

import ..Tangent: reducer_init, reducer_update, reducer_final, reducer_step, reducer_failed, reducer_error, reducer_provenance

export safe_norm, Benettin, Orthonormalizer, CGS, MGS, RGS, orthonormalize, QR_OK, QR_PIVOT, QR_PASSES,
       window_Gw, window_G, orthogonality_threshold, QRSpectrum, HouseholderQR, AugmentedQR, GivensQR, LapackQR, SVDSpectrum, log_singular_values,
       propagate_shadow, propagate_shadow_orbit, separation, displaced

include("norm.jl")
include("benettin.jl")
include("orthonormalize.jl")
include("qr_spectrum.jl")
include("lapack_qr.jl")
include("householder.jl")
include("svd_spectrum.jl")
include("separation.jl")
include("two_trajectory.jl")

end
