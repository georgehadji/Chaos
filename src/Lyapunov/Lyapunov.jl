"""
    Lyapunov

Lyapunov-exponent estimators as reducers of the tangent stream: Benettin (lyap-benettin), the Gram-Schmidt QR spectrum (lyap-gram-schmidt) and their shared helpers, the continuous QR method (continuous-qr) and the covariant Lyapunov vectors of Ginelli's method (clv-ginelli).
"""
module Lyapunov

using ..Foundation
using ..Systems
using ..Integrators: method_name, nominal_order
using ..Tangent: initial_deviations, PropagationResult, _mix
using ..Foundation: fused_product   # frame product with explicit fma (D115), not exported; the same rounding as the tangent stream
using ..Tangent
using StaticArrays
using LinearAlgebra
using StableRNGs

import ..Tangent: reducer_init, reducer_update, reducer_final, reducer_step, reducer_failed, reducer_error, reducer_provenance

export Benettin, Orthonormalizer, CGS, MGS, RGS, orthonormalize, QR_OK, QR_PIVOT, QR_PASSES,
       window_Gw, window_G, orthogonality_threshold, QRSpectrum, HouseholderQR, AugmentedQR, GivensQR, LapackQR, SVDSpectrum, log_singular_values,
       propagate_shadow, propagate_shadow_orbit, separation, displaced,
       ContinuousQR, cqr_rhs, cqr_growth, propagate_cqr_orbit, propagate_cqr_series, propagate_cqr, cqr_bytes, CQ_EMAX,
       ClvGinelli, ClvResult, ClvEnsemble, propagate_clv_orbit, propagate_clv, clv_initial_coefficients, clv_bytes

include("benettin.jl")
include("orthonormalize.jl")
include("qr_spectrum.jl")
include("lapack_qr.jl")
include("householder.jl")
include("svd_spectrum.jl")
include("separation.jl")
include("two_trajectory.jl")
include("continuous_qr.jl")
include("clv_backward.jl")
include("clv_ginelli.jl")

end
