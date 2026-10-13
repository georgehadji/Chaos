"""
    Diagnostics

Residual reducers: det, symplecticity, orthogonality, sum and pairing (diagnostics).
"""
module Diagnostics

using ..Foundation
using ..Systems
using ..Tangent
using StaticArrays
using LinearAlgebra

import ..Tangent: reducer_init, reducer_update, reducer_final

export eps_det, eps_symp, eps_orth, eps_sum, eps_pair
export det_threshold, symp_threshold, sum_threshold, check_residual
export ResidualMonitor, monitor_exponents

include("residuals.jl")
include("monitor.jl")

end
