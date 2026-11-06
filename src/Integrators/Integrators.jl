"""
    Integrators

Symplectic and reference integrators for Hamiltonian flows (integrator-catalog): the generic splitting stepper over cited coefficient tables, the composition
builders, the near-integrable rows, the Gauss collocation methods, Tao's extended phase-space method, the tangent step by the same splitting, the one-step map for
the Tangent layer and ensembles of orbits. Imports Foundation and Systems only; the Taylor reference arrives through the package extension `ChaoslabTaylorExt`.
"""
module Integrators

using ..Foundation
using ..Foundation: NoDevice, DEVICE      # device hook of the ensemble kernel, not exported
using ..Systems
import ..Systems: dimension, forward, inverse, tangent, inverse_tangent
using ..Systems: _sin, _cos
using StaticArrays
using LinearAlgebra
using KernelAbstractions

export AbstractIntegrator, SplittingMethod, CorrectedMethod, Gauss, ImplicitMidpoint, TaoExtended, IntegratorMap
export advance, advance_n, advance_to, advance_checked, tangent_advance, advance_ensemble
export coefficients, stage_count, nominal_order, is_palindromic, table_digits, force_evaluations, method_name, method_source, adjoint_method
export SymplecticEuler, Leapfrog, McLachlan2, Ruth3, ForestRuth, triple_jump, suzuki_jump, PEFRL, CalvoSanzSerna4, Yoshida, KahanLi, SofroniouSpaletta10
export SABA, SBAB, SABAC, SBABC, ABA864, ABAH864, gauss_coefficients, extend_state, physical_state, taylor_reference

include("splitting.jl")
include("tables.jl")
include("tables_near.jl")
include("stepper.jl")
include("gauss.jl")
include("tao.jl")
include("integrator_map.jl")
include("ensemble.jl")

"""
    taylor_reference(flow, z, t)

The state at time `t` by the Taylor series method in extended precision (reference orbits, catalog row R1). Methods are added by the package extension
`ChaoslabTaylorExt` when `TaylorIntegration` is loaded; the library itself does not depend on it.
"""
function taylor_reference end

end
