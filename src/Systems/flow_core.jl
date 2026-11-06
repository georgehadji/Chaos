# Hamiltonian flows (integrator-catalog, open point P1): the flow side of Systems. A separable flow has H = |p|^2/2 + V(q) and supplies the
# gradient `force` of V and its Hessian `force_jacobian`; everything else (the exact partial flows `drift` and `kick`, their elementary
# Jacobians, the vector field) follows once, here. States are `SVector{2n}` in the layout (q_1..q_n, p_1..p_n).

"Supertype of the Hamiltonian flows: a flow supplies `dimension` and `hamiltonian`."
abstract type AbstractFlow end

"Supertype of the separable flows H = |p|^2/2 + V(q): subtypes define `force` and `force_jacobian` (gradient and Hessian of V), `hamiltonian` and the corrector pair."
abstract type AbstractSeparableFlow <: AbstractFlow end

"`hamiltonian(flow, z)`: the energy H(z) in the element type of `z`."
function hamiltonian end
"`force(flow, q)`: the gradient of the potential V at the positions `q` (separable flows)."
function force end
"`force_jacobian(flow, q)`: the Hessian of V at `q`, an `SMatrix`."
function force_jacobian end
"`hamiltonian_gradient(flow, z)`: `(dH/dq, dH/dp)` as one vector."
function hamiltonian_gradient end
"`hamiltonian_hessian(flow, z)`: the Hessian of H in (q, p) (nonseparable flows)."
function hamiltonian_hessian end
"`vector_field(flow, z)`: J grad H, the right-hand side of the equations of motion."
function vector_field end
"`vector_field_jacobian(flow, z)`: the Jacobian of [`vector_field`](@ref)."
function vector_field_jacobian end
"`drift(flow, z, t)`: the exact flow of T = |p|^2/2 for time `t` (either sign): `q += t p`."
function drift end
"`kick(flow, z, t)`: the exact flow of V for time `t` (either sign): `p -= t grad V(q)`."
function kick end
"`drift_increment(flow, z, t)`: the change `drift(flow, z, t) - z` formed without the sum (compensated summation)."
function drift_increment end
"`kick_increment(flow, z, t)`: the change `kick(flow, z, t) - z` formed without the sum."
function kick_increment end
"`drift_jacobian(flow, z, t)`: the Jacobian of [`drift`](@ref)."
function drift_jacobian end
"`kick_jacobian(flow, z, t)`: the Jacobian of [`kick`](@ref) at the pre-kick point `z`."
function kick_jacobian end
"`exact_flow(flow, z, t)`: the exact solution (harmonic oscillator only)."
function exact_flow end
"`exact_flow_matrix(flow, t)`: the matrix of [`exact_flow`](@ref) (harmonic oscillator only)."
function exact_flow_matrix end
"`angular_momentum(flow, z)`: q x p of the planar Kepler flow."
function angular_momentum end
"`runge_lenz(flow, z)`: the Runge-Lenz vector p x L - q/|q| of the planar Kepler flow."
function runge_lenz end
"`perturbation_scale(flow)`: the eps of H = A + eps B (1 for flows without one)."
function perturbation_scale end
"`corrector_force(flow, z)`: gradient in q of W = |grad V|^2, V the unscaled perturbation (the SABAC corrector)."
function corrector_force end
"`corrector_hessian(flow, z)`: Hessian in q of W = |grad V|^2."
function corrector_hessian end

perturbation_scale(::AbstractFlow) = 1

# The partial flows are not inlined into the steppers: a stepper is compiled for every table length, the flow operations once per flow and element type.

@inline _q(z::SVector{D}) where {D} = SVector{D ÷ 2}(ntuple(i -> z[i], Val(D ÷ 2)))
@inline _p(z::SVector{D}) where {D} = SVector{D ÷ 2}(ntuple(i -> z[D ÷ 2 + i], Val(D ÷ 2)))

@noinline drift(::AbstractSeparableFlow, z::SVector, t) = vcat(_q(z) + t * _p(z), _p(z))
@noinline kick(f::AbstractSeparableFlow, z::SVector, t) = vcat(_q(z), _p(z) - t * force(f, _q(z)))
@noinline drift_increment(::AbstractSeparableFlow, z::SVector, t) = vcat(t * _p(z), zero(_p(z)))
@noinline kick_increment(f::AbstractSeparableFlow, z::SVector, t) = vcat(zero(_q(z)), -(t * force(f, _q(z))))

@noinline function drift_jacobian(::AbstractSeparableFlow, z::SVector{D, T}, t) where {D, T}
    I = one(SMatrix{D ÷ 2, D ÷ 2, promote_type(T, typeof(t))})
    return vcat(hcat(I, t * I), hcat(zero(I), I))
end
@noinline function kick_jacobian(f::AbstractSeparableFlow, z::SVector{D, T}, t) where {D, T}
    Hs = force_jacobian(f, _q(z))
    I = one(SMatrix{D ÷ 2, D ÷ 2, promote_type(T, typeof(t))})
    return vcat(hcat(I, zero(I)), hcat(-t * Hs, I))
end

hamiltonian_gradient(f::AbstractSeparableFlow, z::SVector) = vcat(force(f, _q(z)), _p(z))
vector_field(f::AbstractSeparableFlow, z::SVector) = vcat(_p(z), -force(f, _q(z)))
function vector_field_jacobian(f::AbstractSeparableFlow, z::SVector{D, T}) where {D, T}
    Hs = force_jacobian(f, _q(z))
    I = one(SMatrix{D ÷ 2, D ÷ 2, T})
    return vcat(hcat(zero(I), I), hcat(-Hs, zero(I)))
end

"""
    CountingFlow(flow)

A separable flow that counts its force evaluations (`kick_count`), for the cost tests of the integrators; forwards everything else.
Not isbits (a counter), so it is for single orbits only.
"""
struct CountingFlow{F <: AbstractSeparableFlow} <: AbstractSeparableFlow
    flow::F
    n::Base.RefValue{Int}
end
CountingFlow(flow::AbstractSeparableFlow) = CountingFlow(flow, Ref(0))
"Number of force evaluations of a [`CountingFlow`](@ref) since construction or the last [`reset_counts!`](@ref)."
kick_count(c::CountingFlow) = c.n[]
"Reset the counter of a [`CountingFlow`](@ref)."
reset_counts!(c::CountingFlow) = (c.n[] = 0; c)

dimension(c::CountingFlow) = dimension(c.flow)
hamiltonian(c::CountingFlow, z) = hamiltonian(c.flow, z)
perturbation_scale(c::CountingFlow) = perturbation_scale(c.flow)
force(c::CountingFlow, q) = (c.n[] += 1; force(c.flow, q))
force_jacobian(c::CountingFlow, q) = force_jacobian(c.flow, q)
corrector_force(c::CountingFlow, z) = corrector_force(c.flow, z)
corrector_hessian(c::CountingFlow, z) = corrector_hessian(c.flow, z)
