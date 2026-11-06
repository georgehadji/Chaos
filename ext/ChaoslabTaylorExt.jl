"""
    ChaoslabTaylorExt

The Taylor series method as the reference integrator R1 of the integrator catalog (JorbaZou2005 Secs. 2 and 3): adds the methods of `Chaoslab.Integrators.taylor_reference`
through `TaylorIntegration`, loaded in the verification environment only; the library does not depend on it.
"""
module ChaoslabTaylorExt

using Chaoslab, TaylorIntegration, StaticArrays
using Chaoslab.Integrators: Integrators

# The library's `sin` and `cos` of the flows guard non-finite input with `isfinite`; for a Taylor polynomial they are the plain functions.
Chaoslab.Systems._sin(x::TaylorIntegration.Taylor1) = sin(x)
Chaoslab.Systems._cos(x::TaylorIntegration.Taylor1) = cos(x)

"""
    taylor_reference(flow, z, t; order = 40, abstol = 2^-200)

The state at time `t` of the flow from `z` (an `SVector` of `BigFloat`s): Taylor series method of order `order` with the step chosen from the last coefficients so that
the truncation error of a step is at most `abstol`.
"""
function Integrators.taylor_reference(flow::Chaoslab.AbstractFlow, z::SVector{D, T}, t; order::Integer = 40, abstol = ldexp(T(1), -200)) where {D, T <: AbstractFloat}
    f!(dx, x, p, s) = (v = Chaoslab.vector_field(flow, SVector{D}(ntuple(i -> x[i], Val(D)))); for i in 1:D; dx[i] = v[i]; end; nothing)
    sol = taylorinteg(f!, collect(z), zero(T), T(t), Int(order), T(abstol); maxsteps = 10^7, parse_eqs = false)
    return SVector{D, T}(ntuple(i -> sol.x[end, i], Val(D)))
end

end
