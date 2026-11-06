# The one-step map of an integrator as a Systems map, so that the Tangent layer propagates it like any map.

_state_dim(::AbstractIntegrator, flow) = dimension(flow)
_state_dim(::TaoExtended, flow) = 2 * dimension(flow)

"""
    IntegratorMap(m, flow, tau)

The one-step map of the method `m` on `flow` with step `tau`, an `AbstractMap`: `forward` is [`advance`](@ref), `inverse` the adjoint table at `-tau`, `tangent` the
exact Jacobian of the step ([`tangent_advance`](@ref) with the identity frame), `inverse_tangent` that of the inverse step at the image point.
"""
struct IntegratorMap{M <: AbstractIntegrator, F, T} <: AbstractMap
    method::M
    flow::F
    tau::T
end

dimension(im::IntegratorMap) = _state_dim(im.method, im.flow)
forward(im::IntegratorMap, z::SVector) = advance(im.method, im.flow, z, im.tau)
inverse(im::IntegratorMap, z::SVector) = advance(adjoint_method(im.method), im.flow, z, -im.tau)
tangent(im::IntegratorMap, z::SVector{D, T}) where {D, T} = last(tangent_advance(im.method, im.flow, z, one(SMatrix{D, D, T}), im.tau))
inverse_tangent(im::IntegratorMap, z::SVector{D, T}) where {D, T} =
    last(tangent_advance(adjoint_method(im.method), im.flow, z, one(SMatrix{D, D, T}), -im.tau))
