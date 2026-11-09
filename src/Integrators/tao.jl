# Tao's explicit symplectic method for nonseparable H(q, p) in the extended phase space (catalog row G3, Tao2016): the augmented Hamiltonian on (q, p, x, y)
# H_A + H_B + omega H_C with H_A = H(q, y), H_B = H(x, p), H_C = (|q - x|^2 + |p - y|^2)/2; the three partial flows of Eq. (1), the Strang composition Eq. (2),
# higher order by the triple jump with the exponent 1/(l - 1) (Eq. (3) prints 1/(l + 1), which does not raise the order: spec Inconsistencies (d)).
# State layout (q; p; x; y), each block of length n.

"""
    TaoExtended(T, omega; order = 2)

Tao's method with the binding parameter `omega` (required) and the even `order` >= 2 (2: Strang; higher by the triple jump). Acts on the extended state
(see [`extend_state`](@ref)); symplectic in the extended space, not in (q, p).
"""
struct TaoExtended{T, M} <: AbstractIntegrator
    omega::T
    w::NTuple{M, T}          # step fractions of the Strang pieces
    order::Int
end

function TaoExtended(::Type{T}, omega; order::Integer = 2) where {T}
    (order >= 2 && iseven(order)) || throw(ArgumentError("Tao's method needs an even order >= 2, got $order"))
    w = T[one(T)]
    for l in 4:2:order
        g = inv(2 - T(2)^inv(T(l - 1)))
        w = vcat(g .* w, (1 - 2g) .* w, g .* w)
    end
    return TaoExtended{T, length(w)}(T(omega), Tuple(w), Int(order))
end
TaoExtended(omega::Real; order::Integer = 2) = TaoExtended(Float64, omega; order = order)

_coeff_type(::TaoExtended{T}) where {T} = T
nominal_order(m::TaoExtended) = m.order
is_palindromic(::TaoExtended) = true
adjoint_method(m::TaoExtended) = m
method_name(m::TaoExtended) = string("Tao extended phase space, order ", m.order)
method_source(::TaoExtended) = "Tao2016"
stage_count(::TaoExtended{T, M}) where {T, M} = 5 * M

"`extend_state(z)`: the extended state `(q, p, q, p)` of a physical state `z = (q, p)`."
extend_state(z::SVector) = vcat(z, z)
"`physical_state(zext)`: the physical pair `(q, p)` of an extended state."
physical_state(z::SVector{D}) where {D} = SVector{D ÷ 2}(ntuple(i -> z[i], Val(D ÷ 2)))

@inline _blk(z::SVector, k::Int, ::Val{n}) where {n} = SVector{n}(ntuple(i -> z[(k - 1) * n + i], Val(n)))
@inline _blk(H::SMatrix, a::Int, b::Int, ::Val{n}) where {n} = SMatrix{n, n}(ntuple(l -> H[(a - 1) * n + (l - 1) % n + 1, (b - 1) * n + (l - 1) ÷ n + 1], Val(n * n)))

# flow of H_A = H(q, y): p -= d H_q(q, y), x += d H_p(q, y)
function _tao_A(flow, z::SVector{D}, d) where {D}
    v = Val(D ÷ 4); q = _blk(z, 1, v); p = _blk(z, 2, v); x = _blk(z, 3, v); y = _blk(z, 4, v)
    g = hamiltonian_gradient(flow, vcat(q, y))
    return vcat(q, p - d * _blk(g, 1, v), x + d * _blk(g, 2, v), y)
end
# flow of H_B = H(x, p): q += d H_p(x, p), y -= d H_q(x, p)
function _tao_B(flow, z::SVector{D}, d) where {D}
    v = Val(D ÷ 4); q = _blk(z, 1, v); p = _blk(z, 2, v); x = _blk(z, 3, v); y = _blk(z, 4, v)
    g = hamiltonian_gradient(flow, vcat(x, p))
    return vcat(q + d * _blk(g, 2, v), p, x, y - d * _blk(g, 1, v))
end
# flow of omega H_C: (q - x, p - y) rotates by the angle 2 omega d, the means are fixed
function _tao_C(omega, z::SVector{D}, d) where {D}
    v = Val(D ÷ 4); q = _blk(z, 1, v); p = _blk(z, 2, v); x = _blk(z, 3, v); y = _blk(z, 4, v)
    th = 2 * omega * d; c = _cos(th); s = _sin(th)
    u = q - x; w = p - y; mq = (q + x) / 2; mp = (p + y) / 2
    u2 = c * u + s * w; w2 = c * w - s * u
    return vcat(mq + u2 / 2, mp + w2 / 2, mq - u2 / 2, mp - w2 / 2)
end

@inline function _tao_strang(m::TaoExtended, flow, z, d)
    z = _tao_A(flow, z, d / 2); z = _tao_B(flow, z, d / 2); z = _tao_C(m.omega, z, d)
    z = _tao_B(flow, z, d / 2); return _tao_A(flow, z, d / 2)
end

function advance_n(m::TaoExtended{T, M}, flow, z::SVector, tau, n::Integer; compensated::Bool = false) where {T, M}
    _check_types(m, z, tau)
    compensated && throw(ArgumentError("compensated summation is available for the table methods only"))
    for _ in 1:n, j in 1:M
        z = _tao_strang(m, flow, z, m.w[j] * tau)
    end
    return z
end
advance_to(m::TaoExtended, flow, z::SVector, tau, ns; compensated::Bool = false) =
    (_check_types(m, z, tau); out = Vector{typeof(z)}(undef, length(ns)); done = 0; for (j, nj) in enumerate(ns); z = advance_n(m, flow, z, tau, nj - done); done = nj; out[j] = z; end; out)

# Elementary Jacobians of the three partial flows (blocks of size n)
function _tao_AB_jacobian(flow, z::SVector{D, T}, d, which::Symbol) where {D, T}
    v = Val(D ÷ 4); q = _blk(z, 1, v); p = _blk(z, 2, v); x = _blk(z, 3, v); y = _blk(z, 4, v)
    I = one(SMatrix{D ÷ 4, D ÷ 4, T}); Z = zero(I)
    H = hamiltonian_hessian(flow, which === :A ? vcat(q, y) : vcat(x, p))
    Hqq = _blk(H, 1, 1, v); Hqp = _blk(H, 1, 2, v); Hpq = _blk(H, 2, 1, v); Hpp = _blk(H, 2, 2, v)
    which === :A && return vcat(hcat(I, Z, Z, Z), hcat(-d * Hqq, I, Z, -d * Hqp), hcat(d * Hpq, Z, I, d * Hpp), hcat(Z, Z, Z, I))
    return vcat(hcat(I, d * Hpp, d * Hpq, Z), hcat(Z, I, Z, Z), hcat(Z, Z, I, Z), hcat(Z, -d * Hqp, -d * Hqq, I))
end
function _tao_C_jacobian(omega, ::SVector{D, T}, d) where {D, T}
    n = D ÷ 4; I = one(SMatrix{n, n, T}); th = 2 * omega * d; c = _cos(th); s = _sin(th); a = (1 + c) / 2; b = (1 - c) / 2; h = s / 2
    return vcat(hcat(a * I, h * I, b * I, -h * I), hcat(-h * I, a * I, h * I, b * I), hcat(b * I, -h * I, a * I, h * I), hcat(h * I, b * I, -h * I, a * I))
end

function tangent_advance(m::TaoExtended{T, M}, flow, z::SVector, W::SMatrix, tau) where {T, M}
    _check_types(m, z, tau)
    for j in 1:M
        d = m.w[j] * tau
        W = _tao_AB_jacobian(flow, z, d / 2, :A) * W; z = _tao_A(flow, z, d / 2)
        W = _tao_AB_jacobian(flow, z, d / 2, :B) * W; z = _tao_B(flow, z, d / 2)
        W = _tao_C_jacobian(m.omega, z, d) * W; z = _tao_C(m.omega, z, d)
        W = _tao_AB_jacobian(flow, z, d / 2, :B) * W; z = _tao_B(flow, z, d / 2)
        W = _tao_AB_jacobian(flow, z, d / 2, :A) * W; z = _tao_A(flow, z, d / 2)
    end
    return z, W
end
