# The stepper for the table methods and the corrected methods: one step, n steps with merged kicks, stored states, the tangent step. States are immutable `SVector`s
# (no workspace, so nothing moves to the heap under `--check-bounds=yes`, decision D110); the tables are tuples, the loops run over them.

@inline _qv(z::SVector{D}) where {D} = SVector{D ÷ 2}(ntuple(i -> z[i], Val(D ÷ 2)))
@inline _pv(z::SVector{D}) where {D} = SVector{D ÷ 2}(ntuple(i -> z[D ÷ 2 + i], Val(D ÷ 2)))

# D113 (R2): the coefficient type of the method and a floating point step size must be the element type of the state. Resolved at compile time for matching types.
_coeff_type(::SplittingMethod{T}) where {T} = T
@noinline _type_mismatch(what, T, ET) = throw(ArgumentError("$what of type $T does not match the element type $ET of the state: construct the method and the step size in $ET"))
@inline function _check_types(m::AbstractIntegrator, ::Type{ET}, tau) where {ET}
    ET <: AbstractFloat || return nothing                      # automatic-differentiation numbers: not judged
    _coeff_type(m) === ET || _type_mismatch("the coefficients of the method", _coeff_type(m), ET)
    (tau isa AbstractFloat && !(tau isa ET)) && _type_mismatch("the step size", typeof(tau), ET)
    return nothing
end
@inline _check_types(m::AbstractIntegrator, z::SVector, tau) = _check_types(m, eltype(z), tau)

"""
    advance(m, flow, z, tau)

One step of size `tau` (either sign) of the method `m` on the flow `flow` from the state `z`. Non-finite input gives non-finite output.
"""
advance(m::AbstractIntegrator, flow, z::SVector, tau) = advance_n(m, flow, z, tau, 1)

"""
    advance_n(m, flow, z, tau, n; compensated = false)

`n` steps of size `tau`. Kicks of consecutive stages with no drift between them (also across the step boundary) are applied as one kick, so a kick-drift-kick
leapfrog costs `n + 1` force evaluations. `compensated = true` carries the position and momentum updates as a value and a correction term per component
(Hairer2008 Sect. 2); the corrected methods reject it.
"""
function advance_n(m::SplittingMethod{T, K}, flow, z::SVector, tau, n::Integer; compensated::Bool = false) where {T, K}
    _check_types(m, z, tau)
    compensated && return first(_steps_comp(m, flow, z, zero(z), tau, n))
    pend = zero(T)
    for _ in 1:n
        for i in 1:K
            ci = m.c[i]
            if !iszero(ci)
                iszero(pend) || (z = kick(flow, z, pend * tau); pend = zero(T))
                z = drift(flow, z, ci * tau)
            end
            pend += m.d[i]
        end
    end
    return iszero(pend) ? z : kick(flow, z, pend * tau)
end

# Exact two-term sum: s = a + b and the rounding error of s.
@inline function _twosum(a, b)
    s = a + b; bb = s - a
    return s, (a - (s - bb)) + (b - bb)
end
# z + inc with the correction e fed into the increment (Kahan summation with an exact error term)
@inline function _kahan(z, e, inc)
    s, err = _twosum(z, inc + e)
    return s, err
end

# The same loop with the compensated update: increments from the flow, the sum carried as (value, correction).
function _steps_comp(m::SplittingMethod{T, K}, flow, z::SVector, e::SVector, tau, n::Integer) where {T, K}
    pend = zero(T)
    for _ in 1:n
        for i in 1:K
            ci = m.c[i]
            if !iszero(ci)
                if !iszero(pend)
                    z, e = _kahan(z, e, kick_increment(flow, z, pend * tau)); pend = zero(T)
                end
                z, e = _kahan(z, e, drift_increment(flow, z, ci * tau))
            end
            pend += m.d[i]
        end
    end
    iszero(pend) || ((z, e) = _kahan(z, e, kick_increment(flow, z, pend * tau)))
    return z, e
end

"""
    advance_to(m, flow, z, tau, ns; compensated = false) -> Vector

The states after the step counts `ns` (increasing), from one run.
"""
function advance_to(m::SplittingMethod, flow, z::SVector, tau, ns; compensated::Bool = false)
    _check_types(m, z, tau)
    out = Vector{typeof(z)}(undef, length(ns)); e = zero(z); done = 0
    for (j, nj) in enumerate(ns)
        nj >= done || throw(ArgumentError("the step counts must be increasing"))
        if compensated
            z, e = _steps_comp(m, flow, z, e, tau, nj - done)
        else
            z = advance_n(m, flow, z, tau, nj - done)
        end
        done = nj; out[j] = z
    end
    return out
end

"""
    advance_checked(m, flow, z, tau) -> (znew, status)

One step and its health: `status.ok`, `.iterations` (0 for the explicit methods), `.increment`, `.stages`. The explicit methods are always `ok`.
"""
advance_checked(m::AbstractIntegrator, flow, z::SVector, tau) =
    (advance(m, flow, z, tau), (ok = true, iterations = 0, increment = zero(eltype(z)), stages = ()))

# The elementary tangent operations are not inlined into the stepper (compiled once per flow and element type, not once per table length).
@noinline _tangent_drift(flow, z, W, t) = (drift(flow, z, t), drift_jacobian(flow, z, t) * W)
@noinline _tangent_kick(flow, z, W, t) = (kick(flow, z, t), kick_jacobian(flow, z, t) * W)

"""
    tangent_advance(m, flow, z, W, tau) -> (znew, Wnew)

The deviation frame `W` (`2n x k`, an `SMatrix` or a `Matrix`) advanced through the elementary Jacobians in stage order, with the state: the exact derivative of
the one-step map (spec item 17). No renormalization.
"""
function tangent_advance(m::SplittingMethod{T, K}, flow, z::SVector, W::SMatrix, tau) where {T, K}
    _check_types(m, z, tau)
    for i in 1:K
        iszero(m.c[i]) || ((z, W) = _tangent_drift(flow, z, W, m.c[i] * tau))
        iszero(m.d[i]) || ((z, W) = _tangent_kick(flow, z, W, m.d[i] * tau))
    end
    return z, W
end

function tangent_advance(m::AbstractIntegrator, flow, z::SVector, W::AbstractMatrix, tau)
    zn, Wn = tangent_advance(m, flow, z, SMatrix{length(z), size(W, 2)}(W), tau)
    return zn, W isa SMatrix ? Wn : Matrix(Wn)
end

# ---- corrected methods
@inline function _corrector_kappa(m::CorrectedMethod, flow, tau)
    eps = perturbation_scale(flow)
    return tau * tau * tau * (eps * eps) * m.c / 2
end
@noinline _corrector_kick(flow, z::SVector, kappa) = vcat(_qv(z), _pv(z) + kappa * corrector_force(flow, z))
@noinline function _corrector_jacobian(flow, z::SVector{D, T}, kappa) where {D, T}
    I = one(SMatrix{D ÷ 2, D ÷ 2, T})
    return vcat(hcat(I, zero(I)), hcat(kappa * corrector_hessian(flow, z), I))
end

function advance_n(m::CorrectedMethod, flow, z::SVector, tau, n::Integer; compensated::Bool = false)
    _check_types(m, z, tau)
    compensated && throw(ArgumentError("compensated summation is available for the table methods only"))
    kappa = _corrector_kappa(m, flow, tau)
    for _ in 1:n
        z = _corrector_kick(flow, z, kappa)
        z = advance_n(m.base, flow, z, tau, 1)
        z = _corrector_kick(flow, z, kappa)
    end
    return z
end

function tangent_advance(m::CorrectedMethod, flow, z::SVector, W::SMatrix, tau)
    _check_types(m, z, tau)
    kappa = _corrector_kappa(m, flow, tau)
    W = _corrector_jacobian(flow, z, kappa) * W; z = _corrector_kick(flow, z, kappa)
    z, W = tangent_advance(m.base, flow, z, W, tau)
    W = _corrector_jacobian(flow, z, kappa) * W; z = _corrector_kick(flow, z, kappa)
    return z, W
end

_coeff_type(m::CorrectedMethod) = _coeff_type(m.base)
nominal_order(m::CorrectedMethod) = m.order
is_palindromic(::CorrectedMethod) = true
adjoint_method(m::CorrectedMethod) = m
stage_count(m::CorrectedMethod) = stage_count(m.base)
force_evaluations(m::CorrectedMethod) = isolated_force_evaluations(m.base) + 2       # two corrector gradients per step, no merging across them
method_name(m::CorrectedMethod) = string(method_name(m.base)[1:4], "C", method_name(m.base)[5:end])
method_source(m::CorrectedMethod) = method_source(m.base)
table_digits(m::CorrectedMethod) = table_digits(m.base)
coefficients(m::CorrectedMethod) = coefficients(m.base)
