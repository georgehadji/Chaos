# Implicit midpoint and the Gauss collocation methods (catalog rows G1, G2), after Hairer2008: Runge-Kutta form Eqs. (2.3), (2.4); fixed-point iteration with the stopping
# rule "iterate until the largest stage increment is 0 or no longer smaller" and the update from f at the stages of the previous iteration (Sect. 3, 'Iteration until
# convergence'); coefficients split into a machine-number part and a small remainder (Eq. (3.1)). Coefficients are computed, not copied: the s-stage method of order 2s.

_tobig(x::BigFloat) = x
_tobig(x::Double64) = BigFloat(x.hi) + BigFloat(x.lo)
_tobig(x::AbstractFloat) = BigFloat(x)

# Gauss-Legendre coefficients in 320-bit arithmetic: nodes c (roots of the shifted Legendre polynomial), weights b, A from the collocation conditions sum_j A_ij c_j^(k-1) = c_i^k/k.
function _gauss_big(s::Integer)
    return setprecision(BigFloat, max(320, precision(BigFloat) + 64)) do
        xs, ws = _legendre_nodes(s, false)
        c = (xs .+ 1) ./ 2; b = ws ./ 2
        V = [c[j]^(k - 1) for k in 1:s, j in 1:s]
        R = [c[i]^k / k for k in 1:s, i in 1:s]
        A = permutedims(V \ R)
        return A, b, c
    end
end

# Round a 320-bit value into T and its remainder (`hi + lo` is the value to about twice the precision of T; lo = 0 for BigFloat).
function _split(::Type{T}, x::BigFloat, prec) where {T}
    T === BigFloat && return BigFloat(x; precision = prec), zero(BigFloat)
    hi = T(x)
    return hi, T(x - _tobig(hi))
end

"""
    gauss_coefficients(T, s) -> (A, b, c)

The coefficients of the s-stage Gauss method (order 2s, symplectic, symmetric) in the element type `T`, computed in extended precision from the Legendre nodes.
"""
function gauss_coefficients(::Type{T}, s::Integer) where {T}
    s >= 1 || throw(ArgumentError("Gauss needs s >= 1 stages"))
    prec = precision(BigFloat)
    A, b, c = _gauss_big(s)
    cv(x) = _split(T, x, prec)[1]
    return SMatrix{s, s, T}(ntuple(l -> cv(A[(l - 1) % s + 1, (l - 1) ÷ s + 1]), Val(s * s))), SVector{s, T}(ntuple(i -> cv(b[i]), Val(s))), SVector{s, T}(ntuple(i -> cv(c[i]), Val(s)))
end

"""
    Gauss(T, s; maxiter = nothing, tol = nothing, standard = false)

The s-stage Gauss method (G2; `ImplicitMidpoint` is s = 1) with the fixed-point solver of Hairer2008 Sect. 3. A step is accepted when the final stage increment is at most
`8 eps max(1, |Y|)` (`tol` replaces this bound by an absolute one); otherwise, after `maxiter` iterations without a stop (default 50, and the precision in bits for `BigFloat`:
a contraction of 0.05 gains about 4 bits per iteration), or for a non-finite value, the step is a failure and returns an all-NaN state. `standard = true` is the implementation of Hairer2008 Fig. 2.1 (stop at an increment below 2 eps, rounded coefficients, update from f at the final
stages), a negative control.
"""
struct Gauss{T, S, L} <: AbstractIntegrator
    Ahi::SMatrix{S, S, T, L}
    Alo::SMatrix{S, S, T, L}
    bhi::SVector{S, T}
    blo::SVector{S, T}
    maxiter::Int
    tol::Float64
    standard::Bool
end

function Gauss(::Type{T}, s::Integer; maxiter::Union{Nothing, Integer} = nothing, tol::Union{Nothing, Real} = nothing, standard::Bool = false) where {T}
    s >= 1 || throw(ArgumentError("Gauss needs s >= 1 stages"))
    maxiter === nothing && (maxiter = T === BigFloat ? max(50, precision(BigFloat)) : 50)
    maxiter >= 2 || throw(ArgumentError("maxiter must be at least 2"))
    prec = precision(BigFloat)
    A, b, _ = _gauss_big(s)
    sp(x) = standard ? (_split(T, x, prec)[1], zero(T)) : _split(T, x, prec)
    Ah = [sp(A[i, j]) for i in 1:s, j in 1:s]; bh = [sp(x) for x in b]
    return Gauss{T, s, s * s}(SMatrix{s, s, T}(first.(Ah)), SMatrix{s, s, T}(last.(Ah)), SVector{s, T}(first.(bh)), SVector{s, T}(last.(bh)),
                              Int(maxiter), tol === nothing ? NaN : Float64(tol), standard)
end
ImplicitMidpoint(::Type{T}) where {T} = Gauss(T, 1)
ImplicitMidpoint() = ImplicitMidpoint(Float64)
Gauss(s::Integer; kw...) = Gauss(Float64, s; kw...)

_coeff_type(::Gauss{T}) where {T} = T
nominal_order(::Gauss{T, S}) where {T, S} = 2S
is_palindromic(::Gauss) = true
adjoint_method(m::Gauss) = m
stage_count(::Gauss{T, S}) where {T, S} = S
method_name(::Gauss{T, S}) where {T, S} = S == 1 ? "implicit midpoint" : string("Gauss ", S)
method_source(::Gauss) = "Hairer2008"
table_digits(::Gauss) = Inf

@inline _inf_norm(v) = maximum(abs, v)
@inline _nan_state(::SVector{D, T}) where {D, T} = SVector{D, T}(ntuple(_ -> T(NaN), Val(D)))

# Error-free transformations (decision D113): the stage sums and the update sum of the Hairer2008 implementation are evaluated with a compensated dot product, so that
# `sum_j (hi_j + lo_j) F_j` is rounded once, as if the products and the additions were done in twice the working precision. The plain sums round every product and every
# partial sum; the rounding errors of the update sum are not independent of those of the stage sums, and over 3 x 10^6 steps they leave a bias in the energy (finding F1
# of the stage 04 report). BigFloat (remainder 0, high precision) and the negative control `standard = true` keep the plain sums.
@inline _prod_err(w, x::SVector) = (p = w * x; (p, fma.(w, x, -p)))
@inline function _dot_comp(wh::SVector{S}, wl::SVector{S}, F::NTuple{S}) where {S}
    s = zero(F[1]); c = zero(F[1])
    for j in 1:S
        for w in (wh[j], wl[j])
            p, pe = _prod_err(w, F[j])
            s, se = _twosum(s, p)
            c = c + (se + pe)
        end
    end
    return s + c
end
@inline _compensated(m::Gauss{T}) where {T} = !m.standard && T !== BigFloat

# sum_j (Ahi_ij + Alo_ij) F_j for every stage i, h included
@inline function _stage_sum(m::Gauss{T, S}, F, h, i) where {T, S}
    _compensated(m) && return h * _dot_comp(m.Ahi[i, :], m.Alo[i, :], F)
    return h * sum(ntuple(j -> m.Ahi[i, j] * F[j], Val(S))) + h * sum(ntuple(j -> m.Alo[i, j] * F[j], Val(S)))
end
@inline function _update_sum(m::Gauss{T, S}, F, h) where {T, S}
    _compensated(m) && return h * _dot_comp(m.bhi, m.blo, F)
    return h * sum(ntuple(j -> m.bhi[j] * F[j], Val(S))) + h * sum(ntuple(j -> m.blo[j] * F[j], Val(S)))
end

# Helpers take the stage tuples as arguments: a closure over a variable that the iteration reassigns would box it (an allocation per step, and not compilable for a kernel).
# The stages are y + Z_i with y = z + e the state carried as value and correction (compensated summation): the correction enters the stage values, not only the update.
@inline _stage_f(flow, Y::NTuple) = map(y -> vector_field(flow, y), Y)
@inline _stage_y(m::Gauss{T, S}, z, e, h, F) where {T, S} = ntuple(i -> z + (e + _stage_sum(m, F, h, i)), Val(S))
@inline _max_diff(A::NTuple{S}, B::NTuple{S}) where {S} = maximum(ntuple(i -> _inf_norm(A[i] - B[i]), Val(S)))
@inline _max_norm(Y::NTuple{S}) where {S} = maximum(ntuple(i -> _inf_norm(Y[i]), Val(S)))

# Fixed-point iteration; returns the update increment, the status and the final stages.
function _gauss_core(m::Gauss{T, S}, flow, z::SVector{D, ET}, h, e::SVector{D, ET} = zero(z)) where {T, S, D, ET}
    f0 = vector_field(flow, z)
    Y = ntuple(i -> z + (h * sum(ntuple(j -> m.Ahi[i, j], Val(S)))) * f0, Val(S))
    F = ntuple(i -> f0, Val(S))
    epsT = ET(effective_eps(ET)); delta = zero(ET); prev = ET(Inf); stopped = false; k = 0
    while k < m.maxiter && !stopped
        k += 1
        F = _stage_f(flow, Y)
        Yn = _stage_y(m, z, e, h, F)
        delta = _max_diff(Yn, Y)
        Y = Yn
        isfinite(delta) || break
        stopped = m.standard ? delta <= 2 * epsT : (iszero(delta) || (k >= 2 && delta >= prev))
        prev = delta
    end
    m.standard && stopped && (F = _stage_f(flow, Y))
    bound = isnan(m.tol) ? 8 * epsT * max(one(ET), _max_norm(Y)) : ET(m.tol)
    ok = stopped && isfinite(delta) && delta <= bound
    return _update_sum(m, F, h), (ok = ok, iterations = k, increment = delta, stages = Y)
end

function advance_checked(m::Gauss, flow, z::SVector, tau)
    _check_types(m, z, tau)
    inc, st = _gauss_core(m, flow, z, tau)
    return st.ok ? z + inc : _nan_state(z), st
end

function advance_n(m::Gauss, flow, z::SVector, tau, n::Integer; compensated::Bool = false)
    _check_types(m, z, tau)
    e = zero(z)
    for _ in 1:n
        inc, st = _gauss_core(m, flow, z, tau, e)
        st.ok || return _nan_state(z)
        z, e = compensated ? _kahan(z, e, inc) : (z + inc, e)
    end
    return z
end

function advance_to(m::Gauss, flow, z::SVector, tau, ns; compensated::Bool = false)
    _check_types(m, z, tau)
    out = Vector{typeof(z)}(undef, length(ns)); e = zero(z); done = 0
    for (j, nj) in enumerate(ns)
        nj >= done || throw(ArgumentError("the step counts must be increasing"))
        for _ in 1:(nj - done)
            inc, st = _gauss_core(m, flow, z, tau, e)
            st.ok || (z = _nan_state(z); break)
            z, e = compensated ? _kahan(z, e, inc) : (z + inc, e)
        end
        done = nj; out[j] = z
    end
    return out
end

# Tangent: differentiate the converged stage equations, dY_i = W + h sum_j a_ij J_j dY_j, then dy' = W + h sum_i b_i J_i dY_i.
function tangent_advance(m::Gauss{T, S}, flow, z::SVector{D, ET}, W::SMatrix{D, K}, tau) where {T, S, D, ET, K}
    _check_types(m, z, tau)
    inc, st = _gauss_core(m, flow, z, tau)
    st.ok || return _nan_state(z), W * ET(NaN)
    Js = ntuple(i -> Matrix(vector_field_jacobian(flow, st.stages[i])), Val(S))
    A = Matrix{ET}(m.Ahi) + Matrix{ET}(m.Alo)
    M = Matrix{ET}(I, S * D, S * D)
    for i in 1:S, j in 1:S
        M[((i - 1) * D + 1):(i * D), ((j - 1) * D + 1):(j * D)] .-= tau * A[i, j] .* Js[j]
    end
    dY = M \ repeat(Matrix(W), S, 1)
    Wn = Matrix(W)
    for i in 1:S
        Wn .+= tau * (m.bhi[i] + m.blo[i]) .* (Js[i] * dY[((i - 1) * D + 1):(i * D), :])
    end
    return z + inc, SMatrix{D, K}(Wn)
end
