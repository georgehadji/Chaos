# The generic splitting method: a coefficient table (c_i, d_i), i = 1..k, one step = prod_i kick(d_i tau) drift(c_i tau) with the first stage applied
# first (drift with weight c_i, then kick with weight d_i: Yoshida1990 Eq. (2.8); spec Algorithm item 1). The stepper is written once (stepper.jl);
# a method is a cited table. Name and source are type parameters so that a method stays isbits (one kernel for CPU and GPU).

"Supertype of the integrators: the splitting methods, the corrected methods, Gauss and Tao's method."
abstract type AbstractIntegrator end

"""
    SplittingMethod(c, d; order, palindromic, source, name, digits = Inf)

A table method: drift weights `c`, kick weights `d` (k stages, drift then kick, first stage first). `source` is a comma separated list of keys of
`references.bib`, `digits` the digit count of the tabulated constants (`Inf` for algebraic tables). Rejected with `ArgumentError` at construction
for k = 0, a sum condition `|sum c - 1|`, `|sum d - 1|` above `4 k eps_tab`, or a table declared palindromic without `c_i = c_(k+1-i)`,
`d_i = d_(k-i)` (exactly) and `d_k = 0`.
"""
struct SplittingMethod{T, K, N, S} <: AbstractIntegrator
    c::NTuple{K, T}
    d::NTuple{K, T}
    order::Int
    palindromic::Bool
    digits::Float64
end

function SplittingMethod(c::Union{AbstractVector, Tuple}, d::Union{AbstractVector, Tuple}; order::Integer, palindromic::Bool, source, name, digits::Real = Inf)
    k = length(c)
    (k >= 1 && length(d) == k) || throw(ArgumentError("a splitting table needs k >= 1 stages and as many kick as drift weights, got $k and $(length(d))"))
    T = float(promote_type(eltype(c), eltype(d)))
    cv = T[c...]; dv = T[d...]
    _validate_table(cv, dv, palindromic, digits, name)
    return _build_method(Val(k), T, Tuple(cv), Tuple(dv), Int(order), palindromic, Float64(digits), Symbol(name), Symbol(source))
end

_build_method(::Val{K}, ::Type{T}, c, d, order, pal, digits, name, source) where {K, T} =
    SplittingMethod{T, K, name, source}(c, d, order, pal, digits)

function _validate_table(c, d, palindromic, digits, name)
    T = eltype(c); k = length(c)
    (all(isfinite, c) && all(isfinite, d)) || throw(ArgumentError("table $name has a non-finite weight"))
    tol = 4 * k * max(Float64(effective_eps(T)), 10.0^(-Float64(digits)))
    (Float64(abs(sum(c) - one(T))) <= tol && Float64(abs(sum(d) - one(T))) <= tol) ||
        throw(ArgumentError("table $name violates the sum conditions sum c = sum d = 1 (bound $tol)"))
    if palindromic
        (all(i -> c[i] == c[k + 1 - i], 1:k) && all(i -> d[i] == d[k - i], 1:(k - 1)) && iszero(d[k])) ||
            throw(ArgumentError("table $name is declared palindromic but is not (c_i = c_(k+1-i), d_i = d_(k-i), d_k = 0)"))
    end
    return nothing
end

"`coefficients(m) -> (c, d)`: the drift and kick weights of a table method, as tuples."
coefficients(m::SplittingMethod) = (m.c, m.d)
"`stage_count(m)`: the number k of stages of the table."
stage_count(::SplittingMethod{T, K}) where {T, K} = K
"`nominal_order(m)`: the order of the method (for the near-integrable rows the order at eps = 1)."
nominal_order(m::SplittingMethod) = m.order
"`is_palindromic(m)`: whether the table is symmetric, S(tau) S(-tau) = id."
is_palindromic(m::SplittingMethod) = m.palindromic
"`table_digits(m)`: the digit count of the tabulated constants (`Inf` for algebraic tables)."
table_digits(m::SplittingMethod) = m.digits
"`method_name(m)`: the name of the method."
method_name(::SplittingMethod{T, K, N}) where {T, K, N} = String(N)
"`method_source(m)`: the keys of `references.bib` the table comes from, comma separated."
method_source(::SplittingMethod{T, K, N, S}) where {T, K, N, S} = String(S)

"""
    force_evaluations(m) -> Int

Force (kick) evaluations of one step of `m`: the number of maximal runs of kicks with no drift between them (adjacent kicks share one evaluation).
"""
function force_evaluations(m::SplittingMethod{T, K}) where {T, K}
    n = 0; open = false
    for i in 1:K
        iszero(m.c[i]) || (open = false)
        if !iszero(m.d[i])
            open || (n += 1)
            open = true
        end
    end
    return n
end

# Merge the drifts between consecutive kicks, drop empty stages; `trailing` appends the stage (0, 0) when the table ends with a kick (palindromic tables).
function _normalize(cs, ds, trailing::Bool)
    c = eltype(cs)[]; d = eltype(ds)[]; acc = zero(eltype(cs))
    for i in eachindex(cs)
        acc += cs[i]
        if !iszero(ds[i])
            push!(c, acc); push!(d, ds[i]); acc = zero(acc)
        end
    end
    iszero(acc) || (push!(c, acc); push!(d, zero(acc)))
    (trailing && !iszero(d[end])) && (push!(c, zero(acc)); push!(d, zero(acc)))
    return c, d
end

"""
    adjoint_method(m)

The adjoint table (the stages in reverse order): `m` itself for a palindromic table. `inverse(S(tau)) = S*(-tau)`.
"""
function adjoint_method(m::SplittingMethod{T, K, N, S}) where {T, K, N, S}
    m.palindromic && return m
    cs = T[zero(T); (m.c[K + 2 - i] for i in 2:K)...; m.c[1]]
    ds = T[m.d[K]; (m.d[K + 1 - i] for i in 2:K)...; zero(T)]
    c, d = _normalize(cs, ds, false)
    return SplittingMethod(c, d; order = m.order, palindromic = false, source = S, name = string(N, " (adjoint)"), digits = m.digits)
end

# Composition of the pieces `m` with step fractions `gammas` (the first one applied first), neighbouring drifts merged.
function _compose(m::SplittingMethod{T, K, N, S}, gammas, order, name; source = S, digits = m.digits) where {T, K, N, S}
    cs = T[]; ds = T[]
    for g in gammas, i in 1:K
        push!(cs, g * m.c[i]); push!(ds, g * m.d[i])
    end
    c, d = _normalize(cs, ds, m.palindromic)
    return SplittingMethod(c, d; order = order, palindromic = m.palindromic, source = source, name = name, digits = digits)
end

function _check_raise(m, order)
    m.palindromic || throw(ArgumentError("composition raises the order of a symmetric table only"))
    (order > m.order && iseven(order - m.order)) || throw(ArgumentError("target order $order must exceed $(m.order) by an even number"))
    return nothing
end

"""
    triple_jump(m, order)

Raise the symmetric table `m` to `order` (by 2 per level) with the triple jump S(z1 tau) S(z0 tau) S(z1 tau), z0 + 2 z1 = 1,
z0^(p+1) + 2 z1^(p+1) = 0 (Yoshida1990 Eqs. (4.12) to (4.14)); the roots are computed in the element type of `m`.
"""
function triple_jump(m::SplittingMethod{T}, order::Integer) where {T}
    _check_raise(m, order)
    name = string("triple jump ", order, " of ", method_name(m))
    while m.order < order
        z1 = inv(2 - T(2)^inv(T(m.order + 1))); z0 = 1 - 2 * z1
        m = _compose(m, (z1, z0, z1), m.order + 2, name)
    end
    return m
end

"""
    suzuki_jump(m, order)

Raise the symmetric table `m` to `order` with Suzuki's five pieces S(s tau)^2 S((1 - 4 s) tau) S(s tau)^2, `s = 1/(4 - 4^(1/(p+1)))`
(scheme s5odr4 of KahanLi1997 Sect. 3; the condition 4 s^(p+1) + (1 - 4 s)^(p+1) = 0 is the form of Yoshida1990 Eq. (4.13)).
"""
function suzuki_jump(m::SplittingMethod{T}, order::Integer) where {T}
    _check_raise(m, order)
    name = string("Suzuki ", order, " of ", method_name(m))
    while m.order < order
        s = inv(4 - T(4)^inv(T(m.order + 1)))
        m = _compose(m, (s, s, 1 - 4 * s, s, s), m.order + 2, name)
    end
    return m
end

# Palindromic composition of the leapfrog with the step fractions `half` (the last entry the middle one), mirrored.
_palindromic_leapfrog(::Type{T}, half, order, name, source, digits) where {T} =
    _compose(Leapfrog(T, :DKD), vcat(half, reverse(half[1:(end - 1)])), order, name; source = source, digits = digits)
