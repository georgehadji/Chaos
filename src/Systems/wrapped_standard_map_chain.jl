"""
    WrappedStandardMapChain(N, K, eps; boundary = :open)

[`StandardMapChain`](@ref) on the torus with winding counts (decision D28), state a vector of length `4N`:
`(q_1..q_N, p_1..p_N, count q_1..q_N, count p_1..p_N)`, coordinates in `[0, 2pi_T)` and exact integer counts.
`dimension` is `2N` (the tangent size); `tangent` and `inverse_tangent` read the first `2N` entries. If a count reaches
`maxintfloat(T)` a step returns NaN in all `4N` components. The in-place forms allocate nothing.
"""
struct WrappedStandardMapChain{T} <: AbstractMap
    chain::StandardMapChain{T}
end

function WrappedStandardMapChain(N::Integer, K, eps; boundary::Symbol = :open)
    return WrappedStandardMapChain(StandardMapChain(N, K, eps; boundary))
end

dimension(w::WrappedStandardMapChain) = dimension(w.chain)

# Counts out[c0+1:c0+n] at the exactness limit (or NaN) poison the whole state, so the health check fires.
function _guard_counts!(out, c0, n)
    lim = maxintfloat(eltype(out))
    for i in 1:n
        abs(out[c0 + i]) < lim || return fill!(out, oftype(out[1], NaN))
    end
    return out
end

# W-E2 per site: momenta first (force at the pre-step reduced positions), then positions; m' = m + n' + b, n' = n + a
function forward!(out, w::WrappedStandardMapChain, s)
    m = w.chain
    N = m.N
    tp = 2 * _pi(s[1])
    for j in 1:N
        out[N + j], a = _reduce(s[N + j] + _force(m, s, j), tp)
        out[3N + j] = s[3N + j] + a
    end
    for j in 1:N
        out[j], b = _reduce(s[j] + out[N + j], tp)
        out[2N + j] = s[2N + j] + out[3N + j] + b
    end
    return _guard_counts!(out, 2N, 2N)
end

# W-E4: positions x = x' - p' (b), then p = p' - F(x) (a); m = m' - n' + b, n = n' + a with the image count n'
function inverse!(out, w::WrappedStandardMapChain, s)
    m = w.chain
    N = m.N
    tp = 2 * _pi(s[1])
    for j in 1:N
        out[j], b = _reduce(s[j] - s[N + j], tp)
        out[2N + j] = s[2N + j] - s[3N + j] + b
    end
    for j in 1:N
        out[N + j], a = _reduce(s[N + j] - _force(m, out, j), tp)
        out[3N + j] = s[3N + j] + a
    end
    return _guard_counts!(out, 2N, 2N)
end

forward(w::WrappedStandardMapChain, s::AbstractVector) = forward!(similar(s, _promoted(w.chain, s)), w, s)
inverse(w::WrappedStandardMapChain, s::AbstractVector) = inverse!(similar(s, _promoted(w.chain, s)), w, s)

# W-E3: the tangent depends on the reduced positions only
tangent!(M, w::WrappedStandardMapChain, s) = tangent!(M, w.chain, s)
inverse_tangent!(M, w::WrappedStandardMapChain, s) = inverse_tangent!(M, w.chain, s)
tangent(w::WrappedStandardMapChain, s::AbstractVector) = tangent!(_newmat(w.chain, s), w.chain, s)
inverse_tangent(w::WrappedStandardMapChain, s::AbstractVector) = inverse_tangent!(_newmat(w.chain, s), w.chain, s)

"""
    torus_state(sys::WrappedStandardMapChain, z) -> Vector

State of the wrapped chain for the lifted point `z` (length `2N`): each coordinate reduced into `[0, 2pi_T)` in `T`
with its integer count, laid out `(reduced coordinates, counts)`.
"""
function torus_state(w::WrappedStandardMapChain, z::AbstractVector)
    d = dimension(w)
    tp = 2 * _pi(z[1])
    s = similar(z, 2d)
    for j in 1:d
        s[j], s[d + j] = _reduce(z[j], tp)
    end
    return s
end

"""
    lift(sys::WrappedStandardMapChain, s) -> Vector

Lifted point `z_j + 2pi_T c_j` of the wrapped state `s` (length `4N`), computed in its element type.
"""
function lift(w::WrappedStandardMapChain, s::AbstractVector)
    d = dimension(w)
    tp = 2 * _pi(s[1])
    return [s[j] + tp * s[d + j] for j in 1:d]
end

# Error messages here are constant: string interpolation allocates and does not compile in a GPU kernel (D105d).
# Fixed-size states (SVector of length 4N, amendment D105): the steps as tuples, same expressions and order as the in-place forms.
function _check_wrapped(m, S)
    (S % 4 == 0 && m.N == S ÷ 4) || throw(DimensionMismatch("wrapped state length is not 4N for this chain"))
    return S ÷ 4
end

function _guard_counts(x::SVector{S, T}, N) where {S, T}
    lim = maxintfloat(T)
    for i in (2N + 1):(4N)
        abs(x[i]) < lim || return SVector{S, T}(ntuple(_ -> T(NaN), Val(S)))
    end
    return x
end

# Reduction of the entries v (an SVector) into [0, tp): the reduced values and the counts.
@inline function _reduced(v::SVector{N}, tp) where {N}
    R = ntuple(j -> _reduce(v[j], tp), Val(N))
    return SVector{N}(ntuple(j -> R[j][1], Val(N))), SVector{N}(ntuple(j -> R[j][2], Val(N)))
end

function forward(w::WrappedStandardMapChain, s::SVector{S}) where {S}
    m, N = w.chain, _check_wrapped(w.chain, S)
    T, V, tp = _promoted(m, s), Val(N), 2 * _pi(s[1])
    pp, a = _reduced(SVector{N, T}(ntuple(j -> s[N + j] + _force(m, s, j), V)) , tp)
    nn = _blk(s, 3N, V) + a
    qq, b = _reduced(_blk(s, 0, V) + pp, tp)
    return _guard_counts(vcat(qq, pp, _blk(s, 2N, V) + nn + b, nn), N)
end

function inverse(w::WrappedStandardMapChain, s::SVector{S}) where {S}
    m, N = w.chain, _check_wrapped(w.chain, S)
    T, V, tp = _promoted(m, s), Val(N), 2 * _pi(s[1])
    q, b = _reduced(_blk(s, 0, V) - _blk(s, N, V), tp)
    mm = _blk(s, 2N, V) - _blk(s, 3N, V) + b
    pp, a = _reduced(SVector{N, T}(ntuple(j -> s[N + j] - _force(m, q, j), V)), tp)
    return _guard_counts(vcat(q, pp, mm, _blk(s, 3N, V) + a), N)
end

function tangent(w::WrappedStandardMapChain, s::SVector{S}) where {S}
    return _chain_tangent(w.chain, s, Val(_check_wrapped(w.chain, S)), _promoted(w.chain, s), false)
end
function inverse_tangent(w::WrappedStandardMapChain, s::SVector{S}) where {S}
    N = _check_wrapped(w.chain, S)
    return _chain_tangent(w.chain, _Preimage(s, N), Val(N), _promoted(w.chain, s), true)
end

function torus_state(w::WrappedStandardMapChain, z::SVector{d, T}) where {d, T}
    d == dimension(w) || throw(DimensionMismatch("lifted state length is not 2N for this chain"))
    tp = 2 * _pi(z[1])
    return vcat(_reduced(z, tp)...)
end

function lift(w::WrappedStandardMapChain, s::SVector{S}) where {S}
    d = S ÷ 2
    (iseven(S) && d == dimension(w)) || throw(DimensionMismatch("wrapped state length is not 4N for this chain"))
    tp = 2 * _pi(s[1])
    return SVector{d}(ntuple(j -> s[j] + tp * s[d + j], Val(d)))
end
