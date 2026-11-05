# GALI reducer (spec gali, Algorithm 2 to 4, B1 to B5, B12): K deviation vectors, each renormalised at every handoff; GALI_k for
# every order 2..K from the nested factorisation of the unit matrix; the log form at K = d.

"""
    GALI()

Reducer for `propagate` and `propagate_orbit` with `k = K` columns, `2 <= K <= d` (`reducer_init` throws `ArgumentError` outside
that range). At each handoff, for every column `i`, adds `log(alpha_i)`, `alpha_i` the scaled norm of column `i`, to a compensated
sum `S_i`, replaces `W` by its column-normalised matrix and evaluates `G_k = GALI_k`, `k = 2..K`, from the stored unit columns
([`gali_values`](@ref)). The reducer replaces the columns, so it is the last of its chain (D71). `reducer_final` gives
`(n, gali, L, sum, comp, floor)`: `gali` the `SVector{K - 1}` `(G_2, ..., G_K)` at the last handoff `n`; `L = -(S_1 + ... + S_d)` the
log form `ln GALI_d` of an orthonormal initial frame when `K = d` (for `d = 2` bitwise the [`SALI`](@ref) `L`), `nothing` when
`K < d`; `sum` and `comp` the `K` sums and their compensation terms; `floor` the measurement floors [`gali_floor`](@ref) of the
orders `2..K`. A zero or non-finite column norm is recorded as a health failure ("deviation matrix", the first failing column) and
leaves the sums and `W` unchanged. Isbits state, no allocation after warm-up, no library factorisation (the stream kernel compiles).
"""
struct GALI end

# `code` 0 = ok, 1 = zero or non-finite norm `value` of column `column` at handoff `step`; `n` is the last successful handoff.
struct GALIState{T, D, K, M}
    s::SVector{K, T}
    c::SVector{K, T}
    gali::SVector{M, T}
    value::T
    n::Int
    step::Int
    code::Int8
    column::Int
end

"""
    gali_floor(T, k) -> T

Measurement floor `B_meas(T, k) = 2 k effective_eps(T)` of the direct form of `GALI_k` (spec B9): a computed value below it is
unresolved, not wrong.
"""
gali_floor(::Type{T}, k::Integer) where {T} = 2 * k * effective_eps(T)

function reducer_init(::GALI, ::Type{T}, ::Val{D}, ::Val{K}) where {T, D, K}
    (2 <= K <= D) || throw(ArgumentError("GALI needs 2 <= k <= d"))  # constant message: the kernel can reach this
    return GALIState{T, D, K, K - 1}(zeros(SVector{K, T}), zeros(SVector{K, T}), zeros(SVector{K - 1, T}), zero(T), 0, 0, Int8(0), 0)
end

# Scaled Euclidean norm of column j of a matrix, the arithmetic of `safe_norm` as plain loops (a fixed-size loop, not an unrolled
# expression: the unrolled form at d = 200 takes minutes to compile for Double64).
@inline Base.@propagate_inbounds function _norm_col(M, j::Int, ::Val{D}) where {D}
    m = abs(M[1, j])
    for i in 2:D
        m = max(m, abs(M[i, j]))
    end
    (isfinite(m) && m > 0) || return m
    s = abs2(M[1, j] / m)
    for i in 2:D
        s += abs2(M[i, j] / m)
    end
    return m * sqrt(s)
end

# Work matrix: a stack buffer for isbits element types, a heap matrix otherwise (BigFloat, CPU only).
@inline _buffer(::Type{T}, ::Val{D}, ::Val{K}) where {T, D, K} = isbitstype(T) ? MMatrix{D, K, T}(undef) : Matrix{T}(undef, D, K)

# Column-normalised copy of W (division, as the sali reducer).
@inline function _unit_buffer(W::SMatrix{D, K, T}, a) where {D, K, T}
    U = _buffer(T, Val(D), Val(K))
    @inbounds for j in 1:K, i in 1:D
        U[i, j] = W[i, j] / a[j]
    end
    return U
end

# Nested QR of the unit columns in place: `g_k = prod_{i = 2..k} abs(R_ii)` is GALI_k (spec B2 iv). Modified Gram-Schmidt, each
# column orthogonalised twice against the earlier ones, so the first k values do not depend on K (nesting is bitwise). Destroys `Q`.
# `@inbounds` is required: a bounds-error path that holds `Q` defeats the stack allocation (zero allocation per step at d = 200).
@inline function _gram_inplace!(Q::AbstractMatrix{T}, ::Val{D}, ::Val{K}) where {T, D, K}
    G, g = zeros(SVector{K, T}), one(T)
    @inbounds for j in 1:K
        for _ in 1:2, p in 1:(j - 1)
            c = zero(T)
            for i in 1:D
                c += Q[i, p] * Q[i, j]
            end
            for i in 1:D
                Q[i, j] -= c * Q[i, p]
            end
        end
        r = _norm_col(Q, j, Val(D))
        s = iszero(r) ? one(T) : r    # r = 0 means a zero column: it stays zero (a ternary inside the loop defeats stack allocation)
        for i in 1:D
            Q[i, j] /= s
        end
        j > 1 && (g *= r; G = Base.setindex(G, g, j))
    end
    return _tail(G)
end
@inline _tail(G::SVector{K, T}) where {K, T} = SVector{K - 1, T}(ntuple(i -> G[i + 1], Val(K - 1)))

"""
    gali_values(W) -> SVector{K - 1}

`(G_2, ..., G_K)` of the `d x K` matrix `W` (`SMatrix` or `Matrix`, any element type, BigFloat included): each column is divided by
its scaled norm, then a nested QR (modified Gram-Schmidt, each column orthogonalised twice) gives `G_k = prod_{i <= k} abs(R_ii)` for
every `k` at once (spec B2 iv): the volume of the parallelepiped of the first `k` unit columns, `0` for dependent columns, within
[`gali_floor`](@ref) of the true value. A non-finite column gives non-finite values. Written for fixed sizes, no library call.
"""
function gali_values(W::SMatrix{D, K}) where {D, K}
    K >= 2 || throw(ArgumentError("GALI needs a matrix with at least 2 columns, got $K"))
    return _gram_inplace!(_unit_buffer(W, ntuple(j -> _norm_col(W, j, Val(D)), Val(K))), Val(D), Val(K))
end
gali_values(W::AbstractMatrix) = gali_values(SMatrix{size(W, 1), size(W, 2)}(W))

"""
    gali_value(W, k) -> T

`G_k` of the first `k` columns of `W`; `ArgumentError` for `k` outside `2..size(W, 2)`.
"""
function gali_value(W::AbstractMatrix, k::Integer)
    2 <= k <= size(W, 2) || throw(ArgumentError("GALI order k must lie in 2..$(size(W, 2)), got $k"))
    return gali_values(Matrix(W)[:, 1:k])[k - 1]
end

# First column with a zero or non-finite norm, 0 if none.
@inline function _first_bad(a::NTuple{K}) where {K}
    for i in 1:K
        (isfinite(a[i]) && a[i] > 0) || return i
    end
    return 0
end

function reducer_step(::GALI, st::GALIState{T, D, K}, n::Int, x::SVector, W::SMatrix{D, K}) where {T, D, K}
    st.code == 0 || return st, W
    a = ntuple(j -> _norm_col(W, j, Val(D)), Val(K))
    bad = _first_bad(a)
    if bad != 0
        return GALIState{T, D, K, K - 1}(st.s, st.c, st.gali, T(a[bad]), st.n, n, Int8(1), bad), W
    end
    accs = ntuple(i -> add(CompensatedSum{T}(st.s[i], st.c[i]), log(a[i])), Val(K))
    Q = _unit_buffer(W, a)
    U = SMatrix{D, K, T}(Q)
    st1 = GALIState{T, D, K, K - 1}(SVector(map(q -> q.s, accs)), SVector(map(q -> q.c, accs)), _gram_inplace!(Q, Val(D), Val(K)), st.value, n, st.step, st.code, st.column)
    return st1, U
end

# Direct use throws at once (orbit unknown, 0); inside `propagate` the body uses `reducer_step`.
function reducer_update(r::GALI, st::GALIState, n::Int, x::SVector, W::SMatrix)
    st1, W1 = reducer_step(r, st, n, x, W)
    reducer_failed(r, st1) && throw(reducer_error(r, st1, 0, :cpu))
    return st1, W1
end

reducer_failed(::GALI, st::GALIState) = st.code != 0

function reducer_error(::GALI, st::GALIState, orbit, backend)
    value = Float64(st.value)
    return NumericalHealthError("deviation matrix", value, isfinite(value) ? 0.0 : Inf, st.step, orbit, backend, st.column)
end

# Left-to-right sums of the K sums and of the K compensation terms, then the negative of their total (d = 2: the sali form).
function _log_form(st::GALIState{T, D, K}) where {T, D, K}
    s, c = st.s[1], st.c[1]
    for i in 2:K
        s += st.s[i]
        c += st.c[i]
    end
    return -(s + c)
end

function reducer_final(::GALI, st::GALIState{T, D, K}) where {T, D, K}
    L = K == D ? (st.n == 0 ? zero(T) : _log_form(st)) : nothing
    floors = SVector{K - 1, T}(ntuple(i -> gali_floor(T, i + 1), Val(K - 1)))
    return (n = st.n, gali = st.gali, L = L, sum = st.s, comp = st.c, floor = floors)
end

reducer_provenance(::GALI) = (indicator = "gali",)

"""
    check_frame(GALI(), W0; orbit = 0) -> G_K(0)

Frame check of a supplied initial deviation frame (spec Algorithm 1, Amendment D107e): returns `G_K(0)` of the `K` columns if it is at
least `10 gali_floor(T, K)`, else throws `NumericalHealthError("initial frame", ...)` at step 0 (dependent columns would give an
absolute zero and read as chaos). Because the volumes are nested the check covers every lower order. The caller runs it on every
supplied frame before the stream starts; `propagate` does not.
"""
function check_frame(::GALI, W0::AbstractMatrix{T}; orbit::Integer = 0) where {T}
    size(W0, 2) >= 2 || throw(ArgumentError("GALI needs a frame with at least 2 columns, got $(size(W0, 2))"))
    g0 = gali_values(W0)[end]
    threshold = 10 * gali_floor(T, size(W0, 2))
    g0 >= threshold || throw(NumericalHealthError("initial frame", Float64(g0), Float64(threshold), 0, orbit, :cpu, 0))
    return g0
end
