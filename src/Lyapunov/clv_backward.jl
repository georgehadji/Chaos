# Backward half of the CLV method of Ginelli (catalog L12; spec clv-ginelli V6, V8, V9, V12): the packed storage of the R factors, the
# back-substitution with column normalisation, and the delivery of one sample (C_n, local exponents, cosines, V_n, angles, sink).
# Every per-handoff routine works on preallocated arrays and allocates nothing (spec I6); Float64 and Double64 loops are scalar.

# Packed upper triangle: entry (i, j), i <= j, of the factor of handoff c sits at row j (j - 1) / 2 + i of column c.
@inline _tri(i, j) = (j * (j - 1)) ÷ 2 + i

@inline function _pack!(P::Matrix{T}, c::Int, R::SMatrix{k, k, T}) where {k, T}
    for j in 1:k, i in 1:j
        P[_tri(i, j), c] = R[i, j]
    end
    return nothing
end

@inline _unpack(P::Matrix{T}, c::Int, ::Val{k}) where {T, k} =
    SMatrix{k, k, T}(ntuple(l -> (i = (l - 1) % k + 1; j = (l - 1) ÷ k + 1; i <= j ? P[_tri(i, j), c] : zero(T)), Val(k * k)))

@inline function _putframe!(Gf::Array{T, 3}, c::Int, G::SMatrix{d, k, T}) where {d, k, T}
    for j in 1:k, i in 1:d
        Gf[i, j, c] = G[i, j]
    end
    return nothing
end

@inline _frame(Gf::Array{T, 3}, c::Int, ::Val{d}, ::Val{k}) where {d, k, T} = SMatrix{d, k, T}(ntuple(l -> Gf[(l - 1) % d + 1, (l - 1) ÷ d + 1, c], Val(d * k)))

# Working memory and the outputs of one orbit. `C` is the carried matrix C_{n+1}, `Cn` the next one, `lam` the local exponents of the last step.
mutable struct ClvOut{T, S}
    C::Matrix{T}
    Cn::Matrix{T}
    lam::Vector{T}
    Dt::T                                         # tau h in the element type
    sink::S                                       # nothing or sink(n, C_n, lambda_n)
    sets::Union{Nothing, Tuple{Vector{Int}, Vector{Int}}}
    Cs::Array{T, 3}
    lams::Matrix{T}
    cosines::Array{T, 3}
    V::Array{T, 3}
    G::Array{T, 3}
    R::Array{T, 3}
    angles::Matrix{T}
    wantV::Bool
    wantG::Bool
    wantR::Bool
end

# One backward step V6 + V8 for the interval whose packed factor is column c of P: back-substitution column by column on the top-left
# j x j block (C is upper triangular), overflow-safe norm, unit columns; column 1 is e_1 exactly (spec I1, test plan U13).
# A column with zero or non-finite norm is a "backward breakdown" at the global handoff `step` (P14).
function _bstep!(o::ClvOut{T}, P::Matrix{T}, c::Int, ::Val{k}, step::Int, orbit::Int) where {T, k}
    C, Cn = o.C, o.Cn
    for j in 1:k
        for i in j:-1:1
            s = C[i, j]
            for l in (i + 1):j
                s -= P[_tri(i, l), c] * Cn[l, j]
            end
            Cn[i, j] = s / P[_tri(i, i), c]
        end
        for i in (j + 1):k
            Cn[i, j] = zero(T)
        end
        nr = safe_norm(SVector{k, T}(ntuple(i -> Cn[i, j], Val(k))))     # static vector: no allocation
        (isfinite(nr) && nr > 0) || throw(NumericalHealthError("backward breakdown", Float64(nr), 0.0, step, orbit, :cpu, j))
        for i in 1:j
            Cn[i, j] /= nr
        end
        o.lam[j] = -log(nr) / o.Dt
    end
    Cn[1, 1] = one(T)
    o.C, o.Cn = Cn, C
    return nothing
end

# Principal angles (ascending, in the element type) between the spans of the columns `sets[1]` and `sets[2]` of C_n: V = G C with G orthonormal,
# so the angles of the V sets are those of the C sets (V12 b). Uses `qr` and `svdvals` of the element type; not part of the allocation-free path.
function _angles!(o::ClvOut{T}, n::Int) where {T}
    r1, r2 = o.sets
    Q1 = Matrix(qr(o.C[:, r1]).Q)
    Q2 = Matrix(qr(o.C[:, r2]).Q)
    for (q, s) in enumerate(svdvals(Q1' * Q2))
        o.angles[q, n + 1] = acos(clamp(s, zero(T), one(T)))
    end
    return nothing
end

# Sample n of the window: with a sink only the callback runs; otherwise C_n, lambda_n, the cosines |C^T C| (V12 a), V_n = G_n C_n, G_n, R_n
# (P === nothing for the sample n = n_s when n_b = 0: R and lambda stay NaN) and the angles are written at index n + 1.
function _deliver!(o::ClvOut{T}, n::Int, Gn::SMatrix{d, k, T}, P::Union{Nothing, Matrix{T}}, c::Int) where {T, d, k}
    o.sink === nothing && _store!(o, n, Gn, P, c, Val(k))
    o.sink === nothing || o.sink(n, o.C, o.lam)
    return nothing
end

function _store!(o::ClvOut{T}, n::Int, Gn::SMatrix{d, k, T}, P, c::Int, kv::Val{k}) where {T, d, k}
    C = o.C
    for j in 1:k
        o.lams[j, n + 1] = o.lam[j]
        for i in 1:k
            s = zero(T)
            for l in 1:k
                s += C[l, i] * C[l, j]
            end
            o.cosines[i, j, n + 1] = abs(s)
            o.Cs[i, j, n + 1] = C[i, j]
        end
    end
    if o.wantV
        for j in 1:k, i in 1:d
            s = zero(T)
            for l in 1:j
                s += Gn[i, l] * C[l, j]
            end
            o.V[i, j, n + 1] = s
        end
    end
    if o.wantG
        _putframe!(o.G, n + 1, Gn)
    end
    if o.wantR && P !== nothing
        R = _unpack(P, c, kv)
        for j in 1:k, i in 1:k
            o.R[i, j, n + 1] = R[i, j]
        end
    end
    o.sets === nothing || _angles!(o, n)
    return nothing
end

# Backward over the intervals e - 1, ..., a: column n - off + 1 of P holds the factor of interval n, column n - off + 1 of Gf the frame G_n.
# The sample n is delivered when n <= n_s, after the step that produces C_n; `nf + n + 1` is the global handoff of the interval n.
function _bwd_range!(o::ClvOut{T}, P::Matrix{T}, Gf::Array{T, 3}, off::Int, a::Int, e::Int, ns::Int, nf::Int, kv::Val{k}, dv::Val{d}, frames::Bool, orbit::Int) where {T, k, d}
    for n in (e - 1):-1:a
        _bstep!(o, P, n - off + 1, kv, nf + n + 1, orbit)
        n <= ns && _deliver!(o, n, frames ? _frame(Gf, n - off + 1, dv, kv) : zero(SMatrix{d, k, T}), P, n - off + 1)
    end
    return nothing
end
