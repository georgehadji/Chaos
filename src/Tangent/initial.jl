# splitmix64 finalizer; fixed constants make the (seed, orbit) -> stream key stable across Julia versions.
function _mix(z::UInt64)
    z += 0x9e3779b97f4a7c15
    z = (z ⊻ (z >> 30)) * 0xbf58476d1ce4e5b9
    z = (z ⊻ (z >> 27)) * 0x94d049bb133111eb
    return z ⊻ (z >> 31)
end

# Re-orthonormalize the columns in T (two Gram-Schmidt passes), so Double64 and BigFloat reach eps(T).
function _orthonormalize(Q::SMatrix{d, k, T}) where {d, k, T}
    M = Matrix(Q)
    for j in 1:k, _ in 1:2
        for i in 1:(j - 1)
            M[:, j] -= dot(view(M, :, i), view(M, :, j)) * view(M, :, i)
        end
        M[:, j] ./= norm(view(M, :, j))
    end
    return SMatrix{d, k, T}(M)
end

"""
    initial_deviations(seed, orbit, Val(d), Val(k), T) -> SMatrix{d,k,T}

Haar-random orthonormal `d x k` frame for orbit `orbit` of run `seed`: a standard normal matrix from a
`StableRNG` keyed only by `(seed, orbit)`, QR with positive `R` diagonal, first `k` columns (spec, initial vectors).
"""
function initial_deviations(seed::Integer, orbit::Integer, ::Val{d}, ::Val{k}, ::Type{T}) where {d, k, T}
    rng = StableRNG(_mix(_mix(seed % UInt64) + orbit % UInt64))
    g = randn(rng, d * d)
    Q, R = qr(SMatrix{d, d, Float64}(ntuple(i -> g[i], Val(d * d))))
    Qs = Q * Diagonal(sign.(diag(R)))
    Qk = SMatrix{d, k, T}(ntuple(i -> T(Qs[i]), Val(d * k)))
    return T === Float64 ? Qk : _orthonormalize(Qk)
end
