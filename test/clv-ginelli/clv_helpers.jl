# Helpers for the clv-ginelli suite (prefix clv_). Definitions only: nothing here calls a new Chaoslab name (ClvGinelli, propagate_clv_orbit, propagate_clv, clv_bytes,
# clv_initial_coefficients) at load time; every call sits inside a function body or a testset, so a missing module errors the testsets that use it and nothing else.
# Spec: workspaces/module-build/stages/01-spec/output/clv-ginelli-spec.md (D122); API and conventions: workspaces/module-build/stages/02-tests/output/clv-ginelli-test-plan.md.
using Test, Chaoslab, StaticArrays, LinearAlgebra, StableRNGs, DoubleFloats

const CLV_SEED = 20261012                 # _config/numerics-defaults.toml, key `seed`
const CLV_BASE = CLV_SEED + 45000         # this slug uses 20261012 + 45000 .. 45899; testset AC<n> uses 50 n + sub (sub < 50); 45900+ reserved for verify/
clv_seed(ac::Integer, sub::Integer = 0) = CLV_BASE + 50 * ac + sub
const CLV_BITS = 256                      # numerics-defaults.toml, key `bigfloat_bits`
# Package sizes (D38, 3 minute ceiling for the folder) or the full sizes of the spec: verify/clv-ginelli.jl sets CLV_FULL_RUN = true before including runtests.jl.
const CLV_FULL = isdefined(Main, :CLV_FULL_RUN) && Main.CLV_FULL_RUN
clv_size(small, full) = CLV_FULL ? full : small

clv_eps(::Type{T}) where {T} = Float64(effective_eps(T))
clv_big(x::Double64) = BigFloat(x.hi) + BigFloat(x.lo)     # exact
clv_big(x::AbstractFloat) = BigFloat(x)
clv_big(x::Integer) = BigFloat(x)
clv_conv(::Type{T}, x::BigFloat) where {T} = T(x)
clv_conv(::Type{Double64}, x::BigFloat) = (hi = Float64(x); Double64(hi, Float64(x - hi)))
clv_conv(::Type{T}, x) where {T} = T(x)
const CLV_CACHE = Dict{Any, Any}()

# ------------------------------------------------------------------ angles and distances between directions, computed in 256 bits from the stored values
# delta(u, w): the sine of the angle between two directions (spec acceptance criteria preamble), as the norm of the component of u/|u| perpendicular to w/|w| (no cancellation).
function clv_delta(u::AbstractVector, w::AbstractVector)
    ub = map(clv_big, u); wb = map(clv_big, w); ub = ub / norm(ub); wb = wb / norm(wb)
    return Float64(norm(ub - dot(ub, wb) * wb))
end
# Pairwise angle (atan2 form, accurate for tiny angles).
function clv_pair_angle(u::AbstractVector, w::AbstractVector)
    ub = map(clv_big, u); wb = map(clv_big, w); ub = ub / norm(ub); wb = wb / norm(wb)
    c = abs(dot(ub, wb)); s = norm(ub - dot(ub, wb) * wb)
    return Float64(atan(s, c))
end
# Sine of the largest principal angle of span(U) from span(W): spectral norm of (I - P_W) Q_U (the projection formed in 256 bits, the norm in Float64).
function clv_sin_sub(U::AbstractMatrix, W::AbstractMatrix)
    Ub = Matrix{BigFloat}(map(clv_big, U)); Wb = Matrix{BigFloat}(map(clv_big, W))
    QU = Matrix(qr(Ub).Q); QW = Matrix(qr(Wb).Q)
    return opnorm(Float64.(QU - QW * (QW' * QU)))
end
clv_fl(x) = Float64.(x)

# ------------------------------------------------------------------ test fixtures (never in src/): a constant linear map, and wrappers of maps used as doubles
struct clv_LinMap{D, T, L} <: Chaoslab.AbstractMap
    A::SMatrix{D, D, T, L}
    Ainv::SMatrix{D, D, T, L}
end
clv_LinMap(A::SMatrix{D, D, T, L}) where {D, T, L} = clv_LinMap{D, T, L}(A, inv(A))
Chaoslab.dimension(::clv_LinMap{D}) where {D} = D
Chaoslab.forward(m::clv_LinMap{D, T}, z::SVector{D}) where {D, T} = m.A * z
Chaoslab.inverse(m::clv_LinMap{D, T}, z::SVector{D}) where {D, T} = m.Ainv * z
Chaoslab.tangent(m::clv_LinMap, ::SVector) = m.A
Chaoslab.inverse_tangent(m::clv_LinMap, ::SVector) = m.Ainv
# The inverse of an invertible map as a map (forward = inverse, tangent = inverse tangent at the image point), for the time-reversal criterion.
struct clv_Rev{M} <: Chaoslab.AbstractMap
    map::M
end
Chaoslab.dimension(r::clv_Rev) = dimension(r.map)
Chaoslab.forward(r::clv_Rev, z::SVector) = inverse(r.map, z)
Chaoslab.tangent(r::clv_Rev, z::SVector) = inverse_tangent(r.map, z)
# Counting wrapper: calls of `forward` and of `tangent` (single orbits only).
struct clv_CountMap{M} <: Chaoslab.AbstractMap
    map::M
    nf::Base.RefValue{Int}
    nt::Base.RefValue{Int}
end
clv_CountMap(m) = clv_CountMap(m, Ref(0), Ref(0))
Chaoslab.dimension(c::clv_CountMap) = dimension(c.map)
Chaoslab.forward(c::clv_CountMap, z::SVector) = (c.nf[] += 1; forward(c.map, z))
Chaoslab.tangent(c::clv_CountMap, z::SVector) = (c.nt[] += 1; tangent(c.map, z))
# A map that returns NaN from its `from`-th forward call on (fault injection for criterion 14 (a)).
struct clv_NaNMap{M} <: Chaoslab.AbstractMap
    map::M
    from::Int
    n::Base.RefValue{Int}
end
clv_NaNMap(m, from) = clv_NaNMap(m, from, Ref(0))
Chaoslab.dimension(c::clv_NaNMap) = dimension(c.map)
function Chaoslab.forward(c::clv_NaNMap, z::SVector)
    c.n[] += 1
    w = forward(c.map, z)
    return c.n[] >= c.from ? w * oftype(first(w), NaN) : w
end
Chaoslab.tangent(c::clv_NaNMap, z::SVector) = tangent(c.map, z)
# Separable linear flow H = |p|^2/2 - q^T K q / 2 (criterion 1 (d)): force = grad V = -K q, so the vector field is (p, K q) with the matrix A = [[0, I], [K, 0]].
struct clv_QuadFlow{N, T, L} <: Chaoslab.AbstractSeparableFlow
    K::SMatrix{N, N, T, L}
end
Chaoslab.dimension(::clv_QuadFlow{N}) where {N} = 2N
Chaoslab.force(f::clv_QuadFlow, q) = -(f.K * q)
Chaoslab.force_jacobian(f::clv_QuadFlow, q) = -f.K
Chaoslab.hamiltonian(f::clv_QuadFlow, z::SVector{D}) where {D} = (N = D ÷ 2; q = SVector{N}(ntuple(i -> z[i], Val(N))); p = SVector{N}(ntuple(i -> z[N + i], Val(N))); dot(p, p) / 2 - dot(q, f.K * q) / 2)

# ------------------------------------------------------------------ 256-bit exact fixtures (spec V16): exact CLV matrix (unit columns), exponents per unit time, gaps
function clv_expm(A::AbstractMatrix{BigFloat})
    nA = maximum(sum(abs, A; dims = 2)); s = max(0, ceil(Int, log2(Float64(nA) + 1.0e-300)) + 8)
    As = A / BigFloat(2)^s; E = Matrix{BigFloat}(I, size(A)...); term = copy(E)
    for j in 1:60
        term = term * As / j; E += term
    end
    for _ in 1:s
        E = E * E
    end
    return E
end
clv_unit(M::AbstractMatrix) = M ./ sqrt.(sum(abs2, M; dims = 1))
const CLV_B = [1.0 2.0; 2.0 1.0]
# Each fixture: map (element type T), exact V (BigFloat d x d, columns in decreasing-exponent order), exact exponents per unit time (BigFloat), the handoff interval tau, the step h,
# and gamma = g_min tau h (V22). :a is KuptsovParlitz2012 Eq. (69) as a map (M = exp(J)); :b the non-normal symplectic map; :c the Arnold cat map; :d the separable flow (order-4 table).
function clv_fixture(name::Symbol, ::Type{T}; tau::Integer = name === :d ? 20 : 1) where {T}
    return setprecision(BigFloat, CLV_BITS) do
        if name === :a
            J = BigFloat[1 -2 0; 0 -1 0; 0 2 -3]
            M = clv_expm(J); V = clv_unit(BigFloat[1 1 0; 0 1 0; 0 1 1]); lam = BigFloat[1, -1, -3]
            map = clv_LinMap(SMatrix{3, 3, T}(clv_conv.(T, M)))
        elseif name === :b
            a1, a2 = BigFloat(5) / 2, BigFloat(8) / 5; B = BigFloat.(CLV_B)
            S = [Matrix{BigFloat}(I, 2, 2) B; zeros(BigFloat, 2, 2) Matrix{BigFloat}(I, 2, 2)]
            M = S * Diagonal([a1, a2, 1 / a1, 1 / a2]) * inv(S)
            V = clv_unit(S[:, [1, 2, 4, 3]]); lam = [log(a1), log(a2), -log(a2), -log(a1)]
            map = clv_LinMap(SMatrix{4, 4, T}(clv_conv.(T, M)))
        elseif name === :c
            mu = (3 + sqrt(BigFloat(5))) / 2
            M = BigFloat[2 1; 1 1]; V = clv_unit(BigFloat[1 1; mu - 2 1 / mu - 2]); lam = [log(mu), -log(mu)]
            map = clv_LinMap(SMatrix{2, 2, T}(clv_conv.(T, M)))
        elseif name === :d
            c, s = cos(BigFloat(6) / 10), sin(BigFloat(6) / 10); Rot = [c -s; s c]
            k1, k2 = BigFloat(64) / 100, BigFloat(16) / 100; Kb = Rot * Diagonal([k1, k2]) * Rot'
            e1, e2 = Rot[:, 1], Rot[:, 2]; w1, w2 = sqrt(k1), sqrt(k2)
            V = clv_unit([vcat(e1, w1 * e1) vcat(e2, w2 * e2) vcat(e2, -w2 * e2) vcat(e1, -w1 * e1)])
            lam = [w1, w2, -w2, -w1]
            map = IntegratorMap(ForestRuth(T), clv_QuadFlow(SMatrix{2, 2, T}(clv_conv.(T, Kb))), T(1) / 20)
        else
            error("unknown fixture $name")
        end
        h = name === :d ? 1 / 20 : 1.0
        gaps = [Float64(lam[i] - lam[i + 1]) for i in 1:(length(lam) - 1)]
        return (map = map, V = V, lam = lam, tau = Int(tau), h = h, Dt = tau * h, gamma = minimum(gaps) * tau * h, d = length(lam), name = name)
    end
end
clv_B_fix(::Type{T}, gamma, n) where {T} = 1.0e3 * clv_eps(T) + 100 * exp(-gamma * n)       # spec preamble: 1e3 eps + 4 exp(-gamma n), proposal; prefactor 100 calibrated (plan U4)
clv_B_cov(::Type{T}, kappa) where {T} = 64 * clv_eps(T) * max(1.0, kappa)                   # spec preamble, proposal
# Planted tangency d = 2 (V16 (e)): M = Gamma diag(2, 1/2) Gamma^-1, Gamma = [[1, 1], [0, delta]].
function clv_planted(delta, ::Type{T}) where {T}
    return setprecision(BigFloat, CLV_BITS) do
        d = BigFloat(delta); M = [BigFloat(2) -BigFloat(3) / (2 * d); BigFloat(0) BigFloat(1) / 2]
        clv_LinMap(SMatrix{2, 2, T}(clv_conv.(T, M)))
    end
end

# ------------------------------------------------------------------ frames and starting coefficient matrices (test-side)
clv_eye(d::Integer, k::Integer, ::Type{T}) where {T} = SMatrix{d, k, T}(ntuple(l -> (l - 1) % d == (l - 1) ÷ d ? one(T) : zero(T), Val(d * k)))
clv_frame(seed::Integer, orbit::Integer, d::Integer, k::Integer, ::Type{T}) where {T} = initial_deviations(seed, orbit, Val(d), Val(k), T)
# Upper-triangular k x k matrix, positive diagonal, unit columns, from a recorded StableRNGs stream (the test-side C0 handed to the module and to the reference).
function clv_rand_C(rng, k::Integer, ::Type{T}) where {T}
    M = zeros(Float64, k, k)
    for j in 1:k
        v = randn(rng, j); v[j] = abs(v[j]) + 0.1; M[1:j, j] = v / norm(v)
    end
    return SMatrix{k, k, T}(T.(M))
end
clv_Cstart(rng, k, T) = clv_rand_C(rng, k, T)

# ------------------------------------------------------------------ chaotic starts and map objects of the panel
clv_std(::Type{T}, K = 5) where {T} = StandardMap(exact(T, K))
clv_cm(::Type{T}) where {T} = StandardMapChain(2, exact(T, 3), exact(T, 1) / 5)        # Y4, N = 2, K = 3, eps = 0.2 (coupled-maps spec P3, a chaotic panel; d = 4)
function clv_starts(rng, d::Integer, m::Integer, ::Type{T}) where {T}
    return [SVector{d, T}(ntuple(_ -> T(2 * pi * rand(rng)), Val(d))) for _ in 1:m]
end
const CLV_E = 1 / 8
clv_orbit_C() = energy_start(HenonHeilesFlow(), CLV_E, -0.015, 0.0)[1]               # henon-heiles Equation 9, orbit C (sticky chaotic)
clv_table(::Type{T}, key::Symbol) where {T} = key === :o8 ? KahanLi(T, 8) : key === :o4 ? ForestRuth(T) : key === :o2 ? Leapfrog(T) : error("unknown table $key")
const CLV_H = 1 / 20
clv_im(::Type{T}, key::Symbol = :o8, h = CLV_H) where {T} = IntegratorMap(clv_table(T, key), HenonHeilesFlow(), T(h))
# Orbit points x_0 .. x_n by iterating `forward` (bitwise the stream's orbit, L1 S2).
function clv_orbit(map, x0, n::Integer)
    pts = Vector{typeof(x0)}(undef, n + 1); pts[1] = x0; z = x0
    for i in 1:n
        z = forward(map, z); pts[i + 1] = z
    end
    return pts
end
clv_widen(p::SVector{d}, ::Type{S}) where {d, S} = SVector{d, S}(ntuple(i -> clv_conv(S, clv_big(p[i])), Val(d)))

# ------------------------------------------------------------------ the independent reference algorithm (LinearAlgebra only, any element type): Ginelli's four steps (spec V2, V6, V8)
# Per-interval tangent products along given step points (pts has tau * nint + 1 points; later steps on the left, V1).
function clv_interval_tangents(mapS, pts::AbstractVector, tau::Integer, nint::Integer, ::Type{S}) where {S}
    out = Vector{Matrix{S}}(undef, nint)
    for n in 1:nint
        M = Matrix{S}(I, length(pts[1]), length(pts[1]))
        for s in 1:tau
            M = Matrix{S}(tangent(mapS, clv_widen(pts[(n - 1) * tau + s], S))) * M
        end
        out[n] = M
    end
    return out
end
# QR with positive diagonal in the element type.
function clv_qr(M::AbstractMatrix)
    F = qr(M); Q = Matrix(F.Q); R = Matrix(F.R); s = [R[i, i] < 0 ? -one(eltype(R)) : one(eltype(R)) for i in 1:size(R, 2)]
    return Q * Diagonal(s), Diagonal(s) * R[1:size(R, 2), :]
end
clv_colnorms(X) = [norm(view(X, :, j)) for j in 1:size(X, 2)]
# One backward step. mode: :ok (back-substitution with the upper triangle, unit 2-norm), :fwd (control: multiply by R), :lower (control: lower triangle), :norm1 (control: 1-norm).
function clv_back(R::AbstractMatrix, C::AbstractMatrix, mode::Symbol)
    X = mode === :fwd ? R * C : mode === :lower ? LowerTriangular(R) \ C : UpperTriangular(R) \ C
    nr = mode === :norm1 ? [sum(abs, view(X, :, j)) for j in 1:size(X, 2)] : clv_colnorms(X)
    return X ./ transpose(nr), (mode === :fwd ? nr : 1 ./ nr)         # Cn and Lambda (the growth factor of V6: 1 / ||R^-1 c||)
end
# The whole reference: forward transient nf, forward recording N = ns + nb intervals, backward from C0 (= C_N). Returns G_n, R_n, C_n for n = 0..ns (V = G C) and Lambda_i(n) for n = 0..ns
# (the last column NaN when nb = 0), all in the element type of the tangents.
function clv_ref(Ms::Vector{<:AbstractMatrix{S}}, G0::AbstractMatrix, C0::AbstractMatrix, nf::Integer, ns::Integer, nb::Integer; mode::Symbol = :ok) where {S}
    k = size(G0, 2); N = ns + nb; G = Matrix{S}(G0)
    for i in 1:nf
        G = clv_qr(Ms[i] * G)[1]
    end
    Gs = Vector{Matrix{S}}(undef, ns + 1); Rs = Vector{Matrix{S}}(undef, ns + 1); Rall = Vector{Matrix{S}}(undef, N)
    for n in 0:(N - 1)
        n <= ns && (Gs[n + 1] = G)
        Q, R = clv_qr(Ms[nf + n + 1] * G); Rall[n + 1] = R; G = Q
    end
    C = Matrix{S}(C0); Cs = Vector{Matrix{S}}(undef, ns + 1); Lam = fill(S(NaN), k, ns + 1)
    nb == 0 && (Cs[ns + 1] = copy(C))
    for n in (N - 1):-1:0
        C, lam = clv_back(Rall[n + 1], C, mode)
        n <= ns && (Cs[n + 1] = C; Lam[:, n + 1] = lam)
    end
    for n in 0:min(ns, N - 1)
        Rs[n + 1] = Rall[n + 1]
    end
    nb == 0 && (Rs[ns + 1] = fill(S(NaN), k, k))
    # forward exponents of the window n = 0 .. ns - 1 (per interval)
    lnR = [sum(log(Rall[n + 1][i, i]) for n in 0:(ns - 1)) for i in 1:k]
    return (G = Gs, R = Rs, C = Cs, Lambda = Lam, lnR = lnR)
end
# The reference on a map: orbit from x0 (any element type), tangent products in S along the realized orbit, W0 and C0 widened.
function clv_ref_run(mapS, pts, tau, nf, ns, nb, G0, C0, ::Type{S}; mode::Symbol = :ok) where {S}
    Ms = clv_interval_tangents(mapS, pts, tau, nf + ns + nb, S)
    return clv_ref(Ms, map(x -> clv_conv(S, clv_big(x)), Matrix(G0)), map(x -> clv_conv(S, clv_big(x)), Matrix(C0)), nf, ns, nb; mode = mode)
end
# V_n = G_n C_n for the reference.
clv_refV(ref, n) = ref.G[n + 1] * ref.C[n + 1]
# Slice of a stored 3-array (samples n = 0 .. n_s are the third index n + 1).
clv_at(A::AbstractArray{<:Any, 3}, n::Integer) = A[:, :, n + 1]

# ------------------------------------------------------------------ health error of a call, or nothing; allocation measure
function clv_health_error(f)
    try
        f()
        return nothing
    catch e
        e isa NumericalHealthError ? e : rethrow()
    end
end
function clv_alloc_bytes(f, args...; kw...)
    f(args...; kw...)
    return @allocated f(args...; kw...)
end
clv_J(d::Integer) = (N = d ÷ 2; [zeros(N, N) Matrix{Float64}(I, N, N); -Matrix{Float64}(I, N, N) zeros(N, N)])   # canonical matrix, ordering (q, p)
clv_omega(u, w, d::Integer) = (N = d ÷ 2; sum(clv_big(u[i]) * clv_big(w[N + i]) - clv_big(u[N + i]) * clv_big(w[i]) for i in 1:N))

# ------------------------------------------------------------------ Henon-Heiles chaotic starts at E = 1/8 (henon-heiles criterion 10: chaotic when SALI < 1e-8 at 2000; package horizon 1000)
const CLV_W2 = SMatrix{4, 2, Float64}(1, 0, 0, 0, 0, 0, 1, 0)           # SALI frame of the printed runs (henon-heiles Equation 9)
function clv_sali(z0, Tend)
    im = clv_im(Float64); n = round(Int, Tend / CLV_H)
    res = propagate(im, [z0], n; k = 2, tau = 20, seed = 1, reducers = (SALI(),), W0 = [CLV_W2])
    return res.finals[1][1].sali
end
# The first `m` chaotic starts drawn by random_energy_start from the recorded stream clv_seed(0, seedoff), memoized.
function clv_chaotic_starts(m::Integer, seedoff::Integer)
    return get!(CLV_CACHE, (:chaotic, m, seedoff)) do
        flow = HenonHeilesFlow(); a, b = section_interval(flow, CLV_E); w = sqrt(2 * CLV_E)
        rng = StableRNG(clv_seed(0, seedoff)); Tc = clv_size(1000.0, 2000.0); out = SVector{4, Float64}[]; tries = 0
        while length(out) < m && tries < 50 * m + 200
            tries += 1
            z, st = random_energy_start(flow, CLV_E, (a, b), (-w, w), rng)
            st === :ok || continue
            clv_sali(z, Tc) < 1.0e-8 && push!(out, z)
        end
        out
    end
end

# ------------------------------------------------------------------ module drivers (called inside testsets only; they use the proposed API names)
# One orbit with explicit start data. `cfg` is a NamedTuple of ClvGinelli keywords (storage, block); other keywords go to propagate_clv_orbit.
function clv_go(map, x0, nf, ns, nb; tau = 1, h = 1.0, k = nothing, W0 = nothing, C0 = nothing, alg = nothing, outputs = (:V, :G, :R), cfg = (;), kw...)
    d = dimension(map)
    clv = alg === nothing ? ClvGinelli(map; tau = tau, h = h, cfg...) : ClvGinelli(map; tau = tau, h = h, orthonormalizer = alg, cfg...)
    kk = k === nothing ? (W0 === nothing ? d : size(W0, 2)) : k
    return propagate_clv_orbit(clv, x0, nf, ns, nb; k = kk, seed = 1, W0 = W0, C0 = C0, outputs = outputs, kw...)
end
# Pilot exponents (per unit time) of a map by the reference algorithm, Float64, window of `ns` handoffs from x0, and the transient n* of the independence criterion.
function clv_pilot(map, x0, tau::Integer, h; ns::Integer = 300)
    d = length(x0); pts = clv_orbit(map, x0, tau * (ns + 1))
    G0 = clv_frame(1, 1, d, d, Float64); C0 = clv_eye(d, d, Float64)
    return clv_ref_run(map, pts, tau, 0, ns, 1, G0, C0, Float64).lnR ./ (ns * tau * h)
end
# n* = ceil(c / (g_min Dt)) with c = 40 (spec proposes 25; test-side reading U3), g_min the smallest gap of the pilot spectrum, the pair (2, 3) of a Hamiltonian flow skipped; capped at nmax.
function clv_nstar(lam, Dt; skip23::Bool = false, c = 40.0, nmax = clv_size(3000, 6000))
    gaps = [lam[i] - lam[i + 1] for i in 1:(length(lam) - 1)]
    skip23 && deleteat!(gaps, 2)
    return min(nmax, ceil(Int, c / (minimum(gaps) * Dt)))
end
# Largest distance between two CLV matrices (d x k) over single directions `singles` and the span of the columns `block`.
function clv_vdist(VA, VB; singles = 1:size(VA, 2), block = nothing)
    w = isempty(singles) ? 0.0 : maximum(i -> clv_delta(VA[:, i], VB[:, i]), singles)
    block === nothing || (w = max(w, clv_sin_sub(VA[:, block], VB[:, block])))
    return w
end
# Window conditioning G_w of the forward QR: max over the stored R factors (spec GS G10), Float64 and Double64 factors.
clv_Gmax(R::AbstractArray{<:Any, 3}, n::Integer) = maximum(j -> Float64(window_Gw(SMatrix{size(R, 1), size(R, 2), eltype(R)}(R[:, :, j]))), 1:n)

# BigFloat copy of a package map (the tangents along the realized Float64 or Double64 orbit are recomputed in 256 bits: the device of L1 S3).
clv_bigmap(m::StandardMap) = StandardMap(clv_big(m.K))
clv_bigmap(m::StandardMapChain) = StandardMapChain(m.N, clv_big(m.K), clv_big(m.eps), m.ring)
clv_bigmap(im::IntegratorMap) = IntegratorMap(clv_table(BigFloat, :o8), im.flow, clv_big(im.tau))

map_to(::Type{T}, z::SVector{d}) where {T, d} = SVector{d, T}(ntuple(i -> T(z[i]), Val(d)))
