# Helpers for the gali suite (prefix ga_). Definitions only: nothing here calls a new Indicators name (GALI, gali_values, gali_value, gali_floor)
# at load time or inside a helper, so a missing module errors the testsets that use it and nothing else. Spec:
# workspaces/module-build/stages/01-spec/output/gali-spec.md; API proposal and conventions: stages/02-tests/output/gali-test-plan.md.
using Test, Chaoslab, StaticArrays, LinearAlgebra, StableRNGs, DoubleFloats

const GA_SEED = 20261012                  # base seed, _config/numerics-defaults.toml key `seed`
const GA_BASE = GA_SEED + 40000           # this slug uses 20261012 + 40000 .. 40999; testset AC<n> uses 50 n + sub (sub < 50); AC20 uses n = 19, the API tests n = 0
ga_seed(ac::Integer, sub::Integer = 0) = GA_BASE + 50 * ac + sub
const GA_TYPES = (Float64, Double64)      # Float32 is out (D24); BigFloat only for oracle runs and the oracle-grade module runs
const GA_BITS = 256                       # numerics-defaults.toml key `bigfloat_bits`
# Package suite sizes (D38, 3 minute ceiling) or the full sizes of the spec: verify/gali.jl sets GA_FULL_RUN = true before including runtests.jl.
# Every reduction is marked "D38" at its use. The testsets are top-level code (not functions), as in the other suites.
const GA_FULL = isdefined(Main, :GA_FULL_RUN) && Main.GA_FULL_RUN
ga_size(small, full) = GA_FULL ? full : small

ga_eps(::Type{T}) where {T} = Float64(effective_eps(T))     # Float64 eps, Double64 2^-100, BigFloat eps at the precision in force
ga_meas(::Type{T}, k::Integer) where {T} = 2 * k * ga_eps(T)  # B_meas(T, k) = 2 k eps (B9); the test-side copy, gali_floor is compared to it

# ------------------------------------------------------------------ BigFloat conversions and the exact evaluation

ga_big(x::Double64) = BigFloat(x.hi) + BigFloat(x.lo)       # exact; call inside setprecision(BigFloat, bits)
ga_big(x::AbstractFloat) = BigFloat(x)
ga_big(x::Integer) = BigFloat(x)

# BigFloat into T (nearest): Double64 as hi + lo.
ga_to(::Type{Float64}, b) = Float64(b)
ga_to(::Type{Double64}, b) = (h = Float64(b); Double64(h, Float64(b - h)))
ga_to(::Type{BigFloat}, b) = BigFloat(b; precision = GA_BITS)

# |value - ref| as Float64 (value exact in its type, ref a BigFloat), evaluated at `bits` bits.
ga_err(value, ref; bits = GA_BITS) = setprecision(BigFloat, bits) do
    Float64(abs(ga_big(value) - ref))
end

ga_sm(M::AbstractMatrix) = SMatrix{size(M, 1), size(M, 2), eltype(M)}(M)
ga_sv(v::AbstractVector) = SVector{length(v), eltype(v)}(v)

# Unit columns of a BigFloat matrix at the precision in force.
ga_unitcols(U::AbstractMatrix{BigFloat}) = mapslices(c -> c / sqrt(sum(abs2, c)), U; dims = 1)

# Volumes G_1..G_K of the leading columns of the unit matrix U (derived route B2 iv): the product of the norms of the successive residuals of a
# twice-orthogonalised modified Gram-Schmidt, in BigFloat. Call inside setprecision.
function ga_vols(U::AbstractMatrix{BigFloat})
    K = size(U, 2)
    Q = Vector{BigFloat}[]
    vols = Vector{BigFloat}(undef, K)
    v = one(BigFloat)
    for j in 1:K
        w = U[:, j]
        for _ in 1:2, q in Q
            w = w - dot(q, w) * q
        end
        r = sqrt(sum(abs2, w))
        v *= r
        vols[j] = v
        push!(Q, iszero(r) ? w : w / r)
    end
    return vols
end

# G_2..G_K of the raw d x K matrix W (columns normalised first, as the module does), exact for the given entries.
function ga_exact(W::AbstractMatrix; bits = GA_BITS)
    return setprecision(BigFloat, bits) do
        U = ga_unitcols(BigFloat[ga_big(W[i, j]) for i in 1:size(W, 1), j in 1:size(W, 2)])
        ga_vols(U)[2:end]
    end
end

# The four exact routes of B2 for a unit matrix U (BigFloat, call inside setprecision): minors, Gram determinant, singular values, QR.
function ga_combs(n::Int, k::Int)
    k == 0 && return [Int[]]
    n < k && return Vector{Int}[]
    return vcat([vcat([1], c .+ 1) for c in ga_combs(n - 1, k - 1)], [c .+ 1 for c in ga_combs(n - 1, k)])
end
ga_route_minors(U) = sqrt(sum(rows -> det(U[rows, :])^2, ga_combs(size(U, 1), size(U, 2))))
ga_route_gram(U) = sqrt(det(U' * U))
ga_route_svd(U) = prod(sqrt.(max.(ga_jacobi_eigvals(U' * U), zero(eltype(U)))))
ga_route_qr(U) = ga_vols(U)[end]

# Cyclic Jacobi eigenvalues of a symmetric BigFloat matrix (no GenericLinearAlgebra dependency); call inside setprecision.
function ga_jacobi_eigvals(A0::AbstractMatrix{BigFloat}; digits = 100)
    A = Matrix{BigFloat}(A0)
    n = size(A, 1)
    scale = sum(abs2, A)
    for _ in 1:200
        off = sum(abs2(A[i, j]) for i in 1:n for j in 1:n if i != j)
        off <= scale * BigFloat(10)^(-digits) && break
        for p in 1:(n - 1), q in (p + 1):n
            iszero(A[p, q]) && continue
            th = (A[q, q] - A[p, p]) / (2 * A[p, q])
            t = (th >= 0 ? one(th) : -one(th)) / (abs(th) + sqrt(th^2 + 1))
            c = 1 / sqrt(t^2 + 1)
            s = t * c
            for k in 1:n
                akp, akq = A[k, p], A[k, q]
                A[k, p] = c * akp - s * akq
                A[k, q] = s * akp + c * akq
            end
            for k in 1:n
                apk, aqk = A[p, k], A[q, k]
                A[p, k] = c * apk - s * aqk
                A[q, k] = s * apk + c * aqk
            end
        end
    end
    return sort!([A[i, i] for i in 1:n]; rev = true)
end

# ------------------------------------------------------------------ random sets and constructed sets

# Random vector in T with a full-width mantissa (Double64 gets a nonzero low part).
function ga_randvec(rng, ::Type{T}, d::Int) where {T}
    return T[(x = T(randn(rng)); T === Double64 && (x += ldexp(T(randn(rng)), -55)); x) for _ in 1:d]
end

# d x K set with raw column norms spread over 10^-3 .. 10^3 (criterion 1).
function ga_randset(rng, ::Type{T}, d::Int, K::Int; spread::Bool = true) where {T}
    M = hcat((ga_randvec(rng, T, d) * T(spread ? 10.0^(6 * rand(rng) - 3) : 1.0) for _ in 1:K)...)
    return ga_sm(M)
end

# Random orthogonal d x d matrix in BigFloat (Gram-Schmidt of a Gaussian matrix), call inside setprecision.
function ga_randorth(rng, d::Int)
    M = BigFloat[BigFloat(randn(rng)) for _ in 1:d, _ in 1:d]
    Q = Vector{BigFloat}[]
    for j in 1:d
        w = M[:, j]
        for _ in 1:2, q in Q
            w = w - dot(q, w) * q
        end
        push!(Q, w / sqrt(sum(abs2, w)))
    end
    return hcat(Q...)
end

# Unit-diagonal Gram matrix G (k x k, BigFloat) -> d x k matrix of unit vectors with that Gram matrix, rotated by the orthogonal Q (Cholesky factor
# rows as vectors; a zero or negative pivot, a singular Gram matrix, gives a zero component). Call inside setprecision.
function ga_gram_vectors(G::AbstractMatrix{BigFloat}, d::Int, Q::AbstractMatrix{BigFloat})
    k = size(G, 1)
    L = zeros(BigFloat, k, k)
    for i in 1:k
        for j in 1:i
            s = G[i, j] - sum((L[i, m] * L[j, m] for m in 1:(j - 1)); init = zero(BigFloat))
            if i == j
                L[i, i] = sqrt(max(s, zero(BigFloat)))
            else
                L[i, j] = iszero(L[j, j]) ? zero(BigFloat) : s / L[j, j]
            end
        end
    end
    V = zeros(BigFloat, d, k)
    V[1:k, :] = L'
    return Q * V
end

# Near-dependent unit sets in d = 4 (criterion 6), BigFloat, call inside setprecision: k = 2 angle theta; k = 3 two columns at angle theta and a third
# orthogonal; k = 4 a quadruple u_j = cos(theta) p_j + sin(theta) q_j with p_j in the plane (e1, e2), q_j in the plane (e3, e4), distinct angles.
function ga_neardep(k::Int, theta::BigFloat)
    c, s = cos(theta), sin(theta)
    z, o = zero(BigFloat), one(BigFloat)
    k == 2 && return BigFloat[o c; z s; z z; z z]
    k == 3 && return BigFloat[o c z; z s z; z z o; z z z]
    phis, psis = (0.3, 0.9, 1.5, 2.1), (0.2, 1.1, 2.0, 2.9)
    cols = [[c * cos(BigFloat(a)), c * sin(BigFloat(a)), s * cos(BigFloat(b)), s * sin(BigFloat(b))] for (a, b) in zip(phis, psis)]
    return hcat(cols...)
end

# ------------------------------------------------------------------ the systems (Y3, Y4) and the cases

struct GACase
    kind::Symbol                  # :y3 (p = (nu, kappa, mu)) or :y4 (p = (K, eps, 0), N sites)
    N::Int
    p::NTuple{3, Float64}
    ring::Bool
end
ga_y3(nu, kappa, mu) = GACase(:y3, 2, (Float64(nu), Float64(kappa), Float64(mu)), false)
ga_y4(N, K, e; ring = false) = GACase(:y4, N, (Float64(K), Float64(e), 0.0), ring)
ga_d(c::GACase) = c.kind === :y3 ? 4 : 2 * c.N

function ga_sys(::Type{T}, c::GACase) where {T}
    c.kind === :y3 && return FroeschleMap(T(c.p[1]), T(c.p[2]), T(c.p[3]))
    return StandardMapChain(c.N, T(c.p[1]), T(c.p[2]); boundary = c.ring ? :ring : :open)
end
function ga_wsys(::Type{T}, c::GACase) where {T}
    c.kind === :y3 && return WrappedFroeschleMap(T(c.p[1]), T(c.p[2]), T(c.p[3]))
    return WrappedStandardMapChain(c.N, T(c.p[1]), T(c.p[2]); boundary = c.ring ? :ring : :open)
end
ga_state(::Type{T}, c::GACase, z) where {T} = SVector{ga_d(c), T}(ntuple(i -> T(z[i]), ga_d(c)))
ga_frame(seed, j, ::Type{T}, d::Int, K::Int) where {T} = initial_deviations(seed, j, Val(d), Val(K), T)   # first K columns of one Haar frame (nesting)

const GA_Y3A = ga_y3(0.5, 0.1, 1.0e-3)        # parameters of orbits C and D (coupled-maps item 3), item 1 labelling (positions first)
const GA_Y3B = ga_y3(2.0, 0.5, 0.1)           # the chaotic panel of the spec (criteria 9, 13, 15)
const GA_ORBIT_C = (0.5, 0.5, 0.0, 0.0)
const GA_ORBIT_D = (3.0, 0.5, 0.0, 0.0)

ga_draw_chaotic(rng, c::GACase) = [2 * Float64(pi) * rand(rng) for _ in 1:ga_d(c)]       # uniform in [0, 2 pi)^d
function ga_draw_regular(rng, c::GACase)                                                 # small amplitude (criterion 14)
    d = ga_d(c)
    c.kind === :y3 ? [(2 * rand(rng) - 1) * (i <= 2 ? 0.4 : 0.12) for i in 1:4] : [(2 * rand(rng) - 1) * 0.3 for _ in 1:d]
end

# Orbit points pts[i] = x_{i-1} (the point where step i evaluates the tangent), i = 1..n, and x_n; forward in the element type of x0.
function ga_pts(sys, x0, n::Integer)
    pts = Vector{typeof(x0)}(undef, n)
    z = x0
    for i in 1:n
        pts[i] = z
        z = forward(sys, z)
    end
    return pts, z
end

# Checkpoints n <= N: multiples of 10 at 10 per decade (D38 reduction of "every handoff that is a multiple of 10"); all multiples with GA_FULL.
ga_checkpoints(N::Int) = GA_FULL ? collect(10:10:N) :
                         sort(unique(filter(n -> n % 10 == 0 && n <= N, vcat((d .* (1:10) for d in (1, 10, 100, 1000, 10000))...))))

# ------------------------------------------------------------------ the 256-bit reference along an orbit (O2, B10)

# One-step tangent of the system at the point pt (held in T), evaluated in BigFloat at the precision in force, as a dense matrix.
function ga_tanB(sysB, c::GACase, pt)
    z = c.kind === :y3 ? SVector{4, BigFloat}(ntuple(i -> ga_big(pt[i]), 4)) : BigFloat[ga_big(v) for v in pt]
    return Matrix{BigFloat}(tangent(sysB, z))
end
# Any system (the standard map, the time-dependent test map): the state is widened as it is.
ga_tanB(sysB, ::Nothing, pt) = Matrix{BigFloat}(tangent(sysB, map(ga_big, pt)))

# Products Y_n W0 at `bits` bits (tangents at the widened realised points pts[i]); returns Float64 ln G_k(n) for every step n = 1..N (rows k = 1..K, G_1 = 1)
# and, at the steps in `keep` (ascending), the BigFloat G_k (k = 2..K), the column sums S_i = ln norm(Y w_i) and, for K = d, the log form reference
# L = ln abs(det W0) + sum_i ln abs(det M_i) - sum_i S_i. That is ln G_d of the exact product (det is multiplicative; each M_i is well conditioned, so its determinant
# is accurate), valid to any depth: the volume G itself, taken from unit columns rounded to `bits` bits, is only meaningful above 2^-bits.
function ga_ref(c, pts, W0::AbstractMatrix, keep; bits = GA_BITS, all_steps::Bool = true, sysB = nothing)
    return setprecision(BigFloat, bits) do
        sysB = sysB === nothing ? ga_sys(BigFloat, c) : sysB
        d, K = size(W0)
        Y = BigFloat[ga_big(W0[i, j]) for i in 1:d, j in 1:K]
        N = length(pts)
        lnG = fill(NaN, K, N)
        G, S, Lt = Dict{Int, Vector{BigFloat}}(), Dict{Int, Vector{BigFloat}}(), Dict{Int, BigFloat}()
        ks = Set(keep)
        lndet = K == d ? log(abs(det(Y))) : zero(BigFloat)
        for i in 1:N
            M = ga_tanB(sysB, c, pts[i])
            K == d && (lndet += log(abs(det(M))))
            Y = M * Y
            if all_steps || i in ks
                norms = [sqrt(sum(abs2, view(Y, :, j))) for j in 1:K]
                vols = ga_vols(Y ./ reshape(norms, 1, K))
                lnG[:, i] = Float64.(log.(vols))
                if i in ks
                    G[i] = vols[2:end]
                    S[i] = log.(norms)
                    K == d && (Lt[i] = lndet - sum(S[i]))
                end
            end
        end
        (lnG = lnG, G = G, S = S, L = Lt, K = K, d = d, keep = collect(Int, keep))
    end
end

# ln G_d at the kept steps as BigFloat (the log form reference, K = d); `bits` kept for the callers that pass it.
ga_Lref(ref, n::Int; bits = GA_BITS) = ref.L[n]

# Psi_n^(i) (B10, L1 Amendment D35) for the columns i = 1..K of W0 at the steps ns, from the Float64 one-step tangents Ms[1..N]. The matrix norm is the spectral
# norm (opnorm). Scaled recursions, so no overflow for chaotic orbits. Returns a vector of K-vectors, one per n.
function ga_psis(Ms::AbstractVector{<:AbstractMatrix{Float64}}, W0::AbstractMatrix, ns)
    d, K = size(W0)
    N = maximum(ns)
    lw = zeros(K, N + 1)                                  # lw[i, j + 1] = ln norm(w_j^(i)), w_0 unit
    for i in 1:K
        w = Float64.(W0[:, i])
        a = norm(w)
        w = w / a
        acc = log(a)
        lw[i, 1] = acc
        for k in 1:N
            w = Ms[k] * w
            a = norm(w)
            acc += log(a)
            w = w / a
            lw[i, k + 1] = acc
        end
    end
    return map(ns) do n
        P = Matrix{Float64}(I, d, d)
        lP = 0.0
        lnP = zeros(n)
        for j in n:-1:1
            Q = P * Ms[j]
            f = norm(Q)
            P = Q / f
            lP += log(f)
            lnP[j] = lP + log(opnorm(P))
        end
        [sum(j -> exp(lnP[j] + lw[i, j] - lw[i, n + 1]), 1:n) for i in 1:K]
    end
end

# Lambda_n^(k) = sum_{j=1..n} G_k,ref(n) / G_k,ref(j) (the sali A15 form with k columns); the ratio is clipped at e^700 (a scale, not a value).
ga_Lambda(lnGk::AbstractVector{Float64}, n::Int) = sum(j -> exp(min(lnGk[n] - lnGk[j], 700.0)), 1:n)
# B10: direct form bound B_G^k(n) and log form bound B_L(n); epsT = effective eps of the element type, psi the K-vector Psi_n^(i).
# The multiplier of B_L is 32 here, the spec proposes 8: the reference stub of stage 02 realises error / (8 eps sum Psi) up to 1.54 on strongly chaotic starts
# (the per-handoff rounding of ln alpha_i is a few eps per column, more than the 2 eps of the direction error the 8 was derived from), so the gate carries a
# factor 4; stage 04 reports the realised ratios (P5).
const GA_BL_MULT = 32
ga_BG(epsT, k, Lam, Gref, psi) = 4 * k * epsT * (Lam + Gref * sum(view(psi, 1:k))) + 2 * k * epsT
ga_BL(epsT, psi) = GA_BL_MULT * epsT * sum(psi)

const GA_CACHE = Dict{Any, Any}()

# Reference data of one start for any system: realised orbit in the element type of x0, 256-bit (or `bits`) reference at the checkpoints cps, Psi and Lambda at the
# checkpoints (Lambda for every k = 1..K). `c` is a GACase or nothing (then `sysB` must give the BigFloat system, used for the tangents at `bits` bits).
function ga_bundle_core(sysT, c, x0, W0, N::Int, cps; bits = GA_BITS, sysB = nothing)
    pts, _ = ga_pts(sysT, x0, N)
    ref = ga_ref(c, pts, W0, cps; bits = bits, sysB = sysB)
    Ms = [Float64.(Matrix(tangent(sysT, p))) for p in pts]
    psi = ga_psis(Ms, W0, cps)
    K = size(W0, 2)
    lam = [[ga_Lambda(ref.lnG[k, :], n) for k in 1:K] for n in cps]
    return (x0 = x0, W0 = W0, cps = collect(Int, cps), ref = ref, psi = psi, lam = lam, d = size(W0, 1))
end

# Everything a same-orbit criterion needs for one start of Y3 or Y4: realised orbit in T, seeded frame (K = d columns unless W0 is given), reference, Psi and Lambda
# at the checkpoints (the multiples of 10 at 10 per decade, or `cps`). `z` is the start as a tuple of Float64 (widened exactly into T); `fid` = (frame seed, orbit id).
function ga_bundle(c::GACase, ::Type{T}, z, fid, N::Int; W0 = nothing, cps = nothing) where {T}
    return get!(GA_CACHE, (:bundle, c, T, Tuple(z), fid, N, W0 === nothing ? nothing : Tuple(vec(W0)), cps === nothing ? nothing : Tuple(cps))) do
        d = ga_d(c)
        Wf = W0 === nothing ? ga_frame(fid[1], fid[2], T, d, d) : W0
        b = ga_bundle_core(ga_sys(T, c), c, ga_state(T, c, z), Wf, N, cps === nothing ? ga_checkpoints(N) : cps; bits = T === BigFloat ? 2 * GA_BITS : GA_BITS)
        merge(b, (c = c, T = T))
    end
end

# ------------------------------------------------------------------ reducers used by the tests (protocol of tangent-core)

# Wraps a reducer: forwards init, step, failed, error to `inner`; after the inner step it records reducer_final(inner) merged with `eq12` (are the first
# two columns of the replaced W bitwise equal) at the handoffs in `steps` (or at every handoff when `all`). The stream has no emission channel (L1 Q9),
# so this is how a series is read. CPU only (the state holds a Vector): use it with propagate_orbit.
struct ga_Probe{R}
    inner::R
    steps::Vector{Int}
    all::Bool
end
ga_Probe(inner, steps::AbstractVector{<:Integer} = Int[]; all::Bool = false) = ga_Probe(inner, collect(Int, steps), all)

Chaoslab.reducer_init(r::ga_Probe, ::Type{T}, v::Val, k::Val) where {T} = (Chaoslab.reducer_init(r.inner, T, v, k), 1, Any[])
function Chaoslab.reducer_step(r::ga_Probe, st, n::Int, x, W)
    inner, pos, snaps = st
    inner2, W2 = Chaoslab.reducer_step(r.inner, inner, n, x, W)
    if r.all || (pos <= length(r.steps) && n == r.steps[pos])
        eq12 = size(W2, 2) >= 2 && isequal(W2[:, 1], W2[:, 2])
        push!(snaps, merge(Chaoslab.reducer_final(r.inner, inner2), (eq12 = eq12,)))
        pos += 1
    end
    return (inner2, pos, snaps), W2
end
Chaoslab.reducer_failed(r::ga_Probe, st) = Chaoslab.reducer_failed(r.inner, st[1])
Chaoslab.reducer_error(r::ga_Probe, st, orbit, backend) = Chaoslab.reducer_error(r.inner, st[1], orbit, backend)
Chaoslab.reducer_final(r::ga_Probe, st) = (snaps = st[3], inner = Chaoslab.reducer_final(r.inner, st[1]))

# Read-only reducer recording the orbit point at every handoff (placed before the replacer).
struct ga_Trace end
Chaoslab.reducer_init(::ga_Trace, ::Type, ::Val, ::Val) = Any[]
function Chaoslab.reducer_update(::ga_Trace, st, n::Int, x, W)
    push!(st, x)
    return st, W
end
Chaoslab.reducer_final(::ga_Trace, st) = (xs = st,)

# One orbit through the per-orbit body of the stream with the reducer wrapped in a probe; returns reducer_final(probe) (fields snaps, inner).
function ga_run1(reducer, sys, x0, W0, N::Integer, tau::Integer; steps = Int[], all::Bool = false, orbit::Integer = 1)
    probe = ga_Probe(reducer, steps; all = all)
    _, _, st = propagate_orbit(sys, x0, W0, N, tau, (probe,); orbit = orbit)
    return Chaoslab.reducer_final(probe, st[1])
end

# ------------------------------------------------------------------ test maps (never in src/)

# Fault map on the state (n, y, y, y, flag), dimension 4: the counter advances by one per step; only flagged orbits misbehave.
# kind :zero: zero tangent (every column collapses at the first handoff); :col3: tangent diag(1, 1, 0, 1) (column 3 collapses when the frame is I);
# :nan: NaN in the tangent at counter value `at`.
struct ga_Fault <: Chaoslab.AbstractMap
    kind::Symbol
    at::Int
end
Chaoslab.dimension(::ga_Fault) = 4
Chaoslab.forward(::ga_Fault, z) = SVector(z[1] + 1, z[2], z[3], z[4], z[5])
function Chaoslab.tangent(m::ga_Fault, z::SVector{5, T}) where {T}
    z[5] == 1 || return one(SMatrix{4, 4, T})
    m.kind === :zero && return zero(SMatrix{4, 4, T})
    m.kind === :col3 && return SMatrix{4, 4, T}(1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1)
    return SMatrix{4, 4, T}(z[1] == m.at ? T(NaN) : one(T), 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1)
end

# Map whose forward and tangent throw: proves a rejection happens before any work.
struct ga_Boom <: Chaoslab.AbstractMap
    d::Int
end
Chaoslab.dimension(m::ga_Boom) = m.d
Chaoslab.forward(::ga_Boom, z) = error("ga_Boom: forward was called")
Chaoslab.tangent(::ga_Boom, z) = error("ga_Boom: tangent was called")

# Time-dependent tangent: the state is (z_1..z_d, counter); step i applies mats[i]; the orbit itself is constant. CPU only (a Vector field).
struct ga_Seq{T} <: Chaoslab.AbstractMap
    d::Int
    mats::Vector{Matrix{T}}
end
Chaoslab.dimension(m::ga_Seq) = m.d
Chaoslab.forward(m::ga_Seq, z::SVector{s, T}) where {s, T} = SVector{s, T}(ntuple(i -> i == s ? z[s] + 1 : z[i], s))
Chaoslab.tangent(m::ga_Seq, z::SVector{s, T}) where {s, T} = SMatrix{m.d, m.d, T}(m.mats[Int(round(Float64(z[s]))) + 1])

# Elementary symplectic BigFloat factors of dimension 2n (call inside setprecision): shears with symmetric blocks, block-diagonal [A, A^-T], rotations.
function ga_symp_factor(rng, n::Int; amp = 1.0)
    S = BigFloat[0 for _ in 1:n, _ in 1:n]
    for i in 1:n, j in i:n
        S[i, j] = S[j, i] = BigFloat(2 * rand(rng) - 1) * (1.5 * amp)
    end
    Z, I_ = zeros(BigFloat, n, n), Matrix{BigFloat}(I, n, n)
    kind = rand(rng, 1:4)
    kind == 1 && return [I_ S; Z I_]
    kind == 2 && return [I_ Z; S I_]
    if kind == 3
        A = I_ + BigFloat[BigFloat(2 * rand(rng) - 1) * (0.6 * amp) for _ in 1:n, _ in 1:n]
        return [A Z; Z inv(A)']
    end
    th = BigFloat[BigFloat(2 * rand(rng) - 1) * 3 for _ in 1:n]
    C, Sn = Diagonal(cos.(th)), Diagonal(sin.(th))
    return [Matrix(C) Matrix(Sn); -Matrix(Sn) Matrix(C)]
end

# ------------------------------------------------------------------ negative controls (test-side doubles; they do not use the new module)

# Skip-normalisation double: the volume of the leading k raw columns, not divided by the column norms.
function ga_ctl_noscale(W::AbstractMatrix{T}, k::Int) where {T}
    G = W[:, 1:k]' * W[:, 1:k]
    return sqrt(max(det(G), zero(T)))
end
# Gram-determinant double: sqrt(det(U' U)) of the leading k unit columns, in T.
function ga_ctl_gram(W::AbstractMatrix{T}, k::Int) where {T}
    U = hcat((W[:, j] / norm(W[:, j]) for j in 1:k)...)
    return sqrt(max(det(U' * U), zero(T)))
end
# SALI of two vectors (min of the two indices, unit vectors), the quantity the SALI-as-GALI double returns.
function ga_sali(a::AbstractVector{T}, b::AbstractVector{T}) where {T}
    u, v = a / norm(a), b / norm(b)
    return min(norm(u - v), norm(u + v))
end

# ------------------------------------------------------------------ exact spectra at the fixed points (O3, B6, B8)

# Laplacian eigenvalues of the open chain or ring (closed forms, BigFloat, call inside setprecision).
ga_laplacian(N::Int, ring::Bool) = ring ? [2 - 2 * cos(2 * BigFloat(pi) * k / N) for k in 0:(N - 1)] : [2 - 2 * cos(k * BigFloat(pi) / N) for k in 0:(N - 1)]

# Exponent of the mode with tangent block [[1 + g, 1], [g, 1]]: ln lambda_+ with lambda + 1/lambda = abs(2 + g) (hyperbolic), 0 (elliptic).
function ga_mode_exponent(g::BigFloat)
    t = abs(2 + g)
    return t > 2 ? log((t + sqrt(t^2 - 4)) / 2) : zero(BigFloat)
end

# l_1 >= ... >= l_d at the origin of Y3 (two decoupled blocks with g = -nu, -kappa) or Y4 (g = K + eps mu_k), pairing l_{d+1-i} = -l_i.
function ga_exponents(c::GACase)
    return setprecision(BigFloat, GA_BITS) do
        gs = c.kind === :y3 ? [-BigFloat(c.p[1]), -BigFloat(c.p[2])] : [BigFloat(c.p[1]) + BigFloat(c.p[2]) * mu for mu in ga_laplacian(c.N, c.ring)]
        ls = sort([ga_mode_exponent(g) for g in gs]; rev = true)
        vcat(ls, -reverse(ls))
    end
end
# E_k = (k - 1) l_1 - sum_{j=2..k} l_j (B6, with the pairing; equal leading exponents drop their terms automatically).
ga_Ek(l, k::Int) = (k - 1) * l[1] - sum(l[2:k])

# Log singular values of the product of the BigFloat(GA_BITS) one-step tangents Ms (matrices applied in order) at the steps ns: the product is accumulated at
# `hp` bits (chosen from the Float64 growth so that the smallest singular value survives), the singular values come from Jacobi eigenvalues of Y' Y.
function ga_logsv(Ms::AbstractVector{<:AbstractMatrix{BigFloat}}, ns)
    d = size(Ms[1], 1)
    lg = 0.0                                              # Float64 estimate of ln sigma_1 at the last step
    P = Matrix{Float64}(I, d, d)
    for n in 1:maximum(ns)
        P = Float64.(Ms[n]) * P
        f = norm(P)
        P /= f
        lg += log(f)
    end
    hp = max(GA_BITS, ceil(Int, 6 * lg / log(2)) + 320)
    return setprecision(BigFloat, hp) do
        Y = Matrix{BigFloat}(I, d, d)
        out = Dict{Int, Vector{BigFloat}}()
        for n in 1:maximum(ns)
            Y = Ms[n] * Y
            if n in ns
                ev = ga_jacobi_eigvals(Y' * Y; digits = floor(Int, hp * log10(2.0)) - 12)
                out[n] = log.(ev) / 2
            end
        end
        out
    end
end

# ------------------------------------------------------------------ small utilities

# Run f(); return the NumericalHealthError it throws (nothing if it returns), rethrow anything else.
function ga_health_error(f)
    try
        f()
        return nothing
    catch e
        e isa NumericalHealthError ? e : rethrow()
    end
end

# Bytes allocated by one propagate_orbit call after a warm-up call, measured inside a function.
function ga_alloc_bytes(sys, x0, W0, N, tau, reducers)
    propagate_orbit(sys, x0, W0, N, tau, reducers)
    return @allocated propagate_orbit(sys, x0, W0, N, tau, reducers)
end

ga_mean(v) = sum(v) / length(v)
ga_var(v) = (m = ga_mean(v); sum(x -> (x - m)^2, v) / (length(v) - 1))
ga_median(v) = (s = sort(collect(v)); n = length(s); isodd(n) ? s[(n + 1) ÷ 2] : (s[n ÷ 2] + s[n ÷ 2 + 1]) / 2)
# L1 S5: the difference of the ensemble means is at most 4 combined standard errors (plus an allowance).
ga_agree(a, b; allowance = 0.0) = abs(ga_mean(a) - ga_mean(b)) <= 4 * sqrt(ga_var(a) / length(a) + ga_var(b) / length(b)) + allowance
# Least-squares slope of v against u.
ga_slope(u, v) = sum((u .- ga_mean(u)) .* (v .- ga_mean(v))) / sum((u .- ga_mean(u)) .^ 2)
