# Test helpers for the precision suite (prefix pr_). Reference implementations are independent of Chaoslab.
using StaticArrays, DoubleFloats, StableRNGs

const PR_SEED = 20261012            # base seed, _config/numerics-defaults.toml
const PR_NMAX = 100                 # reference orbit length
const PR_C = 8.0                    # Equations 4: c = 8
const PR_TOL = 1e-4                 # Equations 4: perturbative horizon

# Exact conversion to BigFloat at the precision currently in force (caller sets it).
pr_big(x::Float64) = BigFloat(x)
pr_big(x::Double64) = BigFloat(x.hi) + BigFloat(x.lo)
pr_big(x::BigFloat) = x

pr_norm(v) = sqrt(v[1]^2 + v[2]^2)

# 200 seeded ICs uniform in [0, 2pi)^2 as Float64.
function pr_ics(seed, n)
    rng = StableRNG(seed)
    [(2pi * rand(rng), 2pi * rand(rng)) for _ in 1:n]
end

# Reference orbit in BigFloat at the caller's precision: p' = p + K sin x, x' = x + p'.
# DT = [[A, 1], [C, 1]] with A = 1 + K cos x, C = K cos x, stored as Float64 (only magnitudes matter).
function pr_reference(K::BigFloat, x0::Float64, p0::Float64, n::Int)
    xs = Vector{BigFloat}(undef, n + 1); ps = similar(xs)
    A = zeros(n + 1); C = zeros(n + 1)
    xs[1] = BigFloat(x0); ps[1] = BigFloat(p0)
    for k in 1:n+1
        c = K * cos(xs[k])
        C[k] = Float64(c); A[k] = Float64(1 + c)
        k == n + 1 && break
        ps[k+1] = ps[k] + K * sin(xs[k])
        xs[k+1] = xs[k] + ps[k+1]
    end
    (; xs, ps, A, C)
end

# M_n and Gamma_n of Equation 3 from a reference orbit; vectors indexed n+1 for n = 0..N.
function pr_bounds(ref)
    N = length(ref.xs) - 1
    M = zeros(N + 1); G = zeros(N + 1)
    m = 1.0
    for k in 0:N
        m = max(m, abs(Float64(ref.xs[k+1])), abs(Float64(ref.ps[k+1])))
        M[k+1] = m
    end
    for n in 1:N
        a, b, c, d = 1.0, 0.0, 0.0, 1.0                  # P = identity, term for k = n-1
        g = 1.0
        for k in n-1:-1:1                                 # P <- P * DT(z_k)
            A = ref.A[k+1]; C = ref.C[k+1]
            a, b, c, d = a * A + b * C, a + b, c * A + d * C, c + d
            f = a^2 + b^2 + c^2 + d^2
            s = sqrt(max(f * f - 4 * (a * d - b * c)^2, 0.0))
            g += sqrt((f + s) / 2)
        end
        G[n+1] = g
    end
    (M, G)
end

# n_H(T): largest n with c eps M_n Gamma_n <= 1e-4 (monotone scan), capped at the reference length.
function pr_nH(M, G, epsT)
    N = length(M) - 1
    n = 0
    while n < N && PR_C * epsT * M[n+2] * G[n+2] <= PR_TOL
        n += 1
    end
    n
end

# Library orbit in T from Float64 data widened exactly; returns the n+1 points.
function pr_run(::Type{T}, K::T, x0::Float64, p0::Float64, n::Int) where {T}
    sys = StandardMap{T}(K)
    zs = Vector{SVector{2,T}}(undef, n + 1)
    zs[1] = SVector{2,T}(T(x0), T(p0))
    for k in 1:n
        zs[k+1] = forward(sys, zs[k])
    end
    zs
end

# Max-norm distance of a T point to the BigFloat reference point (caller sets precision).
pr_err(z, x, p) = max(abs(pr_big(z[1]) - x), abs(pr_big(z[2]) - p))

# Unit tangent growth S_n by the library (CompensatedSum, tangent, forward), renormalised every step.
function pr_growth_lib(::Type{T}, K::T, x0, p0, u0, n) where {T}
    sys = StandardMap{T}(K)
    z = SVector{2,T}(T(x0), T(p0)); u = SVector{2,T}(T(u0[1]), T(u0[2]))
    acc = CompensatedSum{T}()
    S = Vector{T}(undef, n + 1); S[1] = zero(T)
    for k in 1:n
        w = tangent(sys, z) * u
        nw = pr_norm(w)
        acc = add(acc, log(nw / pr_norm(u)))
        u = w / nw
        z = forward(sys, z)
        S[k+1] = total(acc)
    end
    S
end

# Same quantity in BigFloat with the reference formulas (caller sets precision).
function pr_growth_ref(K::BigFloat, x0, p0, u0, n)
    x = BigFloat(x0); p = BigFloat(p0); u = (BigFloat(u0[1]), BigFloat(u0[2]))
    S = Vector{BigFloat}(undef, n + 1); S[1] = zero(BigFloat)
    for k in 1:n
        c = K * cos(x)
        w = ((1 + c) * u[1] + u[2], c * u[1] + u[2])
        nw = pr_norm(w)
        S[k+1] = S[k] + log(nw / pr_norm(u))
        u = (w[1] / nw, w[2] / nw)
        p = p + K * sin(x); x = x + p
    end
    S
end

# Double64 with all 106 mantissa bits random: hi uniform in [a, b], low word within half an ulp of hi.
pr_d64_from(rng, h) = Double64(h) + Double64((rand(rng) - 0.5) * eps(h))
pr_rand_d64(rng, a, b) = pr_d64_from(rng, a + (b - a) * rand(rng))

# Probes for the fallback matrix (resolve_backend contract).
pr_probe_good() = (loadable = true, functional = true, device = "pr-test-gpu", compiles = true, gpu_faster = true)
pr_probe_noload() = (loadable = false, functional = false, device = "", compiles = false, gpu_faster = false)
pr_probe_nofunc() = (loadable = true, functional = false, device = "", compiles = true, gpu_faster = true)
pr_probe_nocompile() = (loadable = true, functional = true, device = "pr-test-gpu", compiles = false, gpu_faster = true)
pr_probe_slow() = (loadable = true, functional = true, device = "pr-test-gpu", compiles = true, gpu_faster = false)
pr_probe_boom() = error("probe must not be called for BigFloat on auto or gpu")

# Hot-loop kernels for allocation checks (returns the result so nothing is dead-code eliminated).
# Allocation measured behind a function barrier (testset scope boxes return values on Julia 1.10).
pr_allocs(f::F, args::Vararg{Any,N}) where {F,N} = (f(args...); @allocated f(args...))   # forces specialization on f
pr_loop_forward(sys, z, n) = (for _ in 1:n; z = forward(sys, z); end; z)
pr_loop_tangent(sys, z, n) = (m = tangent(sys, z); for _ in 1:n; m = tangent(sys, z); end; m)
pr_loop_add(acc, x, n) = (for _ in 1:n; acc = add(acc, x); end; acc)
