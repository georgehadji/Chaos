# Helpers for the henon-heiles suite (prefix hh_). Definitions only: nothing here calls a new Chaoslab name (energy starts, section map, escape tools) at load time;
# every call sits inside a function body or a testset, so a missing module errors the testsets that use it and nothing else.
# Spec: workspaces/module-build/stages/01-spec/output/henon-heiles-spec.md (D116); API and conventions: workspaces/module-build/stages/02-tests/output/henon-heiles-test-plan.md.
using Test, Chaoslab, StaticArrays, LinearAlgebra, StableRNGs, DoubleFloats, ForwardDiff, TOML

const HH_SEED = 20261012                  # _config/numerics-defaults.toml, key `seed`
const HH_BASE = HH_SEED + 42000           # this slug uses 20261012 + 42000 .. 42999; testset AC<n> uses 40 n + sub (sub < 40)
hh_seed(ac::Integer, sub::Integer = 0) = HH_BASE + 40 * ac + sub
const HH_BITS = 256                       # numerics-defaults.toml, key `bigfloat_bits`
# Package sizes (D38, 3 minute ceiling) or the full sizes of the spec: verify/henon-heiles.jl sets HH_FULL_RUN = true before including runtests.jl.
const HH_FULL = isdefined(Main, :HH_FULL_RUN) && Main.HH_FULL_RUN
hh_size(small, full) = HH_FULL ? full : small

hh_eps(::Type{T}) where {T} = Float64(effective_eps(T))
hh_big(x::Double64) = BigFloat(x.hi) + BigFloat(x.lo)     # exact; call where the precision in force is HH_BITS
hh_big(x::AbstractFloat) = BigFloat(x)
hh_big(x::Integer) = BigFloat(x)
hh_defaults() = TOML.parsefile(joinpath(@__DIR__, "..", "..", "_config", "numerics-defaults.toml"))

# ------------------------------------------------------------------ test-side physics (the spec's equations, independent of the module)
hh_u(q2) = q2^2 / 2 - q2^3 / 3                                    # Equation 7: V(0, q2)
hh_V(q1, q2) = (q1^2 + q2^2) / 2 + q1^2 * q2 - q2^3 / 3           # Equation 1
hh_H(z) = (z[3]^2 + z[4]^2) / 2 + hh_V(z[1], z[2])
hh_rad(z) = sqrt(z[1]^2 + z[2]^2)
# Equation 8 written out: p1 = sqrt(2 (E - p2^2/2 - u(q2))) >= 0 or `nothing` outside D(E)
function hh_start(E::T, q2::T, p2::T) where {T <: AbstractFloat}
    r = 2 * (E - p2^2 / 2 - hh_u(q2))
    return r >= 0 ? SVector{4, T}(zero(T), q2, sqrt(r), p2) : nothing
end
hh_reflect(z) = SVector(-z[1], z[2], -z[3], z[4])                 # q1 -> -q1 (and p1 -> -p1)
hh_flipp(z) = SVector(z[1], z[2], -z[3], -z[4])                   # time reversal
hh_slope(x, y) = (xm = sum(x) / length(x); ym = sum(y) / length(y); sum((x .- xm) .* (y .- ym)) / sum((x .- xm) .^ 2))
hh_mean(v) = sum(v) / length(v)
hh_std(v) = (m = hh_mean(v); sqrt(sum((x - m)^2 for x in v) / (length(v) - 1)))
hh_B(::Type{T}, M, n, tau) where {T} = 16 * hh_eps(T) * M * n * (1 + n * tau)     # round-off bound B(n) of integrator-catalog-spec, Acceptance criteria preamble

# bisection of f (a sign change on [lo, hi]) to full Float64 resolution
function hh_bisect(f, lo, hi)
    flo = f(lo) > 0
    for _ in 1:100
        mid = (lo + hi) / 2
        (f(mid) > 0) == flo ? (lo = mid) : (hi = mid)
    end
    return (lo + hi) / 2
end
# Equation 7: the roots a < 0 < b of u = E bounding the origin component of the section line; at E = 1/6 the closed form (-1/2, 1) (b is a double root)
hh_interval(E) = E == 1 / 6 ? (-0.5, 1.0) : (hh_bisect(q -> hh_u(q) - E, -1.0, 0.0), hh_bisect(q -> hh_u(q) - E, 0.0, 1.0))
# Equation 6: r_b(E), the root of v(r) = r^2/2 - r^3/3 = E in (0, 1]
hh_rb(E) = E == 1 / 6 ? 1.0 : hh_bisect(r -> r^2 / 2 - r^3 / 3 - E, 0.0, 1.0)
# Equation 8: A(E) = integral over [a, b] of 2 sqrt(2 (E - u)) dq2, by q2 = (a + b)/2 - (b - a)/2 cos(phi) (smooth integrand) and the midpoint rule
function hh_area(E; n = 20000)
    a, b = hh_interval(E)
    s = 0.0
    for k in 1:n
        phi = pi * (k - 0.5) / n
        q2 = (a + b) / 2 - (b - a) / 2 * cos(phi)
        s += 2 * sqrt(max(2 * (E - hh_u(q2)), 0.0)) * (b - a) / 2 * sin(phi)
    end
    return s * pi / n
end

# ------------------------------------------------------------------ the D3 group and rotations in 256 bits (call inside setprecision)
function hh_rot_big(k::Integer)                           # rotation by 2 pi k / 3, exact entries cos = (-1/2 or 1), sin = (+-sqrt(3)/2 or 0)
    s3 = sqrt(BigFloat(3)); h = BigFloat(1) / 2
    k % 3 == 0 && return BigFloat[1 0; 0 1]
    k % 3 == 1 && return BigFloat[-h -s3/2; s3/2 -h]
    return BigFloat[-h s3/2; -s3/2 -h]
end
function hh_apply(G::AbstractMatrix, z)                  # (G q, G p) in the precision of G, rounded to the element type of z (call with G a BigFloat matrix)
    T = eltype(z); zb = BigFloat.(z)
    q = G * zb[1:2]; p = G * zb[3:4]
    return SVector{4, T}(T(q[1]), T(q[2]), T(p[1]), T(p[2]))
end

# ------------------------------------------------------------------ reference orbits (package oracle O3: Double64 order-10 table at tau/4 with the sub-step root of Equation 11;
# the Taylor reference R1 of the spec is stage 04, the package environment has no TaylorIntegration)
hh_refmethod() = SofroniouSpaletta10(Double64)
const HH_REFTAU = 80                                      # reference steps per time unit (tau_ref = 1/80 = tau/4)
# root s in (0, tau] of g(Phi_s(z)) = 0 with g(z) < 0 <= g(Phi_tau(z)), by bisection (110 halvings, below the 2^-104 resolution)
function hh_ref_root(m, flow, z, tau, g)
    lo, hi = zero(tau), tau
    for _ in 1:110
        mid = (lo + hi) / 2
        g(advance(m, flow, z, mid)) < 0 ? (lo = mid) : (hi = mid)
    end
    s = (lo + hi) / 2
    return s, advance(m, flow, z, s)
end
# crossings of q1 = 0 upward (p1 > 0) of an orbit from the Float64 state z0 up to time T: (t, q2, p2) as Float64 vectors; a start on the section is the first point at t = 0
function hh_ref_section(z0, T)
    D = Double64; m = hh_refmethod(); flow = HenonHeilesFlow(); tau = D(1) / HH_REFTAU
    z = SVector{4, D}(D.(z0)); t = D(0)
    ts = Float64[]; q2s = Float64[]; p2s = Float64[]
    if z0[1] == 0 && z0[3] > 0
        push!(ts, 0.0); push!(q2s, z0[2]); push!(p2s, z0[4])
    end
    for _ in 1:round(Int, T * HH_REFTAU)
        zn = advance(m, flow, z, tau)
        if z[1] < 0 && zn[1] >= 0 && zn[3] > 0
            s, zc = hh_ref_root(m, flow, z, tau, w -> w[1])
            push!(ts, Float64(t + s)); push!(q2s, Float64(zc[2])); push!(p2s, Float64(zc[4]))
        end
        z = zn; t += tau
    end
    return (t = ts, q2 = q2s, p2 = p2s)
end
# the next upward crossing from the section point (q2, p2) at energy E (Double64 inputs): (q2', p2'), NaN if none within 200 time units
function hh_ref_return(E::D, q2::D, p2::D) where {D <: AbstractFloat}
    m = hh_refmethod(); flow = HenonHeilesFlow(); tau = D(1) / HH_REFTAU
    z = SVector{4, D}(zero(D), q2, sqrt(2 * (E - p2^2 / 2 - hh_u(q2))), p2)
    for _ in 1:(200 * HH_REFTAU)
        zn = advance(m, flow, z, tau)
        if z[1] < 0 && zn[1] >= 0 && zn[3] > 0
            _, zc = hh_ref_root(m, flow, z, tau, w -> w[1])
            return zc[2], zc[4]
        end
        z = zn
    end
    return D(NaN), D(NaN)
end
# finite-difference Jacobian d(q2', p2')/d(q2, p2) of the reference return map at the section point (q2, p2) of energy E (central, step h)
function hh_ref_jacobian(E, q2, p2; h = 1.0e-6)
    D = Double64; E = D(E); q2 = D(q2); p2 = D(p2); h = D(h)
    cq = (hh_ref_return(E, q2 + h, p2), hh_ref_return(E, q2 - h, p2))
    cp = (hh_ref_return(E, q2, p2 + h), hh_ref_return(E, q2, p2 - h))
    col(a, b) = SVector(Float64((a[1] - b[1]) / (2h)), Float64((a[2] - b[2]) / (2h)))
    return hcat(col(cq...), col(cp...))
end
# the state z_n of the orbit with q1(z_n) < 0 <= q1(z_{n+1}) and p1(z_{n+1}) > 0 (a bracketed upward crossing), or `nothing` within nmax steps
function hh_bracket(m, flow, z, tau, nmax = 100000)
    for _ in 1:nmax
        zn = advance(m, flow, z, tau)
        (z[1] < 0 && zn[1] >= 0 && zn[3] > 0) && return z
        z = zn
    end
    return nothing
end
# first time r > R of the orbit from z0 (Float64 state), located by the sub-step root; `nothing` if none before Tmax or on non-finite values
function hh_ref_escape(z0, R, Tmax)
    D = Double64; m = hh_refmethod(); flow = HenonHeilesFlow(); tau = D(1) / HH_REFTAU
    z = SVector{4, D}(D.(z0)); t = D(0); R = D(R)
    for _ in 1:round(Int, Tmax * HH_REFTAU)
        zn = advance(m, flow, z, tau)
        all(isfinite, zn) || return nothing
        if hh_rad(zn) > R
            s, _ = hh_ref_root(m, flow, z, tau, w -> hh_rad(w) - R)
            return Float64(t + s)
        end
        z = zn; t += tau
    end
    return nothing
end

# ------------------------------------------------------------------ estimators on the one-step map of the integrator (existing modules; exponents per unit time = per step / tau)
# One orbit, the QR spectrum with renormalization every `renorm` steps (time unit = 1/tau steps by default). Returns exponents per time unit and per step, G_max, M, n.
function hh_lyap(m, flow, z0, tau, T; alg = RGS(), W0 = one(SMatrix{4, 4, eltype(z0)}), renorm = round(Int, 1 / tau))
    im = IntegratorMap(m, flow, tau)
    res = propagate(im, [z0], round(Int, T / tau); k = 4, tau = renorm, seed = 1, reducers = (QRSpectrum(alg),), W0 = [W0])
    f = res.finals[1][1]
    return (lambda = f.lambda ./ tau, lstep = f.lambda, sum = f.sum, G_max = f.G_max, n = f.n, M = maximum(abs, tangent(im, z0)))
end
# finite-time exponent per time unit of one seeded deviation vector (Benettin, renormalization every time unit)
function hh_lambda1(m, flow, z0, tau, T; seed = 1)
    im = IntegratorMap(m, flow, tau)
    res = propagate(im, [z0], round(Int, T / tau); k = 1, tau = round(Int, 1 / tau), seed = seed, reducers = (Benettin(),))
    return res.finals[1][1].lambda[1] / tau
end
const HH_W2 = SMatrix{4, 2, Float64}(1, 0, 0, 0, 0, 0, 1, 0)   # the SALI frame of the printed runs: (1, 0, 0, 0) and (0, 0, 1, 0) (Equation 9)
# SALI at time T (steps tau = 1/20, renormalization every 20 steps) of the start z0
function hh_sali(m, flow, z0, T; tau = 1 / 20)
    im = IntegratorMap(m, flow, tau)
    res = propagate(im, [z0], round(Int, T / tau); k = 2, tau = round(Int, 1 / tau), seed = 1, reducers = (SALI(),), W0 = [HH_W2])
    return res.finals[1][1].sali
end
# The portrait data of Algorithm 4 at energy E: an n x n grid of starts on the section domain, SALI at T_f, class by the threshold (data only, P7)
function hh_portrait(m, E, n, Tf; thr = 1.0e-8, tau = 1 / 20, chunk = nothing)
    flow = HenonHeilesFlow(); zs = energy_grid(flow, E, n)
    im = IntegratorMap(m, flow, tau)
    res = propagate(im, zs, round(Int, Tf / tau); k = 2, tau = round(Int, 1 / tau), seed = 1, reducers = (SALI(),), W0 = fill(HH_W2, length(zs)), chunk = chunk)
    sali = [f[1].sali for f in res.finals]
    return (n = length(zs), sali = sali, chaotic = sali .< thr, fraction = count(<(thr), sali) / length(zs), x = res.x, W = res.W, starts = zs)
end

# ------------------------------------------------------------------ allocation probe (function barrier; warm up first)
function hh_alloc(f, args...)
    f(args...); f(args...)
    return @allocated f(args...)
end
