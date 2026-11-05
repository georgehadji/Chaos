# Helpers for the coupled-maps suite (prefix cm_). Definitions only: nothing here calls Chaoslab at load time except
# inside function bodies. Systems: Y3 Froeschle 4D map (spec item 1), Y4 chain of N coupled standard maps (item 8).
# Both use the ordering (q_1..q_n, p_1..p_n): the first d/2 components are positions, the last d/2 momenta.

const CM_SEED = 20261012                       # _config/numerics-defaults.toml, seed
const CM_BITS = 256                            # numerics-defaults.toml, bigfloat_bits
const CM_EPS_DD_FLOOR = ldexp(1.0, -104)       # nominal double-double spacing
const CM_EPS_DD = ldexp(1.0, -100)             # accepted effective double-double eps (criterion 16 pass limit), as sm_helpers

cm_eps(::Type{Float64}) = ldexp(1.0, -52)
cm_eps(::Type{Double64}) = CM_EPS_DD

# A case: kind :y3 (p = (nu, kappa, mu)) or :y4 (p = (K, eps, 0), N sites, ring or open chain)
struct CMCase
    kind::Symbol
    N::Int
    p::NTuple{3, Float64}
    ring::Bool
end

cm_y3(nu, kappa, mu) = CMCase(:y3, 2, (Float64(nu), Float64(kappa), Float64(mu)), false)
cm_y4(N, K, e; ring = false) = CMCase(:y4, N, (Float64(K), Float64(e), 0.0), ring)

cm_d(c::CMCase) = c.kind === :y3 ? 4 : 2 * c.N
cm_h(c::CMCase) = cm_d(c) ÷ 2                  # number of positions
cm_amax(c::CMCase) = c.kind === :y3 ? max(abs(c.p[1]), abs(c.p[2])) + 2 * abs(c.p[3]) : abs(c.p[1]) + 4 * abs(c.p[2])
cm_pscale(c::CMCase) = max(1.0, abs(c.p[1]), abs(c.p[2]), abs(c.p[3]))

const CM_Y3_SETS = ((0.5, 0.1, 1.0e-3), (0.0, 0.0, 0.0), (1.0, 1.0, 0.0), (2.0, 0.5, 0.1), (-1.5, 3.0, 0.5))
# fixed Y4 panel (N, K, eps, ring): includes the published eps = 2 beta values, both boundaries, signs of eps
const CM_Y4_CONFIGS = (
    (2, 3.0, 0.2, false), (2, 3.0, 0.6102, false), (3, 1.0, -0.5, true), (3, 0.5, 1.0, false),
    (5, 5.0, 0.6102, true), (5, 0.0, 1.0, false), (10, 3.0, 0.2, true), (10, 1.0, -0.5, false),
)

# systems built in the element type T (parameters are exactly representable Float64 data)
function cm_sys(::Type{T}, c::CMCase) where {T}
    c.kind === :y3 && return FroeschleMap(T(c.p[1]), T(c.p[2]), T(c.p[3]))
    return StandardMapChain(c.N, T(c.p[1]), T(c.p[2]); boundary = c.ring ? :ring : :open)
end

function cm_wsys(::Type{T}, c::CMCase) where {T}
    c.kind === :y3 && return WrappedFroeschleMap(T(c.p[1]), T(c.p[2]), T(c.p[3]))
    return WrappedStandardMapChain(c.N, T(c.p[1]), T(c.p[2]); boundary = c.ring ? :ring : :open)
end

# lifted state: SVector{4} for Y3, plain Vector for Y4 (the suite never assumes StaticArrays for N > 4, P8)
cm_state(::Type{T}, c::CMCase, z) where {T} = c.kind === :y3 ? SVector{4, T}(ntuple(i -> T(z[i]), 4)) : T[T(v) for v in z]

# wrapped state (z reduced coordinates; all counts zero): length 2 d
function cm_wstate(::Type{T}, c::CMCase, z) where {T}
    d = cm_d(c)
    c.kind === :y3 && return SVector{8, T}(ntuple(i -> i <= 4 ? T(z[i]) : zero(T), 8))
    return T[i <= d ? T(z[i]) : zero(T) for i in 1:(2d)]
end

cm_draw(rng, d::Integer, lim) = [(2 * rand(rng) - 1) * lim for _ in 1:d]
cm_draw01(rng, d::Integer) = [2 * Float64(pi) * rand(rng) for _ in 1:d]     # in [0, 2 pi)

# Seeded cases: m points for each Y3 set, m points each with own uniform parameters in [-5, 5]^3, m points for each Y4 config with N <= maxN.
function cm_cases(rng, m::Integer; maxN::Integer = 10, y3::Bool = true, y4::Bool = true, uniform::Bool = true, lim = 4 * Float64(pi), wrapped::Bool = false)
    out = Tuple{CMCase, Vector{Float64}}[]
    draw(c) = wrapped ? cm_draw01(rng, cm_d(c)) : cm_draw(rng, cm_d(c), lim)
    if y3
        for p in CM_Y3_SETS, _ in 1:m
            c = cm_y3(p...)
            push!(out, (c, draw(c)))
        end
        if uniform
            for _ in 1:m
                c = cm_y3((2 * rand(rng) - 1) * 5, (2 * rand(rng) - 1) * 5, (2 * rand(rng) - 1) * 5)
                push!(out, (c, draw(c)))
            end
        end
    end
    if y4
        for (N, K, e, ring) in CM_Y4_CONFIGS
            N <= maxN || continue
            c = cm_y4(N, K, e; ring)
            for _ in 1:m
                push!(out, (c, draw(c)))
            end
        end
    end
    return out
end

# ---- invariants: sigma, n_H, M_n, tau (spec "Invariants", error budget) ----
cm_s(a) = (a + sqrt(a^2 + 4)) / 2
cm_sigma(c::CMCase) = (1 + sqrt(5.0)) / 2 * cm_s(cm_amax(c))

function cm_nH(c::CMCase)
    s = cm_sigma(c)
    n = 0
    while n < 20 && s^(n + 1) <= 1.0e6
        n += 1
    end
    return max(n, 1)
end

# M_n = max(1, parameter scale, max |component| over zs[1:n+1])
function cm_Mn(c::CMCase, zs, n::Integer)
    m = cm_pscale(c)
    for k in 1:(n + 1), v in zs[k]
        m = max(m, abs(Float64(v)))
    end
    return m
end

cm_tau(epsT::Real, c::CMCase, zs, n::Integer) = 32 * epsT * cm_Mn(c, zs, n) * n * cm_sigma(c)^n

# ---- exact widening and distances (256-bit, returned as Float64) ----
cm_big(x::Double64) = BigFloat(x.hi) + BigFloat(x.lo)
cm_big(x::AbstractFloat) = BigFloat(x)

cm_absdiff(a, b) = setprecision(BigFloat, CM_BITS) do
    Float64(maximum(abs(cm_big(a[i]) - cm_big(b[i])) for i in 1:length(a)))
end

# max over components of |a - b| reduced to (-pi, pi]; compares the first nc components (default all)
function cm_tdist(a, b, nc::Integer = length(a))
    return setprecision(BigFloat, CM_BITS) do
        tp = 2 * BigFloat(pi)
        m = 0.0
        for i in 1:nc
            d = cm_big(a[i]) - cm_big(b[i])
            d -= tp * round(d / tp)
            m = max(m, Float64(abs(d)))
        end
        m
    end
end

# n+1 states of the orbit under f (default forward; pass f = inverse for the time-reversed orbit)
function cm_orbit(sys, z, n::Integer; f = forward)
    zs = [z]
    for _ in 1:n
        push!(zs, f(sys, zs[end]))
    end
    return zs
end

# ---- independent plain-formula references (not Chaoslab): used for orbits, tangents and the 256-bit checks ----
# neighbour lists of the chain (open: bonds (j, j+1); ring: bond (N, 1) added, N >= 3)
function cm_nbrs(N::Integer, ring::Bool)
    nb = [Int[] for _ in 1:N]
    for j in 1:(N - 1)
        push!(nb[j], j + 1)
        push!(nb[j + 1], j)
    end
    if ring && N >= 3
        push!(nb[1], N)
        push!(nb[N], 1)
    end
    return nb
end

# F_j(q) = K sin q_j + eps sum_k sin(q_j - q_k)
function cm_force(K, e, nb, q)
    return [K * sin(q[j]) + e * sum((sin(q[j] - q[k]) for k in nb[j]); init = zero(q[j])) for j in eachindex(q)]
end

# G = dF/dq (symmetric)
function cm_G(K, e, nb, q)
    N = length(q)
    G = zeros(typeof(q[1]), N, N)
    for j in 1:N
        G[j, j] = K * cos(q[j])
        for k in nb[j]
            G[j, j] += e * cos(q[j] - q[k])
            G[j, k] = -e * cos(q[j] - q[k])
        end
    end
    return G
end

# one lifted step, p = (par1, par2, par3) in the element type of z
function cm_ref_step(c::CMCase, p, z)
    if c.kind === :y3
        nu, ka, mu = p
        x1 = z[1] + z[3]
        x2 = z[2] + z[4]
        f = mu * (one(x1) - cos(x1 + x2))
        return [x1, x2, z[3] - nu * sin(x1) - f, z[4] - ka * sin(x2) - f]
    end
    N = c.N
    nb = cm_nbrs(N, c.ring)
    mom = [z[N + j] + F for (j, F) in enumerate(cm_force(p[1], p[2], nb, z[1:N]))]
    return vcat([z[j] + mom[j] for j in 1:N], mom)
end

# inverse step
function cm_ref_inverse(c::CMCase, p, z)
    if c.kind === :y3
        nu, ka, mu = p
        f = mu * (one(z[1]) - cos(z[1] + z[2]))
        x3 = z[3] + nu * sin(z[1]) + f
        x4 = z[4] + ka * sin(z[2]) + f
        return [z[1] - x3, z[2] - x4, x3, x4]
    end
    N = c.N
    nb = cm_nbrs(N, c.ring)
    q = [z[j] - z[N + j] for j in 1:N]
    return vcat(q, [z[N + j] - F for (j, F) in enumerate(cm_force(p[1], p[2], nb, q))])
end

function cm_ref_orbit(c::CMCase, p, z, n::Integer; f = cm_ref_step)
    zs = [collect(z)]
    for _ in 1:n
        push!(zs, f(c, p, zs[end]))
    end
    return zs
end

# one-step tangent, dense Matrix in the element type of z (spec items 4 and 9), Jacobian at the pre-step point z
function cm_ref_tangent(c::CMCase, p, z)
    if c.kind === :y3
        nu, ka, mu = p
        x1 = z[1] + z[3]
        x2 = z[2] + z[4]
        s = x1 + x2
        a, b, cc = -nu * cos(x1) - mu * sin(s), -mu * sin(s), -ka * cos(x2) - mu * sin(s)
        o, Z = one(a), zero(a)
        return [o Z o Z; Z o Z o; a b (o+a) b; b cc b (o+cc)]
    end
    N = c.N
    G = cm_G(p[1], p[2], cm_nbrs(N, c.ring), z[1:N])
    return [(I(N) + G) I(N); G I(N)] .* one(z[1])
end

# Jacobian of the inverse step, evaluated at the image point zi (formulas of items 5 and 9)
function cm_ref_inverse_tangent(c::CMCase, p, zi)
    if c.kind === :y3
        nu, ka, mu = p
        s = zi[1] + zi[2]
        a, b, cc = -nu * cos(zi[1]) - mu * sin(s), -mu * sin(s), -ka * cos(zi[2]) - mu * sin(s)
        o, Z = one(a), zero(a)
        return [(o+a) b -o Z; b (o+cc) Z -o; -a -b o Z; -b -cc Z o]
    end
    N = c.N
    q = [zi[j] - zi[N + j] for j in 1:N]
    G = cm_G(p[1], p[2], cm_nbrs(N, c.ring), q)
    return [I(N) -I(N); -G (I(N) + G)] .* one(zi[1])
end

cm_params(::Type{T}, c::CMCase) where {T} = (T(c.p[1]), T(c.p[2]), T(c.p[3]))

# Omega = [[0, I], [-I, 0]] in element type T
function cm_omega(::Type{T}, d::Integer) where {T}
    h = d ÷ 2
    Om = zeros(T, d, d)
    for j in 1:h
        Om[j, h + j] = one(T)
        Om[h + j, j] = -one(T)
    end
    return Om
end

# a perturbed copy of z (complex) for the complex-step derivative
cm_perturb(z::SVector, k, delta) = Base.setindex(complex.(z), complex(z[k]) + delta, k)
function cm_perturb(z::AbstractVector, k, delta)
    w = complex.(z)
    w[k] += delta
    return w
end

# Fixed-point class and magnitude from a trace (det = 1)
function cm_class(tr)
    if abs(tr) > 2
        return (:hyperbolic, (abs(tr) + sqrt(tr^2 - 4)) / 2)
    elseif abs(tr) < 2
        return (:elliptic, acos(tr / 2))
    else
        return (:parabolic, one(tr))
    end
end

# Laplacian eigenvalues of the chain graph, ascending, in Float64 from a 256-bit evaluation
function cm_laplacian_eigs(N::Integer, ring::Bool)
    return setprecision(BigFloat, CM_BITS) do
        ev = ring ? [2 - 2 * cos(2 * BigFloat(pi) * k / N) for k in 0:(N - 1)] : [2 - 2 * cos(k * BigFloat(pi) / N) for k in 0:(N - 1)]
        sort!(ev)
        Float64.(ev)
    end
end

# Cyclic Jacobi eigenvalues of a symmetric BigFloat matrix (no GenericLinearAlgebra dependency); call inside setprecision
function cm_jacobi_eigvals(A0::AbstractMatrix{BigFloat})
    A = Matrix{BigFloat}(A0)
    n = size(A, 1)
    scale = sum(abs2, A)
    for _ in 1:100
        off = sum(abs2(A[i, j]) for i in 1:n for j in 1:n if i != j)
        off <= scale * BigFloat(10)^(-100) && break
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

# singular values (descending) of a BigFloat matrix, via the eigenvalues of Y' Y
cm_svdvals(Y::AbstractMatrix{BigFloat}) = sqrt.(max.(cm_jacobi_eigvals(Y' * Y), zero(BigFloat)))

# Largest relative error of sin and cos in Double64 on the sample xs against 256-bit values (call inside setprecision)
function cm_dd_trig_error(xs)
    m = 0.0
    for x in xs
        b = cm_big(x)
        for (f, g) in ((sin, sin), (cos, cos))
            r = g(b)
            iszero(r) && continue
            m = max(m, Float64(abs(cm_big(f(x)) - r) / abs(r)))
        end
    end
    return m
end
