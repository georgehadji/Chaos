# Helpers for the continuous-qr suite (prefix cq_). Definitions only: nothing here calls a new Chaoslab name (ContinuousQR, propagate_cqr_orbit, propagate_cqr,
# cqr_rhs, cqr_growth) at load time; every call sits inside a function body or a testset, so a missing module errors the testsets that use it and nothing else.
# Spec: workspaces/module-build/stages/01-spec/output/continuous-qr-spec.md (D120); API and conventions: workspaces/module-build/stages/02-tests/output/continuous-qr-test-plan.md.
using Test, Chaoslab, StaticArrays, LinearAlgebra, StableRNGs, DoubleFloats

const CQ_SEED = 20261012                  # _config/numerics-defaults.toml, key `seed`
const CQ_BASE = CQ_SEED + 44000           # this slug uses 20261012 + 44000 .. 44999; testset AC<n> uses 50 n + sub (sub < 50)
cq_seed(ac::Integer, sub::Integer = 0) = CQ_BASE + 50 * ac + sub
const CQ_BITS = 256                       # numerics-defaults.toml, key `bigfloat_bits`
# Package sizes (D38, 3 minute ceiling for the folder) or the full sizes of the spec: verify/continuous-qr.jl sets CQ_FULL_RUN = true before including runtests.jl.
const CQ_FULL = isdefined(Main, :CQ_FULL_RUN) && Main.CQ_FULL_RUN
cq_size(small, full) = CQ_FULL ? full : small

cq_eps(::Type{T}) where {T} = Float64(effective_eps(T))
cq_big(x::Double64) = BigFloat(x.hi) + BigFloat(x.lo)     # exact; call where the precision in force is the one wanted
cq_big(x::AbstractFloat) = BigFloat(x)
cq_big(x::Integer) = BigFloat(x)
cq_mean(v) = sum(v) / length(v)
cq_var(v) = (m = cq_mean(v); sum(x -> (x - m)^2, v) / (length(v) - 1))
cq_median(v) = (s = sort(collect(v)); n = length(s); isodd(n) ? s[(n + 1) ÷ 2] : (s[n ÷ 2] + s[n ÷ 2 + 1]) / 2)
cq_slope(u, v) = sum((u .- cq_mean(u)) .* (v .- cq_mean(v))) / sum((u .- cq_mean(u)) .^ 2)
cq_agree(a, b; allowance = 0.0) = abs(cq_mean(a) - cq_mean(b)) <= 4 * sqrt(cq_var(a) / length(a) + cq_var(b) / length(b)) + allowance   # L1 S5
const CQ_EMAX = 1.0e-2                                      # calibrated health tolerance E_max (the spec proposes 1e-3, point U6)
const CQ_CACHE = Dict{Any, Any}()                          # memoized computations shared by several criteria

# ------------------------------------------------------------------ test fixtures (never in src/): constant Hamiltonian matrices as flows, the identity step map, counters
# A linear flow dz/dt = A z (constant A): the field Jacobian is A at every point, so no orbit enters (spec C19). The "orbit step" is the identity.
struct cq_LinFlow{D, T, L} <: Chaoslab.AbstractFlow
    A::SMatrix{D, D, T, L}
end
Chaoslab.dimension(::cq_LinFlow{D}) where {D} = D
Chaoslab.vector_field(f::cq_LinFlow, z::SVector) = f.A * z
Chaoslab.vector_field_jacobian(f::cq_LinFlow, ::SVector) = f.A
struct cq_FixedMap{F, T} <: Chaoslab.AbstractMap
    flow::F
    tau::T
end
Chaoslab.dimension(m::cq_FixedMap) = dimension(m.flow)
Chaoslab.forward(::cq_FixedMap, z::SVector) = z
cq_lin(A::AbstractMatrix, ::Type{T}) where {T} = (d = size(A, 1); cq_LinFlow(SMatrix{d, d, T}(T.(A))))
cq_x0(d::Integer, ::Type{T}) where {T} = zero(SVector{d, T})
# The estimator on a constant matrix: order q, step h, optional penalty factor and form.
cq_est(A::AbstractMatrix, ::Type{T}, h, q; kw...) where {T} = (f = cq_lin(A, T); ContinuousQR(cq_FixedMap(f, T(h)), f; order = q, merge((orthogonality_tol = CQ_EMAX,), NamedTuple(kw))...))
# The estimator on an integrator map with the calibrated health tolerance (spec point U6): E_max = 1e-2 (order 4) and 5e-2 (order 2) instead of the proposal 1e-3, which the
# realized equilibrium defect of k < d runs and of the order-2 scheme exceeds (plan, Calibration).
cq_cqr(im, flow; order = 4, kw...) = ContinuousQR(im, flow; order = order, merge((orthogonality_tol = order == 4 ? CQ_EMAX : 5 * CQ_EMAX,), NamedTuple(kw))...)
# Counting wrappers (Jacobian calls through the flow, orbit steps through the map), single orbits only.
struct cq_CountFlow{F} <: Chaoslab.AbstractFlow
    flow::F
    n::Base.RefValue{Int}
end
cq_CountFlow(f) = cq_CountFlow(f, Ref(0))
Chaoslab.dimension(c::cq_CountFlow) = dimension(c.flow)
Chaoslab.vector_field(c::cq_CountFlow, z::SVector) = vector_field(c.flow, z)
Chaoslab.vector_field_jacobian(c::cq_CountFlow, z::SVector) = (c.n[] += 1; vector_field_jacobian(c.flow, z))
struct cq_CountMap{M} <: Chaoslab.AbstractMap
    map::M
    tau::Float64
    n::Base.RefValue{Int}
end
cq_CountMap(m) = cq_CountMap(m, Float64(m.tau), Ref(0))
Chaoslab.dimension(c::cq_CountMap) = dimension(c.map)
Chaoslab.forward(c::cq_CountMap, z::SVector) = (c.n[] += 1; forward(c.map, z))
# A flow whose Jacobian turns NaN from the `from`-th call on (fault injection for criterion 19b), and a flow whose step map multiplies Q by (1 + eps) is not needed:
# the orthogonality fault of criterion 13c is injected through the frame (Q0 scaled), see AC13c.
struct cq_NaNJac{F} <: Chaoslab.AbstractFlow
    flow::F
    from::Int
    n::Base.RefValue{Int}
end
cq_NaNJac(f, from) = cq_NaNJac(f, from, Ref(0))
Chaoslab.dimension(c::cq_NaNJac) = dimension(c.flow)
Chaoslab.vector_field(c::cq_NaNJac, z::SVector) = vector_field(c.flow, z)
function Chaoslab.vector_field_jacobian(c::cq_NaNJac, z::SVector)
    c.n[] += 1
    A = vector_field_jacobian(c.flow, z)
    return c.n[] >= c.from ? A * oftype(first(A), NaN) : A
end

# ------------------------------------------------------------------ fixture matrices (Float64 literals; the estimator and the references use the stored values in the element type)
cq_diag2() = [0.7 0.0; 0.0 -0.7]
cq_diag4() = [0.7 0 0 0; 0 0.3 0 0; 0 0 -0.3 0; 0 0 0 -0.7]
cq_tri2() = [0.5 2.0; 0.0 -0.5]
cq_tri4() = (M = zeros(4, 4); M[1, 1] = 0.7; M[2, 2] = -0.7; M[3:4, 3:4] = cq_tri2(); M)   # direct sum of diag(0.7, -0.7) and the triangular block
cq_inv_osc(lam) = [0.0 1.0; lam^2 0.0]
cq_skew(N::Integer) = (I2 = Matrix{Float64}(I, N, N); [zeros(N, N) I2; -I2 zeros(N, N)])
# A = T^-1 A0 T with A0 a direct sum of inverted oscillators (q_j' = p_j, p_j' = lam_j^2 q_j) in the ordering (q1, q2, p1, p2) and T = [[I, B], [0, I]] (C19 (e)).
function cq_coupled(lam, B)
    A0 = [zeros(2, 2) Matrix{Float64}(I, 2, 2); Diagonal(lam .^ 2) zeros(2, 2)]
    Tm = [Matrix{Float64}(I, 2, 2) B; zeros(2, 2) Matrix{Float64}(I, 2, 2)]
    return Tm \ A0 * Tm
end
const CQ_B = [1.0 2.0; 2.0 1.0]

# ------------------------------------------------------------------ frames (orthonormal in the element type: Gram-Schmidt twice in T)
function cq_orth(M::AbstractMatrix{T}) where {T}
    M = copy(M); d, k = size(M)
    for j in 1:k, _ in 1:2
        for i in 1:(j - 1)
            M[:, j] -= dot(view(M, :, i), view(M, :, j)) * view(M, :, i)
        end
        M[:, j] ./= norm(view(M, :, j))
    end
    return SMatrix{d, k, T}(M)
end
cq_frame(rng, d::Integer, k::Integer, ::Type{T}) where {T} = cq_orth(T.(randn(rng, d, k)))
cq_eye(d::Integer, k::Integer, ::Type{T}) where {T} = SMatrix{d, k, T}(ntuple(l -> (l - 1) % d == (l - 1) ÷ d ? one(T) : zero(T), Val(d * k)))
cq_rot2(theta, ::Type{T}) where {T} = (c = cos(T(theta)); s = sin(T(theta)); SMatrix{2, 2, T}(c, s, -s, c))
cq_defect(Q::SMatrix{d, k, T}) where {d, k, T} = maximum(abs, Q' * Q - one(SMatrix{k, k, T}))   # diagnostics D4: largest absolute entry of Q^T Q - I

# ------------------------------------------------------------------ 256-bit references of the constant-matrix fixtures (C16, C19): call inside setprecision(BigFloat, CQ_BITS)
# exp(A t) of a BigFloat matrix by scaling and squaring with a Taylor series (the exact exponential of the stored matrix to well below the 256-bit rounding).
function cq_expm(A::AbstractMatrix{BigFloat})
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
# Reference l_i(t), i = 1..k, of the frame Q0 (stored values) under dz/dt = A z (C16): l_m = ln vol_m(Y X_m) - ln vol_(m-1) minus the same for X itself (X not exactly orthonormal).
function cq_lref(A::AbstractMatrix, Q0::AbstractMatrix, t)
    Ab = Matrix{BigFloat}(map(cq_big, A)); Qb = Matrix{BigFloat}(map(cq_big, Q0)); tb = cq_big(t)
    Y = cq_expm(Ab * tb); k = size(Q0, 2); out = zeros(BigFloat, k); prev = zero(BigFloat)
    for m in 1:k
        X = Qb[:, 1:m]; YX = Y * X
        lv = (log(det(YX' * YX)) - log(det(X' * X))) / 2
        out[m] = lv - prev; prev = lv
    end
    return out
end
# Closed form of C19 (b): A = diag(lam, -lam), Q0 = rotation by theta: l_1(t) and the rotation angle phi(t).
cq_l1_rot(lam, theta, t) = lam * t + log(cos(theta)^2 + sin(theta)^2 * exp(-4 * lam * t)) / 2
cq_phi_rot(lam, theta, t) = atan(exp(-2 * lam * t) * tan(theta))

# ------------------------------------------------------------------ tolerances of the spec (proposals, calibrated in the test plan)
cq_AF(A) = norm(A)                                                       # Frobenius norm of a matrix
# B_lin(t, h) = C_q ||A||_F t (1 + ||A||_F t) (Delta ||A||_F)^q + 16 eps t, C_2 = C_4 = 1 (spec preamble of the acceptance criteria)
cq_Blin(t, h, q, AF, ::Type{T}) where {T} = AF * t * (1 + AF * t) * ((q / 2) * h * AF)^q + 16 * cq_eps(T) * t

# ------------------------------------------------------------------ the scheme written test-side (the spec's Algorithm item 2 on a function of the node matrices)
# Afun(n) is the matrix at the orbit point n (n = 0 .. N). `mode` carries one defect for the negative controls: :ok, :upper (S from the strictly upper triangle),
# :Sflip (S of the wrong sign), :signu (growth with the wrong sign), :nuQ (penalty acts on Q instead of Q E), :sorted, :polar (Q replaced by its polar factor after each
# Q-step), :f32 (the product Q^T A Q formed in single precision).
function cq_G(Q::SMatrix{d, k, T}, A, nu, form::Symbol, mode::Symbol = :ok) where {d, k, T}
    AQ = A * Q
    U = mode === :f32 ? SMatrix{k, k, T}(SMatrix{k, k, Float32}(SMatrix{d, k, Float32}(Q)' * SMatrix{d, k, Float32}(AQ))) : Q' * AQ
    S = SMatrix{k, k, T}(ntuple(l -> (i = (l - 1) % k + 1; j = (l - 1) ÷ k + 1;
        mode === :upper ? (i < j ? U[i, j] : i > j ? -U[j, i] : zero(T)) : (i > j ? U[i, j] : i < j ? -U[j, i] : zero(T))), Val(k * k)))
    mode === :Sflip && (S = -S)
    Id = one(SMatrix{k, k, T}); E = Q' * Q - Id
    pen = mode === :nuQ ? nu * Q : nu * (Q * E)
    G = form === :f84 ? Q * S - pen : AQ - Q * U + (Q * S) * (Id + E) - pen
    u = SVector{k, T}(ntuple(i -> U[i, i] / (Id + E)[i, i], Val(k)))                       # D121 (M2): the rate is divided by the squared column length, as in cqr_growth
    return G, (mode === :signu ? -u : u)
end
# One Q-step of order q from the nodes: returns (Q+, increment of l, Emax over the stages' result).
function cq_step(Q::SMatrix{d, k, T}, A0, A1, A2, h, q::Integer, cnu, form, mode) where {d, k, T}
    nu(A) = T(cnu) * norm(A)
    if q == 2
        K1, u1 = cq_G(Q, A0, nu(A0), form, mode); Qs = Q + h * K1
        K2, u2 = cq_G(Qs, A1, nu(A1), form, mode)
        return Q + (h / 2) * (K1 + K2), (h / 2) * (u1 + u2)
    end
    D = 2h
    K1, u1 = cq_G(Q, A0, nu(A0), form, mode)
    K2, u2 = cq_G(Q + h * K1, A1, nu(A1), form, mode)
    K3, u3 = cq_G(Q + h * K2, A1, nu(A1), form, mode)
    K4, u4 = cq_G(Q + D * K3, A2, nu(A2), form, mode)
    return Q + (D / 6) * (K1 + 2K2 + 2K3 + K4), (D / 6) * (u1 + 2u2 + 2u3 + u4)
end
# The series at the checkpoints (orbit steps, multiples of q/2): returns (Q, l, E) at the end and `cps` a Dict step => (Q, l). l is accumulated naively in T (for BigFloat
# and for the small fixtures the rounding is far below the tolerances that use it).
function cq_scheme(Afun, Q0::SMatrix{d, k, T}, h, q::Integer, N::Integer; cnu = k == d ? 1.0 : 2.0, form::Symbol = k == d ? :f84 : :f92, mode::Symbol = :ok,
                   checkpoints = ()) where {d, k, T}
    m = q ÷ 2; Q = Q0; l = zeros(T, k); cps = Dict{Int, Any}(); Emax = cq_defect(Q0)
    A = Afun(0); n = 0
    while n < N
        A1 = Afun(n + 1); A2 = q == 4 ? Afun(n + 2) : A1
        Q, inc = cq_step(Q, A, A1, A2, T(h), q, cnu, form, mode)
        mode === :polar && (F = svd(Matrix(Q)); Q = SMatrix{d, k, T}(F.U * F.Vt))
        l .+= inc; A = q == 4 ? A2 : A1; n += m
        Emax = max(Emax, cq_defect(Q))
        n in checkpoints && (cps[n] = (Q = Q, l = copy(l)))
    end
    mode === :sorted && sort!(l; rev = true)
    return (Q = Q, l = l, Emax = Emax, cps = cps)
end

# ------------------------------------------------------------------ Henon-Heiles fixtures (spec: E = 1/8, p1 rebuilt by the module's energy_start; orbits R, C of henon-heiles Equation 9)
const CQ_E = 1 / 8
const CQ_W2 = SMatrix{4, 2, Float64}(1, 0, 0, 0, 0, 0, 1, 0)           # SALI frame of the printed runs (henon-heiles Equation 9)
cq_start(q2, p2; E = CQ_E) = energy_start(HenonHeilesFlow(), E, q2, p2)[1]
cq_orbit_R() = cq_start(0.1, 0.0)                                       # quasi-periodic orbit R
cq_orbit_C() = cq_start(-0.015, 0.0)                                    # sticky chaotic orbit C
cq_table(::Type{T}, key::Symbol) where {T} = key === :o8 ? KahanLi(T, 8) : key === :o4 ? ForestRuth(T) : key === :o2 ? Leapfrog(T) : error("unknown table $key")
const CQ_H = 1 / 20
cq_im(flow, ::Type{T}, key::Symbol = :o8, h = CQ_H) where {T} = IntegratorMap(cq_table(T, key), flow, T(h))
cq_hh(::Type{T} = Float64, key::Symbol = :o8, h = CQ_H) where {T} = (flow = HenonHeilesFlow(); (flow, cq_im(flow, T, key, h)))
# The orbit points x_0 .. x_N (iterating `forward` is bitwise the stream's orbit, L1 S2).
function cq_pts(im, x0, N::Integer)
    pts = Vector{typeof(x0)}(undef, N + 1); pts[1] = x0; z = x0
    for i in 1:N
        z = forward(im, z); pts[i + 1] = z
    end
    return pts
end
# SALI at time T (h = 1/20): the class label of henon-heiles criterion 10 (chaotic when below 1e-8).
function cq_sali(z0, T)
    flow, im = cq_hh(); n = round(Int, T / CQ_H)
    res = propagate(im, [z0], n; k = 2, tau = 20, seed = 1, reducers = (SALI(),), W0 = [CQ_W2])
    return res.finals[1][1].sali
end
# Starts of the two classes at E = 1/8 drawn by random_energy_start (seeded), classified by SALI at Tclass (package 1000, full 2000): the first `m` of each, memoized.
const CQ_NCLASS = cq_size(4, 200)
function cq_classes(m::Integer, seedoff::Integer)
    return get!(CQ_CACHE, (:classes, m, seedoff)) do
        flow = HenonHeilesFlow(); a, b = section_interval(flow, CQ_E); w = sqrt(2 * CQ_E)
        rng = StableRNG(cq_seed(0, seedoff)); Tc = cq_size(1000.0, 2000.0)
        reg = SVector{4, Float64}[]; cha = SVector{4, Float64}[]; tries = 0
        while (length(reg) < m || length(cha) < m) && tries < 50 * m + 200
            tries += 1
            z, st = random_energy_start(flow, CQ_E, (a, b), (-w, w), rng)
            st === :ok || continue
            s = cq_sali(z, Tc)
            (s < 1.0e-8 ? (length(cha) < m && push!(cha, z)) : (length(reg) < m && push!(reg, z)))
        end
        (regular = reg, chaotic = cha)
    end
end
# The discrete QR spectrum (k = 4, tau = 1) and Benettin (k = 1) on the identical orbit and frame, per unit time (per step / h).
function cq_discrete(im, z0, W0, N::Integer, h; alg = RGS())
    res = propagate(im, [z0], N; k = size(W0, 2), tau = 1, seed = 1, reducers = (QRSpectrum(alg),), W0 = [W0])
    return res.finals[1][1].lambda ./ h
end
function cq_benettin(im, z0, w0, N::Integer, h)
    res = propagate(im, [z0], N; k = 1, tau = 1, seed = 1, reducers = (Benettin(),), W0 = [SMatrix{4, 1, Float64}(w0)])
    return res.finals[1][1].lambda[1] / h
end
# Time series of the module on one orbit: runs to N steps, final record (the series is the prefix of the same run, so a run to a smaller N is bitwise a prefix).
function cq_run(cq, x0, Q0, N::Integer; every = 0)
    return every > 0 ? propagate_cqr_series(cq, x0, Q0, N; every = every) : (propagate_cqr_orbit(cq, x0, Q0, N)..., nothing)
end
# Health error of a call, or nothing.
function cq_health_error(f)
    try
        f()
        return nothing
    catch e
        e isa NumericalHealthError ? e : rethrow()
    end
end
# Bytes allocated by one call after a warm-up call, measured inside a function.
function cq_alloc_bytes(f, args...)
    f(args...)
    return @allocated f(args...)
end
