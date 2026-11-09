# Helpers for the integrator-catalog suite (prefix ic_). Definitions only: nothing here calls a new Chaoslab name (Systems flows, Integrators methods)
# at load time; every call sits inside a function body or a testset, so a missing module errors the testsets that use it and nothing else.
# Spec: workspaces/module-build/stages/01-spec/output/integrator-catalog-spec.md; API and conventions:
# workspaces/module-build/stages/02-tests/output/integrator-catalog-test-plan.md.
using Test, Chaoslab, StaticArrays, LinearAlgebra, StableRNGs, DoubleFloats, ForwardDiff

const IC_SEED = 20261012                  # _config/numerics-defaults.toml, key `seed`
const IC_BASE = IC_SEED + 41000           # this slug uses 20261012 + 41000 .. 41999; testset AC<n> uses 40 n + sub (sub < 40); FL (flows) uses n = 0, API n = 24
ic_seed(ac::Integer, sub::Integer = 0) = IC_BASE + 40 * ac + sub
const IC_BITS = 256                       # numerics-defaults.toml, key `bigfloat_bits`
# Package suite sizes (D38, 3 minute ceiling) or the full sizes of the spec: verify/integrator-catalog.jl sets IC_FULL_RUN = true before including runtests.jl.
const IC_FULL = isdefined(Main, :IC_FULL_RUN) && Main.IC_FULL_RUN
ic_size(small, full) = IC_FULL ? full : small

ic_eps(::Type{T}) where {T} = Float64(effective_eps(T))   # Float64 eps, Double64 2^-100, BigFloat eps at the precision in force

# ------------------------------------------------------------------ conversions
ic_big(x::Double64) = BigFloat(x.hi) + BigFloat(x.lo)     # exact; call where the precision in force is IC_BITS
ic_big(x::AbstractFloat) = BigFloat(x)
ic_big(x::Integer) = BigFloat(x)
ic_to(::Type{T}, v) where {T} = SVector{length(v), T}(ntuple(i -> T(v[i]), length(v)))   # exact for Float64 data into Float64, Double64, BigFloat
ic_err(z, ref) = Float64(maximum(abs(ic_big(z[i]) - ref[i]) for i in eachindex(ref)))   # max-norm error of a state against a BigFloat reference
ic_mmax(z) = max(1.0, maximum(abs, Float64.(z)))

# ------------------------------------------------------------------ flows and starting points
# kinds: :f1_1 (F1, n = 1), :f1_3 (F1, n = 3, incommensurate), :f4 (F4 pendulum, H = p^2/2 - eps cos q), :f2 (F2 Kepler, planar), :t1 (Tao test flow),
# :t2 (Henon-Heiles of Hairer2008). State layout (q_1..q_n, p_1..p_n); F2 and T2 are (q1, q2, p1, p2).
const IC_F1W = (1.0, 0.7071067811865476, 0.5773502691896257)   # frequencies <= 1; 1/sqrt(2), 1/sqrt(3) rounded to binary64 (exactly representable data)
ic_f1w(kind::Symbol) = kind === :f1_1 ? IC_F1W[1:1] : IC_F1W
ic_dim(kind::Symbol) = kind === :f1_1 ? 2 : kind === :f1_3 ? 6 : kind in (:f4, :t1) ? 2 : 4
function ic_flow(::Type{T}, kind::Symbol; eps = 1.0) where {T}
    kind === :f1_1 && return HarmonicOscillator(T[T(IC_F1W[1])])
    kind === :f1_3 && return HarmonicOscillator(T[T(w) for w in IC_F1W])
    kind === :f4 && return Pendulum(T(eps))
    kind === :f2 && return Kepler()
    kind === :t1 && return TaoNonseparable()
    kind === :t2 && return HenonHeilesFlow()
    error("unknown flow kind $kind")
end

# Kepler state (q1, q2, p1, p2): GM = 1, a = 1, eccentricity e, true anomaly f (conic formulas; derived)
function ic_kepler_state(e, f)
    p = 1 - e^2
    r = p / (1 + e * cos(f))
    s = 1 / sqrt(p)
    return [r * cos(f), r * sin(f), -s * sin(f), s * (cos(f) + e)]
end
function ic_draw(rng, kind::Symbol)                      # Float64 data; spec "Acceptance criteria" preamble
    kind === :f2 && return ic_kepler_state((0.1, 0.5, 0.9)[rand(rng, 1:3)], 2 * Float64(pi) * rand(rng))
    kind === :t1 && return 3 .* rand(rng, 2) .- 1.5
    return 4 .* rand(rng, ic_dim(kind)) .- 2
end
const IC_START = Dict(:f4 => [1.0, 0.0], :f2 => ic_kepler_state(0.5, 0.7), :f1_3 => [1.0, -0.5, 0.25, 0.3, 0.1, -0.2],
                      :f1_1 => [1.0, 0.3], :t1 => [1.0, 0.5], :t2 => [0.0, 0.3, 0.371, 0.2])   # t2: p1 is a placeholder (AC15 draws p1 from H0)

# test-side Hamiltonians (generic in the element type; independent of the module)
function ic_H(kind::Symbol, z; eps = 1.0)
    if kind === :f1_1 || kind === :f1_3
        ws = ic_f1w(kind); n = length(ws); T = eltype(z)   # w^2 in the element type of z (D112): a Float64 w^2 would carry 1e-16 relative error into a Double64 or BigFloat reference
        return sum((z[n + j]^2 + T(ws[j])^2 * z[j]^2) / 2 for j in 1:n)
    end
    kind === :f4 && return z[2]^2 / 2 - eps * cos(z[1])
    kind === :f2 && return (z[3]^2 + z[4]^2) / 2 - 1 / sqrt(z[1]^2 + z[2]^2)
    kind === :t1 && return (z[1]^2 + 1) * (z[2]^2 + 1) / 2
    kind === :t2 && return (z[3]^2 + z[4]^2) / 2 + (z[1]^2 + z[2]^2) / 2 + z[1]^2 * z[2] - z[2]^3 / 3
    error("unknown flow kind $kind")
end
ic_L(z) = z[1] * z[4] - z[2] * z[3]                      # F2 angular momentum q x p
# test-side vector field J grad H by ForwardDiff of the test-side Hamiltonian
function ic_vf(kind::Symbol, z; eps = 1.0)
    g = ForwardDiff.gradient(w -> ic_H(kind, w; eps = eps), collect(z)); h = length(z) ÷ 2
    return vcat(g[(h + 1):end], -g[1:h])
end

# exact solution of F1 (rotation in each plane), BigFloat; call inside setprecision
function ic_f1_exact(kind::Symbol, z0, t)
    ws = ic_f1w(kind); n = length(ws)
    out = Vector{BigFloat}(undef, 2n)
    for j in 1:n
        w = BigFloat(ws[j]); q = BigFloat(z0[j]); p = BigFloat(z0[n + j]); c = cos(w * t); s = sin(w * t)
        out[j] = q * c + p / w * s
        out[n + j] = -q * w * s + p * c
    end
    return out
end

# Oracle O1 (spec): F1 exact; F2, F4, T1 by self-convergence of Gauss s = 8 at step 1/64 in 256 bits (iteration cap raised to 500 and the stopping tolerance set to 1e-40 absolute: the defaults, 50 iterations and 8 eps at 256 bits, are too tight for a reference) (the Taylor reference needs the extension, AC19/AC21).
const IC_REFS = Dict{Any, Vector{BigFloat}}()
function ic_ref(kind::Symbol, z0::AbstractVector{Float64}, T::Real; eps = 1.0)
    return get!(IC_REFS, (kind, Tuple(z0), Float64(T), eps)) do
        setprecision(BigFloat, IC_BITS) do
            kind in (:f1_1, :f1_3) && return ic_f1_exact(kind, z0, BigFloat(T))
            flow = ic_flow(BigFloat, kind; eps = eps)
            return Vector(advance_n(Gauss(BigFloat, 8; maxiter = 500, tol = 1e-40), flow, ic_to(BigFloat, z0), BigFloat(1) / 64, round(Int, T * 64)))
        end
    end
end

# ------------------------------------------------------------------ method catalog (thunks: nothing is called at load time)
struct ICM
    id::String
    mk::Function          # T -> method in element type T
    p::Int                # nominal order (for P rows: the order at eps = 1, F1 linear)
    sym::Union{Bool, Nothing}   # declared palindromic; nothing = not asserted (S9)
    dtab::Float64         # digits of the tabulated constants, Inf if algebraic or rational
end
const IC_S = (
    ICM("S1DK", T -> SymplecticEuler(T, :DK), 1, false, Inf), ICM("S1KD", T -> SymplecticEuler(T, :KD), 1, false, Inf),
    ICM("S2DKD", T -> Leapfrog(T, :DKD), 2, true, Inf), ICM("S2KDK", T -> Leapfrog(T, :KDK), 2, true, Inf),
    ICM("S3", T -> McLachlan2(T), 2, true, Inf), ICM("S4", T -> Ruth3(T), 3, false, Inf),
    ICM("S5", T -> ForestRuth(T), 4, true, Inf),
    ICM("S6-4", T -> triple_jump(Leapfrog(T), 4), 4, true, Inf), ICM("S6-6", T -> triple_jump(Leapfrog(T), 6), 6, true, Inf),
    ICM("S6-8", T -> triple_jump(Leapfrog(T), 8), 8, true, Inf),
    ICM("S7-4", T -> suzuki_jump(Leapfrog(T), 4), 4, true, Inf), ICM("S7-6", T -> suzuki_jump(Leapfrog(T), 6), 6, true, Inf),
    ICM("S7-8", T -> suzuki_jump(Leapfrog(T), 8), 8, true, Inf),
    ICM("S8", T -> PEFRL(T), 4, true, 16), ICM("S9", T -> CalvoSanzSerna4(T), 4, nothing, 18),
    ICM("S11-6A", T -> Yoshida(T, 6, :A), 6, true, 15), ICM("S11-8D", T -> Yoshida(T, 8, :D), 8, true, 15),
    ICM("S11-8A", T -> Yoshida(T, 8, :A), 8, true, 15),
    ICM("S12-6", T -> KahanLi(T, 6), 6, true, 20), ICM("S12-8", T -> KahanLi(T, 8), 8, true, 20),
    ICM("S12-s7odr6", T -> KahanLi(T, :s7odr6), 6, true, 20),
    ICM("S13", T -> SofroniouSpaletta10(T), 10, true, 32),
)
ic_psaba(n) = ICM("SABA$n", T -> SABA(T, n), 2, true, 40)
ic_psbab(n) = ICM("SBAB$n", T -> SBAB(T, n), 2, true, 40)
const IC_P = (ic_psaba(1), ic_psaba(2), ic_psaba(3), ic_psaba(4), ic_psbab(1), ic_psbab(2), ic_psbab(3), ic_psbab(4),
              ICM("ABA864", T -> ABA864(T), 4, true, 40), ICM("ABAH864", T -> ABAH864(T), 4, true, 40))
const IC_PBIG = (ic_psaba(6), ic_psaba(10), ic_psbab(6), ic_psbab(10))      # BigFloat only (spec criterion 6)
const IC_PC = (ICM("SABAC2", T -> SABAC(T, 2), 4, true, 40), ICM("SABAC3", T -> SABAC(T, 3), 4, true, 40), ICM("SABAC4", T -> SABAC(T, 4), 4, true, 40),
               ICM("SBABC2", T -> SBABC(T, 2), 4, true, 40), ICM("SBABC3", T -> SBABC(T, 3), 4, true, 40), ICM("SBABC4", T -> SBABC(T, 4), 4, true, 40),
               ICM("SABAC1", T -> SABAC(T, 1), 2, true, 40), ICM("SBABC1", T -> SBABC(T, 1), 2, true, 40))
const IC_EXPL = (IC_S..., IC_P...)                       # explicit table methods (corrector methods are IC_PC)
const IC_ALLX = (IC_EXPL..., IC_PC...)                   # every explicit method with a one-step map on a separable flow
# implicit: Gauss s = 1..4 (G1 = s = 1)
const IC_GAUSS = ((1, T -> ImplicitMidpoint(T)), (2, T -> Gauss(T, 2)), (3, T -> Gauss(T, 3)), (4, T -> Gauss(T, 4)))
ic_slope_type(p::Integer) = p <= 2 ? Float64 : (p <= 8 ? Double64 : BigFloat)   # D112: orders 3 and 4 in double-double too (a Float64 scan of S4 on F2 reaches the round-off band before its slope settles)

# expected force (kick) evaluations of one step from a stored table: adjacent kicks (no drift between) share one evaluation
function ic_expected_evals(c, d)
    cnt = 0; prev = 0
    for i in eachindex(d)
        iszero(d[i]) && continue
        merged = prev > 0 && all(iszero, c[(prev + 1):i])
        merged || (cnt += 1)
        prev = i
    end
    return cnt
end

# D113: force evaluations per step in a long run, from the stored table: maximal runs of kicks in the cyclic sequence (drift c_i, then kick d_i, stage by stage)
function ic_expected_cyclic(c, d)
    ev = Symbol[]
    for i in eachindex(d)
        iszero(c[i]) || push!(ev, :D)
        iszero(d[i]) || push!(ev, :K)
    end
    :D in ev || return :K in ev ? 1 : 0
    return count(i -> ev[i] === :K && ev[mod1(i - 1, length(ev))] !== :K, eachindex(ev))
end

# ------------------------------------------------------------------ polynomial step matrices (F1, n = 1, omega = 1: q' = p, p' = -q), 256 bits
# M(tau) = prod_i K_i D_i as a matrix polynomial truncated at degree deg; D_i = I + c_i tau E12, K_i = I - d_i tau E21 (drift then kick, first stage rightmost)
function ic_steppoly(c, d, deg)
    M = [zeros(BigFloat, 2, 2) for _ in 0:deg]; M[1] = Matrix{BigFloat}(I, 2, 2)
    E12 = BigFloat[0 1; 0 0]; E21 = BigFloat[0 0; 1 0]
    for i in eachindex(c)
        for (a, X) in ((BigFloat(c[i]), E12), (-BigFloat(d[i]), E21))
            iszero(a) && continue
            new = [copy(m) for m in M]
            for j in 2:(deg + 1)
                new[j] .+= a .* (X * M[j - 1])
            end
            M = new
        end
    end
    return M
end
function ic_rotpoly(deg)
    A = BigFloat[0 1; -1 0]; P = Matrix{BigFloat}(I, 2, 2)
    out = Matrix{BigFloat}[]
    for j in 0:deg
        push!(out, P / factorial(big(j))); P = P * A
    end
    return out
end
ic_polymul(A, B, deg) = [sum((A[i] * B[j - i + 1] for i in 1:j), init = zeros(BigFloat, 2, 2)) for j in 1:(deg + 1)]
ic_polyscale(A, w) = [A[j] * w^(j - 1) for j in eachindex(A)]
ic_leapfrogpoly(deg) = ic_steppoly(BigFloat[1 // 2, 1 // 2], BigFloat[1, 0], deg)
# composition of leapfrog pieces with step fractions w (applied w[1] first); the product of the scaled base polynomials
function ic_composepoly(w, deg)
    base = ic_leapfrogpoly(deg); out = nothing
    for wj in w
        s = ic_polyscale(base, BigFloat(wj))
        out = out === nothing ? s : ic_polymul(s, out, deg)
    end
    return out
end
# step fractions of the triple jump (Yoshida1990 Eqs. (4.12) to (4.14)) and of the 5-piece Suzuki rule (KahanLi1997 Sect. 3), built from leapfrog (order 2), in 256 bits
function ic_tj_weights(order)
    w = BigFloat[1]; n = 2
    while n < order
        z1 = 1 / (2 - BigFloat(2)^(BigFloat(1) / (n + 1))); z0 = 1 - 2 * z1
        w = vcat((z1 .* w), (z0 .* w), (z1 .* w)); n += 2
    end
    return w
end
function ic_su_weights(order)
    w = BigFloat[1]; n = 2
    while n < order
        s = 1 / (4 - BigFloat(4)^(BigFloat(1) / (n + 1)))
        w = vcat((s .* w), (s .* w), ((1 - 4s) .* w), (s .* w), (s .* w)); n += 2
    end
    return w
end

# ------------------------------------------------------------------ orders and slopes
# Round-off budget B(n) (spec preamble). eps_tab = max(eps_T, 10^-dtab).
ic_B(::Type{T}, M, n, tau) where {T} = 16 * ic_eps(T) * M * n * (1 + n * tau)
ic_etab(::Type{T}, dtab) where {T} = max(ic_eps(T), 10.0^(-dtab))
# Final-state errors for tau = tau0/2^j, j = 0..J, T = Tend. Returns (taus, errs, ns).
function ic_scan(m, flow, z0, Tend, ref, tau0, ::Type{T}; J = 4) where {T}
    taus = [tau0 / 2^j for j in 0:J]
    ns = [round(Int, Tend / t) for t in taus]
    errs = [ic_err(advance_n(m, flow, ic_to(T, z0), T(taus[j + 1]), ns[j + 1]), ref) for j in 0:J]
    return taus, errs, ns
end
# largest power-of-two tau0 <= 1 whose first error is below 1e-2 (spec criterion 4)
function ic_pick_tau0(m, flow, z0, Tend, ref, ::Type{T}) where {T}
    for j in 0:8
        tau = ldexp(1.0, -j)
        err = ic_err(advance_n(m, flow, ic_to(T, z0), T(tau), round(Int, Tend / tau)), ref)
        err < 1e-2 && return tau
    end
    return ldexp(1.0, -8)
end
# local log2 slopes of consecutive halvings, only where both errors exceed `lims`
# Asymptotic range (D112): the tables of high order or large coefficients are pre-asymptotic at the first halvings (the first error below 1e-2 does not put them
# in the range of the leading term), so the order is asserted on the two finest resolved halvings, at least three halvings must be resolved; earlier slopes are reported.
ic_asym(sl, p, tol) = length(sl) >= 3 && all(abs(s - p) <= tol for s in sl[(end - 1):end])
ic_slopes(errs, lims) = [log2(errs[j] / errs[j + 1]) for j in 1:(length(errs) - 1) if errs[j] > lims[j] && errs[j + 1] > lims[j + 1]]
# single-step local error slope on a flow in BigFloat: step h against one step of Gauss s = 8 (error ~ h^17); returns log2 slopes over h = hs
function ic_local_slopes(mk, kind, z0, hs; eps = 1.0)
    return setprecision(BigFloat, IC_BITS) do
        flow = ic_flow(BigFloat, kind; eps = eps); m = mk(BigFloat); zb = ic_to(BigFloat, z0); g = Gauss(BigFloat, 8; maxiter = 500, tol = 1e-40)
        errs = Float64[]
        for h in hs
            hb = BigFloat(h)
            push!(errs, ic_err(advance(m, flow, zb, hb), Vector(advance(g, flow, zb, hb))))
        end
        return [log2(errs[j] / errs[j + 1]) for j in 1:(length(errs) - 1)], errs
    end
end

# ------------------------------------------------------------------ Tao extended space: layout (q; p; x; y), form Omega = blockdiag(J_n, J_n)
function ic_tao_omega(::Type{T}, n::Integer) where {T}
    J = [zeros(T, n, n) Matrix{T}(I, n, n); -Matrix{T}(I, n, n) zeros(T, n, n)]
    return [J zeros(T, 2n, 2n); zeros(T, 2n, 2n) J]
end

# ------------------------------------------------------------------ allocation probe (function barrier; warm up first)
function ic_alloc(f, args...)
    f(args...); f(args...)
    return @allocated f(args...)
end

# D113 (F4 of the stage 04 report): the device hook of the ensemble kernel. A fake device installed in `Foundation.DEVICE` stands for the CUDA device of
# ext/ChaoslabCUDAExt.jl (the extension defines `Integrators.ensemble_backend(::CUDADevice)`; its presence is checked in verify/integrator-catalog.jl, where CUDA loads).
struct IcFakeDevice <: Chaoslab.Foundation.AbstractDevice end
Chaoslab.Integrators.ensemble_backend(::IcFakeDevice) = Chaoslab.Integrators.ensemble_backend(Chaoslab.Foundation.NoDevice())
struct IcUnhookedDevice <: Chaoslab.Foundation.AbstractDevice end
