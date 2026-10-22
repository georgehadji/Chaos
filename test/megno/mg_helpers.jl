# Helpers for the megno suite (prefix mg_). Definitions only: nothing here calls an Indicators name at load time, so a missing
# Chaoslab.Indicators errors the testsets that use it and nothing else. Spec: workspaces/module-build/stages/01-spec/output/megno-spec.md
# API proposal: workspaces/module-build/stages/02-tests/output/megno-test-plan.md (top section).
# lb_regular (the K = 0.5 regular-class screening, L1 S7), lb_Psi / lb_tangents (Psi_n, L1 Amendment D35) are the single definitions in the lyap-benettin
# folder, included once behind a guard as the sali suite does: this folder sorts after lyap-benettin, so under Pkg.test they are already loaded
# (unlike fli-ofli, which sorts before it and keeps copies), and verify/megno.jl includes this suite alone.
using Test, Chaoslab, StaticArrays, LinearAlgebra, StableRNGs, DoubleFloats
@isdefined(lb_regular) || include(joinpath(@__DIR__, "..", "lyap-benettin", "lb_helpers.jl"))

const MG_SEED = 20261012                  # base seed, _config/numerics-defaults.toml key `seed`
const MG_BASE = MG_SEED + 35000           # this slug uses 20261012 + 35000 .. 35999; testset AC<n> uses 50 n + sub (sub < 50)
mg_seed(ac::Integer, sub::Integer = 0) = MG_BASE + 50 * ac + sub
const MG_TYPES = (Float64, Double64)      # Float32 is out (D24); BigFloat only for oracle runs
const MG_BITS = 256                       # numerics-defaults.toml key `bigfloat_bits`
# Package suite sizes (D38, 3 minute ceiling) or the full sizes of the spec: verify/megno.jl defines MG_FULL_RUN = true before including
# runtests.jl. Every reduction is marked "D38" at its use. The testsets are top-level code (not functions), as in the other suites.
const MG_FULL = isdefined(Main, :MG_FULL_RUN) && Main.MG_FULL_RUN

mg_eps(::Type{T}) where {T} = Float64(effective_eps(T))     # Float64 eps, Double64 2^-100, BigFloat eps at the current precision
const MG_CACHE = Dict{Any, Any}()         # memoized computations shared by several criteria

# The convention tag of the spec (Algorithm, Representation), recorded in the run record.
const MG_CONVENTION = "weights at the end of the window, mean over the handoffs 1..K"

# ------------------------------------------------------------------ BigFloat conversions

mg_big(x::Double64) = BigFloat(x.hi) + BigFloat(x.lo)       # exact; call inside setprecision(BigFloat, bits)
mg_big(x::AbstractFloat) = BigFloat(x)
mg_err(value, ref; bits = MG_BITS) = setprecision(BigFloat, bits) do
    Float64(abs(mg_big(value) - ref))
end
mg_diff(value, ref; bits = MG_BITS) = setprecision(BigFloat, bits) do
    Float64(mg_big(value) - ref)
end

# ------------------------------------------------------------------ the panel (S7) and the cases

mg_K(::Type{T}, label::Symbol) where {T} =
    label === :two_pi ? exact(T, pi) * 2 : label === :half ? exact(T, 0.5) : label === :ten ? exact(T, 10) :
    label === :five ? exact(T, 5) : error("unknown K $label")

mg_starts(rng, m::Integer) = [SVector(2pi * rand(rng), 2pi * rand(rng)) for _ in 1:m]    # Float64, uniform in [0, 2 pi)^2 (S7)
mg_cast(::Type{T}, z::SVector{2}) where {T} = SVector{2, T}(T(z[1]), T(z[2]))
mg_frame(seed, j, ::Type{T}) where {T} = initial_deviations(seed, j, Val(2), Val(2), T)   # the pair: first two seeded columns

# m starts of the panel class `label` (:half = K = 0.5 regular class; :two_pi; :ten; :five) widened exactly into T, with their orbit ids and seeded frames.
function mg_case(label::Symbol, ::Type{T}, m::Int) where {T}
    return get!(MG_CACHE, (:case, label, T, m)) do
        K = mg_K(T, label)
        if label === :half
            r = lb_regular()                                  # ids, starts and the frame seed of the L1 screening (one definition)
            ids, starts, fseed = r.ids[1:m], r.starts[1:m], r.seed
        else
            fseed = mg_seed(0, findfirst(==(label), (:two_pi, :ten, :five)))
            ids = collect(1:m)
            starts = mg_starts(StableRNG(fseed), m)
        end
        x0 = [mg_cast(T, s) for s in starts]
        (sys = StandardMap(K), K = K, ids = ids, x0 = x0, fseed = fseed, W0 = [mg_frame(fseed, id, T) for id in ids])
    end
end

# Checkpoints n <= N: multiples of 10 at 10 per decade (D38 reduction of "every handoff that is a multiple of 10"); all multiples with MG_FULL.
mg_checkpoints(N::Int) = MG_FULL ? collect(10:10:N) :
                         sort(unique(filter(n -> n % 10 == 0 && n <= N, vcat((d .* (1:10) for d in (1, 10, 100, 1000))...))))

# Panel of the same-orbit criteria 6 and 10 (class, N, starts): 20 starts per class and type, n up to 1e4, as the spec says.
mg_panel() = ((:half, 10_000, 20), (:two_pi, 10_000, 20), (:ten, 10_000, 20))

# Orbit of the standard map in T: pts[i] = x_{i-1} (the point before step i, where M_i is evaluated).
function mg_pts(sys, x0, n::Integer)
    pts = Vector{typeof(x0)}(undef, n)
    z = x0
    for i in 1:n
        pts[i] = z
        z = forward(sys, z)
    end
    return pts
end

# ------------------------------------------------------------------ BigFloat reference: the exact sums s_i(k) with s_i(0) = 0 (E2, E3, E8)

# Step-i tangent of the standard map at the realized point pts[i]: [[1 + c, 1], [c, 1]], c = K cos x (L1 E1), as (a, b, c, d).
function mg_mats_std(Kt, pts)
    return function (i)
        c = mg_big(Kt) * cos(mg_big(pts[i][1]))
        return (1 + c, one(c), c, one(c))
    end
end
mg_mats_const(a, b, c, d) = i -> (BigFloat(a), BigFloat(b), BigFloat(c), BigFloat(d))   # BigFloat at the precision in force when called

# s_i(k) = ln(norm(Y_k w_i(0)) / norm(w_i(0))) for k = 1..N, i = 1, 2, from the exact widened frame W0 and the tangent maps mats(i), at `bits` bits.
# Dividing by the norm of the stored frame column makes the convention s(0) = 0 exact (the stored frame is unit only to roundoff).
function mg_ref_s(mats, W0, N::Int; bits = MG_BITS)
    return setprecision(BigFloat, bits) do
        w1 = (mg_big(W0[1, 1]), mg_big(W0[2, 1]))
        w2 = (mg_big(W0[1, 2]), mg_big(W0[2, 2]))
        l1, l2 = log(w1[1]^2 + w1[2]^2) / 2, log(w2[1]^2 + w2[2]^2) / 2
        s1, s2 = Vector{BigFloat}(undef, N), Vector{BigFloat}(undef, N)
        for i in 1:N
            ma, mb, mc, md = mats(i)
            w1 = (ma * w1[1] + mb * w1[2], mc * w1[1] + md * w1[2])
            w2 = (ma * w2[1] + mb * w2[2], mc * w2[1] + md * w2[2])
            s1[i] = log(w1[1]^2 + w1[2]^2) / 2 - l1
            s2[i] = log(w2[1]^2 + w2[2]^2) / 2 - l2
        end
        (s1 = s1, s2 = s2)
    end
end

# ------------------------------------------------------------------ the three series from the sums (M5, M6), BigFloat

# M5, the weighted-sum definition, on the handoff sums sh[k] = s(k tau), k = 1..K (s(0) = 0): a_k = sh[k] - sh[k-1],
#   Y(t_K) = (2 / t_K) sum t_k a_k, Ybar = (1/K) sum_{j<=K} Y(t_j), Ybar20 = tau sum_{j<=K} 3 sum_{k<=j} t_k^2 a_k / t_K^3, t_k = k tau.
# Keyword switches build the faulty doubles of the negative controls (criterion 1): `noweight` (weight 1 instead of t_k), `nofactor` (no factor 2),
# `mean0` (the mean over the handoffs 0..K-1 instead of 1..K).
function mg_M5(sh::AbstractVector{BigFloat}, tau::Integer; noweight::Bool = false, nofactor::Bool = false, mean0::Bool = false, bits = MG_BITS)
    return setprecision(BigFloat, bits) do
        K = length(sh)
        Y, Yb, Y20b = Vector{BigFloat}(undef, K), Vector{BigFloat}(undef, K), Vector{BigFloat}(undef, K)
        R1 = R2 = Q = Q2 = zero(BigFloat)
        prev = zero(BigFloat)
        for k in 1:K
            t = BigFloat(k) * tau
            a = sh[k] - prev
            prev = sh[k]
            R1 += (noweight ? one(t) : t) * a
            R2 += t^2 * a
            Y[k] = (nofactor ? one(t) : 2) * R1 / t
            Qold = Q
            Q += Y[k]
            Q2 += 3 * R2
            Yb[k] = (mean0 ? Qold : Q) / k
            Y20b[k] = tau * Q2 / t^3
        end
        (Y = Y, Ybar = Yb, Ybar20 = Y20b)
    end
end

# M6, the Abel form, cumulative (all K at once): Y(n) = 2 (s(K) - (1/K) sum_{k<K} s(k)), Y20(t_K) = 3 tau^2 (K^2 s(K) - sum_{1<=k<K} (2k + 1) s(k)),
# Ybar = (1/K) sum_j Y(j), Ybar20 = tau sum_j Y20(j) / n^3. This is the BigFloat post-processing of a stored series s (criterion 7(b)).
function mg_M6(sh::AbstractVector{BigFloat}, tau::Integer; bits = MG_BITS)
    return setprecision(BigFloat, bits) do
        K = length(sh)
        Y, Yb, Y20b = Vector{BigFloat}(undef, K), Vector{BigFloat}(undef, K), Vector{BigFloat}(undef, K)
        P = W = Q = Q2 = zero(BigFloat)               # P = sum_{k<K} s(k), W = sum_{1<=k<K} (2k + 1) s(k)
        for k in 1:K
            sk = sh[k]
            Y[k] = 2 * (sk - P / k)
            y20 = 3 * BigFloat(tau)^2 * (BigFloat(k)^2 * sk - W)
            Q += Y[k]
            Q2 += y20
            Yb[k] = Q / k
            Y20b[k] = tau * Q2 / (BigFloat(k) * tau)^3
            P += sk
            W += (2k + 1) * sk
        end
        (Y = Y, Ybar = Yb, Ybar20 = Y20b)
    end
end

# The harmonic-weight form of Ybar (M6): (2/K) sum_{k=1..K} (1 - (H_K - H_k)) s(k), single K, independent of the recurrences.
function mg_Ybar_harmonic(sh::AbstractVector{BigFloat}, K::Int)
    return setprecision(BigFloat, MG_BITS) do
        H = zeros(BigFloat, K + 1)                    # H[j + 1] = H_j
        for j in 1:K
            H[j + 1] = H[j] + inv(BigFloat(j))
        end
        2 * sum((1 - (H[K + 1] - H[k + 1])) * sh[k] for k in 1:K) / K
    end
end

# Exact forms for constant increments lambda per map step (M8(a)), at the handoff K with step tau, n = K tau.
function mg_const_forms(lam::BigFloat, tau::Integer, K::Integer)
    n = K * tau
    return (Y = lam * (n + tau), Ybar = lam * tau * (K + 3) / 2, Ybar20 = lam * tau * (K + 1)^2 * (K + 2) / (4 * BigFloat(K)^2))
end

# The handoff sums sh[k] = s[k tau] from a step-indexed vector s[1..N].
mg_handoffs(s::AbstractVector, tau::Integer) = s[tau:tau:length(s)]

# ------------------------------------------------------------------ bounds (M11, L1 S3 and Amendment D35)

# Squared Frobenius norm of the first-window product (L1 S3) for a constant tangent M.
mg_Gfirst(M::SMatrix, tau::Int) = sum(abs2, M^tau)
# 6 B_S(n) with B_S = 8 eps (n + G_first) of L1 S3 (constant-tangent cases, K_s = 0; M11: the coefficient sum of the three series is below 6);
# with `psi`, max(6 B_S(n), 48 eps Psi_n^(i)), the per-vector form for an ill-conditioned frame (fli Amendment D72 item 2; test plan, Deviations).
function mg_lin_B(epsT, M::SMatrix, Ms, w0, n::Int, tau::Int, memo, psi::Bool = false)
    b = 6 * 8 * epsT * (n + mg_Gfirst(M, tau))
    psi || return b
    p = get!(memo, n) do
        lb_Psi(view(Ms, 1:n), w0)
    end
    return max(b, 48 * epsT * p)
end

# One orbit of a same-orbit reference: the Float64 tangents (for Psi), the frame, the memo tables, the BigFloat reference rows per tau.
struct mg_Orbit{M}
    Ms::Vector{SMatrix{2, 2, Float64, 4}}
    W0::M
    psi::Dict{Tuple{Int, Int}, Float64}
    ref::Dict{Int, Vector{NamedTuple}}         # tau => rows at the checkpoints that are handoffs of tau
end

# Psi_j^(i) (L1 Amendment D35, lb_Psi) of column i at step j.
mg_Psi(o::mg_Orbit, i::Int, j::Int) = get!(o.psi, (i, j)) do
    lb_Psi(view(o.Ms, 1:j), o.W0[:, i])
end
# Psi_max^(i)(n): the maximum of Psi_g^(i) over the grid g (the checkpoints), g <= n (M11, P8: a grid maximum).
mg_Psimax(o::mg_Orbit, i::Int, n::Int, grid) = maximum(mg_Psi(o, i, g) for g in grid if g <= n)
# B_M(n) = 48 eps Psi_max (M11).
mg_BM(epsT, o::mg_Orbit, i::Int, n::Int, grid) = 48 * epsT * mg_Psimax(o, i, n, grid)

# Reference rows of one start: BigFloat(256) (or `bits`) exact sums at the realized points of the T orbit, then M5 per tau, kept at the checkpoints.
function mg_ref_rows(K, pts, W0, N::Int, taus, cps; bits = MG_BITS)
    sv = mg_ref_s(mg_mats_std(K, pts), W0, N; bits = bits)
    out = Dict{Int, Vector{NamedTuple}}()
    for tau in taus
        sh1, sh2 = mg_handoffs(sv.s1, tau), mg_handoffs(sv.s2, tau)
        r1, r2 = mg_M5(sh1, tau; bits = bits), mg_M5(sh2, tau; bits = bits)
        out[tau] = [(n = n, s = (sh1[n ÷ tau], sh2[n ÷ tau]), Y = (r1.Y[n ÷ tau], r2.Y[n ÷ tau]),
                     Ybar = (r1.Ybar[n ÷ tau], r2.Ybar[n ÷ tau]), Ybar20 = (r1.Ybar20[n ÷ tau], r2.Ybar20[n ÷ tau]))
                    for n in filter(n -> n % tau == 0, cps)]
    end
    return out
end

# Taus of the reference rows: (1, 4) for the long panel of criterion 6, (1, 2, 4, 5, 10) for the tau study of criterion 9 (N = 1e3).
mg_ref_taus(N::Int) = N > 1000 ? (1, 4) : (1, 2, 4, 5, 10)

# Same-orbit reference for m starts of a panel class (standard map): BigFloat products at the realized points of the T orbit.
function mg_orbit_ref(label::Symbol, ::Type{T}, m::Int, N::Int; bits = MG_BITS) where {T}
    return get!(MG_CACHE, (:ref, label, T, m, N, bits)) do
        c = mg_case(label, T, m)
        cps = mg_checkpoints(N)
        per = map(1:m) do j
            pts = mg_pts(c.sys, c.x0[j], N)
            ref = mg_ref_rows(c.K, pts, c.W0[j], N, mg_ref_taus(N), cps; bits = bits)
            mg_Orbit(lb_tangents(c.sys, c.x0[j], N), c.W0[j], Dict{Tuple{Int, Int}, Float64}(), ref)
        end
        (cps = cps, per = per, case = c, bits = bits)
    end
end

# ------------------------------------------------------------------ reducers used by the tests (protocol of tangent-core)

# Wraps a reducer: forwards init, step, failed, error to `inner`; after the inner step it records reducer_final(inner) at the handoffs in `steps` (or
# at every handoff when `all`). The stream has no emission channel (L1 Q9), so this is how a series is read. CPU only (the state holds a Vector):
# use it with propagate_orbit.
struct mg_Probe{R}
    inner::R
    steps::Vector{Int}
    all::Bool
end
mg_Probe(inner, steps::AbstractVector{<:Integer} = Int[]; all::Bool = false) = mg_Probe(inner, collect(Int, steps), all)

Chaoslab.reducer_init(r::mg_Probe, ::Type{T}, v::Val, k::Val) where {T} = (Chaoslab.reducer_init(r.inner, T, v, k), 1, Any[])
function Chaoslab.reducer_step(r::mg_Probe, st, n::Int, x, W)
    inner, pos, snaps = st
    inner2, W2 = Chaoslab.reducer_step(r.inner, inner, n, x, W)
    if r.all || (pos <= length(r.steps) && n == r.steps[pos])
        push!(snaps, Chaoslab.reducer_final(r.inner, inner2))
        pos += 1
    end
    return (inner2, pos, snaps), W2
end
Chaoslab.reducer_failed(r::mg_Probe, st) = Chaoslab.reducer_failed(r.inner, st[1])
Chaoslab.reducer_error(r::mg_Probe, st, orbit, backend) = Chaoslab.reducer_error(r.inner, st[1], orbit, backend)
Chaoslab.reducer_final(r::mg_Probe, st) = (snaps = st[3], inner = Chaoslab.reducer_final(r.inner, st[1]))

# Read-only reducer recording the orbit point at every handoff.
struct mg_Trace end
Chaoslab.reducer_init(::mg_Trace, ::Type, ::Val, ::Val) = Any[]
function Chaoslab.reducer_update(::mg_Trace, st, n::Int, x, W)
    push!(st, x)
    return st, W
end
Chaoslab.reducer_final(::mg_Trace, st) = (xs = st,)

# Serial stand-in for `propagate` on the CPU: the per-orbit body (propagate_orbit) over orbits 1:length(x0s), the chain `rs` with every reducer wrapped in
# a probe; returns, per orbit, the tuple of reducer_final(probe) (fields snaps, inner), one entry per reducer of the chain. W0s[j] is the d x k frame of
# orbit j; `all` is a Bool or a tuple of Bools (one per reducer).
function mg_run(rs::Tuple, sys, x0s, W0s, N::Integer, tau::Integer; steps = Int[], all = false)
    alls = all isa Tuple ? all : ntuple(_ -> all, length(rs))
    return map(eachindex(x0s)) do j
        probes = map((r, a) -> mg_Probe(r, steps; all = a), rs, alls)
        _, _, st = propagate_orbit(sys, x0s[j], W0s[j], N, tau, probes; orbit = j)
        map(Chaoslab.reducer_final, probes, st)
    end
end

# Memoized MEGNO runs for mg_orbit_ref's starts: snapshots at the checkpoints that are handoffs of tau.
function mg_orbit_runs(label::Symbol, ::Type{T}, tau::Int, m::Int, N::Int) where {T}
    return get!(MG_CACHE, (:run, label, T, tau, m, N)) do
        c = mg_case(label, T, m)
        cps = filter(n -> n % tau == 0, mg_checkpoints(N))
        [r[1] for r in mg_run((MEGNO(),), c.sys, c.x0, c.W0, N, tau; steps = cps)]
    end
end

# One row per (start, checkpoint) of a standard-map case: the run, the BigFloat reference, B_M of M11 for both columns and the ratios error / bound
# (0 when both are 0) per series. Errors are formed in BigFloat at `bits` bits against the reference. grid = the Psi_max grid (the checkpoints).
function mg_rows(label::Symbol, ::Type{T}, tau::Int, m::Int, N::Int; bits = MG_BITS) where {T}
    return get!(MG_CACHE, (:rows, label, T, tau, m, N, bits)) do
        O = mg_orbit_ref(label, T, m, N; bits = bits)
        R = mg_orbit_runs(label, T, tau, m, N)
        epsT = mg_eps(T)
        rows = NamedTuple[]
        for j in 1:m
            o = O.per[j]
            refs = o.ref[tau]
            for (s, r) in zip(R[j].snaps, refs)
                s.n == r.n || error("mg_rows: run and reference disagree on the step ($(s.n) vs $(r.n))")
                B = (mg_BM(epsT, o, 1, s.n, O.cps), mg_BM(epsT, o, 2, s.n, O.cps))
                e = (Y = ntuple(i -> mg_err(s.Y[i], r.Y[i]; bits = bits), 2), Ybar = ntuple(i -> mg_err(s.Ybar[i], r.Ybar[i]; bits = bits), 2),
                     Ybar20 = ntuple(i -> mg_err(s.Ybar20[i], r.Ybar20[i]; bits = bits), 2),
                     S = ntuple(i -> mg_err(s.sum[i] + s.comp[i], r.s[i]; bits = bits), 2))
                ratio(x, b) = b == 0 ? (x == 0 ? 0.0 : Inf) : x / b
                worst = maximum(ratio(e[f][i], B[i]) for f in (:Y, :Ybar, :Ybar20), i in 1:2)
                push!(rows, (j = j, n = s.n, Y = s.Y, Ybar = s.Ybar, Ybar20 = s.Ybar20, S = (s.sum[1] + s.comp[1], s.sum[2] + s.comp[2]),
                             e = e, B = B, ratio = worst, s = s, r = r))
            end
        end
        rows
    end
end

# ------------------------------------------------------------------ test maps (never in src/)

# Fault map on the state (n, y, flag), dimension 2: the counter advances by one per step; only flagged orbits misbehave.
# kind :zero: zero tangent (both columns collapse at the first handoff); :col2: tangent diag(1, 0) (column 2 collapses when the frame is I);
# :nan: NaN in the tangent at counter value `at`.
struct mg_Fault <: Chaoslab.AbstractMap
    kind::Symbol
    at::Int
end
Chaoslab.dimension(::mg_Fault) = 2
Chaoslab.forward(::mg_Fault, z) = SVector(z[1] + 1, z[2], z[3])
function Chaoslab.tangent(m::mg_Fault, z::SVector{3, T}) where {T}
    z[3] == 1 || return one(SMatrix{2, 2, T})
    m.kind === :zero && return zero(SMatrix{2, 2, T})
    m.kind === :col2 && return SMatrix{2, 2, T}(1, 0, 0, 0)
    return SMatrix{2, 2, T}(z[1] == m.at ? T(NaN) : one(T), 0, 0, 1)
end

# Map whose forward and tangent throw: proves a rejection happens before any work.
struct mg_Boom <: Chaoslab.AbstractMap end
Chaoslab.dimension(::mg_Boom) = 2
Chaoslab.forward(::mg_Boom, z) = error("mg_Boom: forward was called")
Chaoslab.tangent(::mg_Boom, z) = error("mg_Boom: tangent was called")

# Map with the prescribed diagonal tangent diag(a_m, 1/a_m), a_m = exp(c_m), at step m (criterion 1): the state is (counter, 0), dimension 2. Column 1 of the
# identity frame carries +c_m, column 2 carries -c_m. kinds: :const c = 0.05; :log1 c = ln((m + 1)/m), so that s(k) = ln(1 + k); :alt +0.2, -0.2 alternating
# (m odd: +0.2); :pathA c = 0.05; :pathB c = 0.1 for m <= 20 and 0 afterwards (the paths of M7: s(40) = 2 for both); :long c = 0.1 (criterion 1(v)).
mg_c(kind::Symbol, m::Integer) =
    kind === :const ? 0.05 : kind === :log1 ? log((m + 1) / m) : kind === :alt ? (isodd(m) ? 0.2 : -0.2) :
    kind === :pathA ? 0.05 : kind === :pathB ? (m <= 20 ? 0.1 : 0.0) : kind === :long ? 0.1 : error("unknown kind $kind")
struct mg_Diag <: Chaoslab.AbstractMap
    kind::Symbol
end
Chaoslab.dimension(::mg_Diag) = 2
Chaoslab.forward(::mg_Diag, z) = SVector(z[1] + 1, z[2])
function Chaoslab.tangent(m::mg_Diag, z::SVector{2, T}) where {T}
    j = round(Int, Float64(z[1])) + 1
    a = m.kind === :log1 ? T(j + 1) / T(j) : exp(T(mg_c(m.kind, j)))      # :log1 exactly (m + 1)/m in T, so that s(k) = ln(1 + k) to rounding
    return SMatrix{2, 2, T}(a, 0, 0, 1 / a)
end
# G_first of L1 S3 for the diagonal map: the squared Frobenius norm of the product over the first window, exp(2 C) + exp(-2 C), C = sum of c_m, m <= tau.
mg_G_diag(kind::Symbol, tau::Int) = (C = sum(mg_c(kind, m) for m in 1:tau); exp(2C) + exp(-2C))

# ------------------------------------------------------------------ small utilities

# Run f(); return the NumericalHealthError it throws (nothing if it returns), rethrow anything else.
function mg_health_error(f)
    try
        f()
        return nothing
    catch e
        e isa NumericalHealthError ? e : rethrow()
    end
end

# Bytes allocated by one propagate_orbit call after a warm-up call, measured inside a function.
function mg_alloc_bytes(sys, x0, W0, N, tau, reducers)
    propagate_orbit(sys, x0, W0, N, tau, reducers)
    return @allocated propagate_orbit(sys, x0, W0, N, tau, reducers)
end

# Number of elements violating a predicate: `@test mg_nbad(f, xs) == 0` keeps a failure report short (no dump of the data).
mg_nbad(f, xs) = count(x -> !f(x), xs)
mg_mean(v) = sum(v) / length(v)
mg_var(v) = (m = mg_mean(v); sum(x -> (x - m)^2, v) / (length(v) - 1))
mg_median(v) = (s = sort(collect(v)); n = length(s); isodd(n) ? s[(n + 1) ÷ 2] : (s[n ÷ 2] + s[n ÷ 2 + 1]) / 2)
# L1 S5: the difference of the ensemble means is at most 4 combined standard errors (plus an allowance).
mg_agree(a, b; allowance = 0.0) = abs(mg_mean(a) - mg_mean(b)) <= 4 * sqrt(mg_var(a) / length(a) + mg_var(b) / length(b)) + allowance

# Least-squares slope of v against u.
mg_slope(u, v) = sum((u .- mg_mean(u)) .* (v .- mg_mean(v))) / sum((u .- mg_mean(u)) .^ 2)

# E14 of L1: sigma(K) = (S + sqrt(S^2 - 4)) / 2, S^2 = 2 K^2 + 2 |K| + 5 (singular values of the standard-map tangent lie in [1/sigma, sigma]).
function mg_sigma(K::T) where {T}
    S2 = 2 * K^2 + 2 * abs(K) + 5
    return (sqrt(S2) + sqrt(S2 - 4)) / 2
end

# The series of the snapshots `sn` (reducer finals) for column i: field `name` of each snapshot, entry i.
mg_col(sn, name::Symbol, i::Int) = [getproperty(s, name)[i] for s in sn]
# The totals S_i = sum + comp of the snapshots, column i.
mg_S(sn, i::Int) = [s.sum[i] + s.comp[i] for s in sn]

# Ensemble of Y_1/n and Ybar_1/n (column 1) in type T on m starts of a class (the same Float64 starts widened exactly), through the per-orbit body with the
# bare reducer.
function mg_ens(::Type{T}, label::Symbol, m::Int, n::Int) where {T}
    c = mg_case(label, T, m)
    out = map(1:m) do j
        _, _, st = propagate_orbit(c.sys, c.x0[j], c.W0[j], n, 1, (MEGNO(),); orbit = j)
        f = reducer_final(MEGNO(), st[1])
        (Float64(f.Y[1]) / n, Float64(f.Ybar[1]) / n)
    end
    return (Y = [o[1] for o in out], Ybar = [o[2] for o in out])
end

# ------------------------------------------------------------------ Y5 and the hyperbolic fixed point (M10(b), M10(c))

# Eigen data of A = [[2, 1], [1, 1]] in BigFloat (call inside setprecision): mu and the orthonormal v_+, v_- = R90 v_+ (L1 E10).
function mg_y5_eigen()
    mu = (3 + sqrt(BigFloat(5))) / 2
    phi = (1 + sqrt(BigFloat(5))) / 2
    nv = sqrt(phi^2 + 1)
    return (mu = mu, vp = (phi / nv, 1 / nv), vm = (-1 / nv, phi / nv))
end
# Coefficients (c_1, c_2) of the two columns of W0 on (v_+, v_-): c1 = expanding, c2 = contracting (BigFloat, from the exact widened frame).
function mg_y5_coeffs(W0)
    return setprecision(BigFloat, MG_BITS) do
        e = mg_y5_eigen()
        w = ((mg_big(W0[1, 1]), mg_big(W0[2, 1])), (mg_big(W0[1, 2]), mg_big(W0[2, 2])))
        c1 = ntuple(i -> w[i][1] * e.vp[1] + w[i][2] * e.vp[2], 2)
        c2 = ntuple(i -> w[i][1] * e.vm[1] + w[i][2] * e.vm[2], 2)
        (c1 = c1, c2 = c2)
    end
end
# Supplied unit frame whose first column has expanding coefficient c1 (and contracting sqrt(1 - c1^2)); the second column is the first rotated by 90 degrees,
# so its expanding coefficient is -sqrt(1 - c1^2). Built in BigFloat, rounded into T.
function mg_y5_frame(c1::Real, ::Type{T}) where {T}
    return setprecision(BigFloat, MG_BITS) do
        e = mg_y5_eigen()
        c = BigFloat(c1)
        s = sqrt(1 - c^2)
        w1 = (c * e.vp[1] + s * e.vm[1], c * e.vp[2] + s * e.vm[2])
        w2 = (-w1[2], w1[1])
        SMatrix{2, 2, T}(T(w1[1]), T(w1[2]), T(w2[1]), T(w2[2]))
    end
end

# ------------------------------------------------------------------ rho and M_state for criterion 13 (copied in spirit from lb_ref_full: one pass)

# rho_j (sigma_max(Y_j) / norm(Y_j w), L1 criterion 7) for j = 1..n, and M_state (max(1, |K|, largest |coordinate| on the lifted orbit up to j), S1) for the
# standard map at K = Kf (Float64), start x0 and vector w0, in BigFloat(256): returns the running maximum over j <= n of rho_j and M_state at n.
function mg_rho_max(Kf::Float64, x0, w0, n::Int)
    return setprecision(BigFloat, MG_BITS) do
        K = BigFloat(Kf)
        x, p = mg_big(x0[1]), mg_big(x0[2])
        a, b, c_, d = one(K), zero(K), zero(K), one(K)
        M = max(1.0, abs(Kf))
        w1, w2 = mg_big(w0[1]), mg_big(w0[2])
        rmax = 0.0
        for _ in 1:n
            M = max(M, Float64(abs(x)), Float64(abs(p)))
            c = K * cos(x)
            a, b, c_, d = (1 + c) * a + c_, (1 + c) * b + d, c * a + c_, c * b + d
            p = p + K * sin(x)
            x = x + p
            y1, y2 = a * w1 + b * w2, c_ * w1 + d * w2
            F = a^2 + b^2 + c_^2 + d^2
            dt = a * d - b * c_
            smax = sqrt((F + sqrt(F^2 - 4 * dt^2)) / 2)
            rmax = max(rmax, Float64(smax / sqrt(y1^2 + y2^2)))
        end
        M = max(M, Float64(abs(x)), Float64(abs(p)))
        (rmax, M)
    end
end
