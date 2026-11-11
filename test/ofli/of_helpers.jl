# Helpers for the ofli suite (prefix of_). Definitions only: nothing here calls a new Chaoslab name (OFLI, ofli_initial, check_frame on an OFLI) at load time;
# every call sits inside a function body or a testset, so a missing module errors the testsets that use it and nothing else.
# Spec: workspaces/module-build/stages/01-spec/output/ofli-spec.md (D118); API and conventions: workspaces/module-build/stages/02-tests/output/ofli-test-plan.md.
using Test, Chaoslab, StaticArrays, LinearAlgebra, StableRNGs, DoubleFloats

const OF_SEED = 20261012                  # _config/numerics-defaults.toml, key `seed`
const OF_BASE = OF_SEED + 43000           # this slug uses 20261012 + 43000 .. 43999; testset AC<n> uses 50 n + sub (sub < 50)
of_seed(ac::Integer, sub::Integer = 0) = OF_BASE + 50 * ac + sub
const OF_BITS = 256                       # numerics-defaults.toml, key `bigfloat_bits`
# Package sizes (D38, 3 minute ceiling for the folder) or the full sizes of the spec: verify/ofli.jl sets OF_FULL_RUN = true before including runtests.jl.
const OF_FULL = isdefined(Main, :OF_FULL_RUN) && Main.OF_FULL_RUN
of_size(small, full) = OF_FULL ? full : small

of_eps(::Type{T}) where {T} = Float64(effective_eps(T))
of_big(x::Double64) = BigFloat(x.hi) + BigFloat(x.lo)     # exact; call where the precision in force is the one wanted
of_big(x::AbstractFloat) = BigFloat(x)
of_big(x::Integer) = BigFloat(x)
const OF_CACHE = Dict{Any, Any}()                          # memoized computations shared by several criteria

of_mean(v) = sum(v) / length(v)
of_var(v) = (m = of_mean(v); sum(x -> (x - m)^2, v) / (length(v) - 1))
of_median(v) = (s = sort(collect(v)); n = length(s); isodd(n) ? s[(n + 1) ÷ 2] : (s[n ÷ 2] + s[n ÷ 2 + 1]) / 2)
of_slope(u, v) = sum((u .- of_mean(u)) .* (v .- of_mean(v))) / sum((u .- of_mean(u)) .^ 2)
of_agree(a, b; allowance = 0.0) = abs(of_mean(a) - of_mean(b)) <= 4 * sqrt(of_var(a) / length(a) + of_var(b) / length(b)) + allowance   # L1 S5

# ------------------------------------------------------------------ test fixtures (never in src/): the free particle, a rescaled field, fault maps
# Free particle H = |p|^2/2 on N degrees of freedom: separable, force zero; the exact step is the drift (fixture of the spec, no catalog row).
struct of_FreeFlow{N} <: Chaoslab.AbstractSeparableFlow end
Chaoslab.dimension(::of_FreeFlow{N}) where {N} = 2N
Chaoslab.force(::of_FreeFlow, q::SVector) = zero(q)
Chaoslab.force_jacobian(::of_FreeFlow{N}, q::SVector{N, T}) where {N, T} = zero(SMatrix{N, N, T})
Chaoslab.hamiltonian(::of_FreeFlow{N}, z::SVector) where {N} = sum(abs2, z[(N + 1):(2N)]) / 2
# The exact one-step map of step h: q += h p, tangent [[I, h I], [0, I]] (independent of the state).
struct of_FreeMap{F, T} <: Chaoslab.AbstractMap
    flow::F
    h::T
end
Chaoslab.dimension(m::of_FreeMap) = dimension(m.flow)
Chaoslab.forward(m::of_FreeMap, z::SVector) = drift(m.flow, z, m.h)
Chaoslab.tangent(m::of_FreeMap, z::SVector) = drift_jacobian(m.flow, z, m.h)
# The field of `flow` multiplied by c (invariant I4): only the reducer holds it, the map integrates the original flow.
struct of_Scaled{F, T}
    flow::F
    c::T
end
Chaoslab.dimension(s::of_Scaled) = dimension(s.flow)
Chaoslab.vector_field(s::of_Scaled, z::SVector) = s.c * vector_field(s.flow, z)
# Fault maps on the state (q, p) of the one degree of freedom flows (dimension 2): the point is fixed, so the field never vanishes at the harmonic start.
# kind :zero: zero tangent (both columns collapse at the first handoff); :col2: tangent diag(1, 0) (column 2 collapses for the frame I).
struct of_Fault <: Chaoslab.AbstractMap
    kind::Symbol
end
Chaoslab.dimension(::of_Fault) = 2
Chaoslab.forward(::of_Fault, z) = z
Chaoslab.tangent(m::of_Fault, z::SVector{2, T}) where {T} = m.kind === :zero ? zero(SMatrix{2, 2, T}) : SMatrix{2, 2, T}(1, 0, 0, 0)

# ------------------------------------------------------------------ readers of the stream (the stream has no emission channel, so a series is read through a wrapper)
# Wraps a reducer: forwards init, step, failed, error; after the inner step records reducer_final(inner) at every handoff. CPU only: use with propagate_orbit.
struct of_Probe{R}
    inner::R
end
Chaoslab.reducer_init(r::of_Probe, ::Type{T}, v::Val, k::Val) where {T} = (Chaoslab.reducer_init(r.inner, T, v, k), Any[])
function Chaoslab.reducer_step(r::of_Probe, st, n::Int, x, W)
    inner2, W2 = Chaoslab.reducer_step(r.inner, st[1], n, x, W)
    push!(st[2], Chaoslab.reducer_final(r.inner, inner2))
    return (inner2, st[2]), W2
end
Chaoslab.reducer_failed(r::of_Probe, st) = Chaoslab.reducer_failed(r.inner, st[1])
Chaoslab.reducer_error(r::of_Probe, st, orbit, backend) = Chaoslab.reducer_error(r.inner, st[1], orbit, backend)
Chaoslab.reducer_final(r::of_Probe, st) = (snaps = st[2], inner = Chaoslab.reducer_final(r.inner, st[1]))
# Read-only reducer recording the orbit point and the ARRIVING deviation matrix at every handoff (it does not replace: place it before the replacer, or alone
# for the unnormalised cumulative matrix Y W0).
struct of_Trace end
Chaoslab.reducer_init(::of_Trace, ::Type, ::Val, ::Val) = (Any[], Any[])
function Chaoslab.reducer_update(::of_Trace, st, n::Int, x, W)
    push!(st[1], x); push!(st[2], W)
    return st, W
end
Chaoslab.reducer_final(::of_Trace, st) = (xs = st[1], Ws = st[2])

# One orbit through the per-orbit body with the chain `rs` (each wrapped in a probe when `probe`): returns (x, W, per-reducer finals).
of_snaps(r, sys, x0, W0, N::Integer, tau::Integer) = begin
    p = of_Probe(r)
    _, _, st = propagate_orbit(sys, x0, W0, N, tau, (p,))
    Chaoslab.reducer_final(p, st[1]).snaps
end
# The orbit points x_0 .. x_N of the stream (iterating `forward` is bitwise the stream's orbit; pts[i + 1] = x_i).
function of_pts(sys, x0, N::Integer)
    pts = Vector{typeof(x0)}(undef, N + 1); pts[1] = x0; z = x0
    for i in 1:N
        z = forward(sys, z); pts[i + 1] = z
    end
    return pts
end
# The columns of a frame as a list of single-column frames (k = 1 runs read the per column series G_i).
of_cols(W::SMatrix{d, k, T}) where {d, k, T} = [SMatrix{d, 1, T}(W[:, i]) for i in 1:k]
of_frame(seed, j, d::Integer, k::Integer, ::Type{T}) where {T} = initial_deviations(seed, j, Val(d), Val(k), T)
# Orthonormal frame from Gram-Schmidt of random normal columns (seeded), in Float64.
function of_rotated(rng, d::Integer)
    A = randn(rng, d, d); Q = zeros(d, d)
    for j in 1:d
        v = A[:, j]
        for _ in 1:2, i in 1:(j - 1)
            v -= dot(Q[:, i], v) * Q[:, i]
        end
        Q[:, j] = v / norm(v)
    end
    return SMatrix{d, d, Float64}(Q)
end

# Bytes allocated by one propagate_orbit call after a warm-up call, measured inside a function.
function of_alloc_bytes(sys, x0, W0, N, tau, reducers)
    propagate_orbit(sys, x0, W0, N, tau, reducers)
    return @allocated propagate_orbit(sys, x0, W0, N, tau, reducers)
end
function of_health_error(f)
    try
        f()
        return nothing
    catch e
        e isa NumericalHealthError ? e : rethrow()
    end
end

# ------------------------------------------------------------------ Henon-Heiles fixtures (spec Equations 8, 9: E = 1/8, p1 rebuilt by the module's energy_start)
const OF_E = 1 / 8
const OF_W2 = SMatrix{4, 2, Float64}(1, 0, 0, 0, 0, 0, 1, 0)          # SALI frame of the printed runs: (1, 0, 0, 0) and (0, 0, 1, 0) (henon-heiles Equation 9)
of_start(q2, p2; E = OF_E) = energy_start(HenonHeilesFlow(), E, q2, p2)[1]
of_orbit_R() = of_start(0.1, 0.0)                                    # quasi-periodic orbit R (henon-heiles AC10)
of_orbit_C() = of_start(-0.015, 0.0)                                 # sticky chaotic orbit C
of_table(::Type{T}, key::Symbol) where {T} = key === :o8 ? KahanLi(T, 8) : key === :o4 ? ForestRuth(T) : key === :o2 ? Leapfrog(T) : error("unknown table $key")
const OF_H = 1 / 20
of_im(::Type{T}, key::Symbol = :o8, h = OF_H) where {T} = IntegratorMap(of_table(T, key), HenonHeilesFlow(), T(h))
# SALI at time T (steps h = 1/20, handoffs every 20 steps) of the start z0: the class label of henon-heiles criterion 10 (chaotic when below 1e-8)
function of_sali(z0, T)
    im = of_im(Float64); n = round(Int, T / OF_H)
    res = propagate(im, [z0], n; k = 2, tau = 20, seed = 1, reducers = (SALI(),), W0 = [OF_W2])
    return res.finals[1][1].sali
end
# Starts of the two classes at E = 1/8 drawn by random_energy_start (seeded), classified by SALI at Tclass: the first `m` of each, memoized. D38: Tclass = 1000
# in the package (SALI of a chaotic orbit is below 1e-8 well before), 2000 with OF_FULL (the convention of henon-heiles criterion 10).
const OF_NCLASS = of_size(8, 200)                                        # starts per class kept (package: 8, full: 200 of criterion 7)
function of_classes(m::Integer, seedoff::Integer)
    return get!(OF_CACHE, (:classes, m, seedoff)) do
        flow = HenonHeilesFlow(); a, b = section_interval(flow, OF_E); w = sqrt(2 * OF_E)
        rng = StableRNG(of_seed(0, seedoff)); Tc = of_size(1000.0, 2000.0)
        reg = SVector{4, Float64}[]; cha = SVector{4, Float64}[]; tries = 0
        while (length(reg) < m || length(cha) < m) && tries < 50 * m + 200
            tries += 1
            z, st = random_energy_start(flow, OF_E, (a, b), (-w, w), rng)
            st === :ok || continue
            s = of_sali(z, Tc)
            (s < 1.0e-8 ? (length(cha) < m && push!(cha, z)) : (length(reg) < m && push!(reg, z)))
        end
        (regular = reg, chaotic = cha)
    end
end

# ------------------------------------------------------------------ the rounding bound B_G (spec O12) and the same-orbit reference (spec O-b, O9)
# Float64 tangent M_j = tangent of the step j at the realized point x_(j-1) (pts[j]), j = 1 .. N.
of_tangents(im, pts, N::Integer) = [SMatrix{size(pts[1], 1), size(pts[1], 1), Float64}(map(Float64, tangent(im, pts[j]))) for j in 1:N]
# ln ||M_n ... M_j|| for j = 1 .. n by scaled backward products (no overflow for chaotic orbits); returned as a vector indexed by j.
function of_prodnorms(Ms, n::Integer)
    d = size(Ms[1], 1); out = zeros(n)
    P = one(SMatrix{d, d, Float64}); lP = 0.0
    for j in n:-1:1
        Q = P * Ms[j]; f = norm(Q); P = Q / f; lP += log(f)
        out[j] = lP + log(maximum(svdvals(P)))
    end
    return out
end
# Amplification sum Psi_n of one column (lyap-benettin S3, Amendment D35, d dimensional, 2-norm of the products): lw[j + 1] = ln ||w_j|| of the cumulative
# column (lw[1] = 0, unit start), lp = of_prodnorms(Ms, n).
of_Psi(lp, lw, n::Integer) = sum(j -> exp(lp[j] + lw[j] - lw[n + 1]), 1:n)           # lw[j] = ln ||w_(j-1)||
of_BG(epsT, psi, sinth) = 8 * epsT * psi + 8 * epsT / sinth                           # O12 for one maximiser

# Exact widening of a splitting method to BigFloat (the same map as the Float64 or Double64 table, not the digits string rebuilt in 256 bits).
of_widen(m::SplittingMethod{T, K, N, S}) where {T, K, N, S} =
    SplittingMethod{BigFloat, K, N, S}(map(of_big, m.c), map(of_big, m.d), m.order, m.palindromic, m.digits)

# BigFloat cumulative tangent product Y_n = M_n ... M_1 of the map `im` (type T) at the realized points pts, n = 0 .. N: Ys[n + 1]. Call outside setprecision.
function of_refY(im, pts, N::Integer; bits = OF_BITS)
    return setprecision(BigFloat, bits) do
        d = size(pts[1], 1)
        imb = IntegratorMap(of_widen(im.method), im.flow, of_big(im.tau))
        Y = one(SMatrix{d, d, BigFloat}); Ys = Vector{SMatrix{d, d, BigFloat, d * d}}(undef, N + 1); Ys[1] = Y
        for n in 1:N
            Y = tangent(imb, SVector{d, BigFloat}(map(of_big, pts[n]))) * Y; Ys[n + 1] = Y
        end
        Ys
    end
end
# Reference series of the columns of the frame W0 (type T, widened exactly) on the realized points: g (ln ||P w||), S (ln ||w||), sinth (||P w|| / ||w||), as
# k x (N + 1) BigFloat matrices; column n + 1 is the step n. The field is the exact one at the realized point (BigFloat).
function of_refseries(flow, Ys, pts, W0; bits = OF_BITS)
    return setprecision(BigFloat, bits) do
        d, k = size(W0); N = length(pts) - 1
        Wb = SMatrix{d, k, BigFloat}(map(of_big, W0))
        g = Matrix{BigFloat}(undef, k, N + 1); S = similar(g); sn = similar(g)
        for n in 0:N
            f = vector_field(flow, SVector{d, BigFloat}(map(of_big, pts[n + 1]))); ff = dot(f, f)
            Wc = Ys[n + 1] * Wb
            for i in 1:k
                w = Wc[:, i]; pw = w - (dot(f, w) / ff) * f
                a, b = norm(w), norm(pw)
                S[i, n + 1] = log(a); g[i, n + 1] = b > 0 ? log(b) : -BigFloat(Inf); sn[i, n + 1] = b / a
            end
        end
        (g = g, S = S, sinth = sn)
    end
end
# Running maxima over the handoffs n = 0 (the s = 0 term), tau, 2 tau, .. of the reference g: per column G_i(n), the maximiser step, and the basis maximum.
# Returns functions of the handoff step n (multiple of tau): (GB, stepB, colB, G1, step1).
function of_refmax(g, tau::Integer, nmax::Integer)
    k = size(g, 1); n_h = nmax ÷ tau
    GB = Vector{BigFloat}(undef, n_h + 1); stB = zeros(Int, n_h + 1); cB = zeros(Int, n_h + 1)
    G1 = similar(GB); st1 = zeros(Int, n_h + 1)
    b = -BigFloat(Inf); b1 = -BigFloat(Inf); sb = 0; s1 = 0; cb = 0
    for h in 0:n_h
        n = h * tau
        for i in 1:k
            v = g[i, n + 1]
            v > b && ((b, sb, cb) = (v, n, i))
            i == 1 && v > b1 && ((b1, s1) = (v, n))
        end
        GB[h + 1] = b; stB[h + 1] = sb; cB[h + 1] = cb; G1[h + 1] = b1; st1[h + 1] = s1
    end
    return (GB = GB, stepB = stB, colB = cB, G1 = G1, step1 = st1, tau = tau)
end

# An orbit case of henon-heiles in type T for the same-orbit criteria: the stream points, the seeded frame, the Float64 tangents, the BigFloat product Y_n
# and the reference series of the frame. `key` names the start ((:R), (:C), (:reg, i), (:cha, i)); memoized.
function of_case(key, ::Type{T}, table::Symbol, N::Integer; frame = nothing) where {T}
    return get!(OF_CACHE, (:case, key, T, table, N, frame)) do
        z0 = key === :R ? of_orbit_R() : key === :C ? of_orbit_C() :
             key[1] === :reg ? of_classes(OF_NCLASS, 1).regular[key[2]] : of_classes(OF_NCLASS, 1).chaotic[key[2]]
        x0 = SVector{4, T}(map(T, z0)); im = of_im(T, table)
        j = key === :R ? 1 : key === :C ? 2 : key[1] === :reg ? 10 + key[2] : 100 + key[2]
        W0 = frame === nothing ? of_frame(of_seed(0, 40), j, 4, 4, T) : SMatrix{4, 4, T}(map(T, frame))
        pts = of_pts(im, x0, N)
        Ys = of_refY(im, pts, N)
        ref = of_refseries(HenonHeilesFlow(), Ys, pts, W0)
        (im = im, x0 = x0, W0 = W0, pts = pts, Ys = Ys, ref = ref, Ms = of_tangents(im, pts, N), j = j, N = N)
    end
end
# Float64 B_G of the column i at the handoff step m of a case: Psi from the tangents and the reference column norms, sin theta from the reference.
function of_bound(c, T, i::Integer, m::Integer)
    m == 0 && return of_BG(of_eps(T), 0.0, Float64(c.ref.sinth[i, 1]))
    lp = get!(OF_CACHE, (:lp, objectid(c.pts), m)) do
        of_prodnorms(c.Ms, m)
    end
    lw = vcat(0.0, Float64.(view(c.ref.S, i, 2:(m + 1))))
    return of_BG(of_eps(T), of_Psi(lp, lw, m), Float64(c.ref.sinth[i, m + 1]))
end

# Compare a computed run (snapshots at every handoff of tau) with the reference at the checkpoints: returns the rows (n, form, error, bound, ratio, maxdiff).
# `g0` the computed s = 0 terms (ofli_initial). The computed full series is max(g0, snapshot maxima).
function of_rows(c, T, tau::Integer, snaps, g0, checkpoints)
    mx = of_refmax(c.ref.g, tau, c.N); rows = NamedTuple[]
    for s in snaps
        s.n in checkpoints || continue
        h = s.n ÷ tau
        # basis form: the s = 0 terms enter the computed maximum with column index and step 0
        i0 = argmax(g0); gB = max(g0[i0], s.ofli); mB = g0[i0] > s.ofli ? (i0, 0) : (s.colB, s.stepB)
        g1 = max(g0[1], s.ofli1); m1 = g0[1] > s.ofli1 ? 0 : s.step1
        errB = setprecision(BigFloat, OF_BITS) do
            Float64(abs(of_big(gB) - mx.GB[h + 1]))
        end
        err1 = setprecision(BigFloat, OF_BITS) do
            Float64(abs(of_big(g1) - mx.G1[h + 1]))
        end
        bB = max(of_bound(c, T, mB[1], mB[2]), of_bound(c, T, mx.colB[h + 1], mx.stepB[h + 1]))
        b1 = max(of_bound(c, T, 1, m1), of_bound(c, T, 1, mx.step1[h + 1]))
        push!(rows, (n = s.n, form = :B, err = errB, bound = bB, ratio = errB / bB, maxdiff = mB != (mx.colB[h + 1], mx.stepB[h + 1])))
        push!(rows, (n = s.n, form = :one, err = err1, bound = b1, ratio = err1 / b1, maxdiff = m1 != mx.step1[h + 1]))
    end
    return rows
end
# Cheap B_G of a column without the BigFloat reference (the rounding bound is negligible against the other tolerances of criteria 8, 11, 13): lw from the
# Float64 tangents, sin theta supplied (the reducer's own record at the maximiser).
function of_lw(Ms, w0, n::Integer)
    d = size(Ms[1], 1); w = SVector{d, Float64}(map(Float64, w0)); a = norm(w); w = w / a; acc = log(a)
    lw = zeros(n + 1); lw[1] = acc
    for j in 1:n
        w = Ms[j] * w; a = norm(w); acc += log(a); w = w / a; lw[j + 1] = acc
    end
    return lw
end
function of_bound_cheap(Ms, w0, m::Integer, sinth, T, cache)
    m == 0 && return of_BG(of_eps(T), 0.0, sinth)
    lp = get!(cache, m) do
        of_prodnorms(Ms, m)
    end
    return of_BG(of_eps(T), of_Psi(lp, of_lw(Ms, w0, m), m), sinth)
end
of_checkpoints(N::Int, tau::Int) = sort(unique(filter(n -> n % tau == 0 && n % 10 == 0 && n <= N, vcat((d .* (1:10) for d in (1, 10, 100, 1000, 10000))...))))

# ------------------------------------------------------------------ test-side doubles of the reducer (negative controls and the independent reference of the algebra)
# `mode` :ok is the algorithm of the spec written test-side (no compensation: for the algebra of the small systems); the other modes carry one defect each:
# :unprojected (returns the log norm), :parallel (returns the component ALONG the field instead of the orthogonal one), :lagged (field at the previous handoff point),
# :f32 (projection formed in single precision). The double has no s = 0 term of its own: `of_initial(mode, flow, x0, W0)` below.
struct of_Double
    flow::Any
    mode::Symbol
end
Chaoslab.reducer_init(::of_Double, ::Type{T}, ::Val{d}, ::Val{k}) where {T, d, k} = (zeros(T, k), fill(T(-Inf), k), T(-Inf), T(-Inf), Ref(zeros(T, d)))
function Chaoslab.reducer_step(r::of_Double, st, n::Int, x, W::SMatrix{d, k, T}) where {d, k, T}
    S, G, _, _, prev = st
    fcur = vector_field(r.flow, x)
    f = (r.mode === :lagged && !all(iszero, prev[])) ? prev[] : fcur       # the field of the previous handoff point (tau = 1: of x_(n-1))
    r.mode === :lagged && (prev[] = fcur)
    ff = dot(f, f)
    for i in 1:k
        w = W[:, i]
        a = norm(w)
        b = r.mode === :unprojected ? a : r.mode === :parallel ? abs(dot(f, w)) / sqrt(ff) :
            r.mode === :f32 ? T(norm(SVector{d, Float32}(w) - (dot(SVector{d, Float32}(f), SVector{d, Float32}(w)) / Float32(ff)) * SVector{d, Float32}(f))) : norm(w - (dot(f, w) / ff) * f)
        g = b > 0 ? S[i] + log(b) : T(-Inf)
        G[i] = max(G[i], g)
        S[i] += log(a)
    end
    return (S, G, maximum(G), G[1], prev), W ./ SVector{k, T}(ntuple(i -> norm(W[:, i]), Val(k)))'
end
Chaoslab.reducer_final(::of_Double, st) = (ofli = st[3], ofli1 = st[4], sum = copy(st[1]))
function of_initial(mode::Symbol, flow, x0, W0::SMatrix{d, k, T}) where {d, k, T}
    mode === :no_s0 && return SVector{k, T}(ntuple(_ -> T(-Inf), Val(k)))
    f = vector_field(flow, x0); ff = dot(f, f)
    return SVector{k, T}(ntuple(i -> (w = W0[:, i]; b = norm(w - (dot(f, w) / ff) * f); b > 0 ? log(b) : T(-Inf)), Val(k)))
end

# A unit-free vector perpendicular to u (Gram-Schmidt of a fixed vector), test-side.
of_perp(u) = (v = [0.3, -0.5, 0.7, 0.2]; v - dot(v, u) * u)

# Chain finals: one probe per reducer of the chain, the read series per reducer.
function of_chain(rs::Tuple, sys, x0, W0, N, tau)
    ps = map(of_Probe, rs)
    _, _, st = propagate_orbit(sys, x0, W0, N, tau, ps)
    return map((p, s) -> Chaoslab.reducer_final(p, s).snaps, ps, st)
end

