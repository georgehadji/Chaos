# Helpers for the fli-ofli suite (prefix fl_). Definitions only: nothing here calls an Indicators name at load time, so a missing
# Chaoslab.Indicators errors the testsets that use it and nothing else. Spec: workspaces/module-build/stages/01-spec/output/fli-ofli-spec.md
# API proposal: workspaces/module-build/stages/02-tests/output/fli-ofli-test-plan.md (top section).
# Only lb_psi.jl (Psi_n, L1 Amendment D35) is shared with the lyap-benettin folder, and it is included behind its own guard: this folder sorts
# before lyap-benettin, whose runtests.jl includes lb_helpers.jl unguarded, so nothing else of that folder may be loaded here first.
using Test, Chaoslab, StaticArrays, LinearAlgebra, StableRNGs, DoubleFloats
@isdefined(LB_PSI_LOADED) || include(joinpath(@__DIR__, "..", "lyap-benettin", "lb_psi.jl"))   # lb_tangents, lb_Psi

const FL_SEED = 20261012                  # base seed, _config/numerics-defaults.toml key `seed`
const FL_BASE = FL_SEED + 33000           # this slug uses 20261012 + 33000 .. 33999; testset AC<n> uses 50 n + sub (sub < 50)
fl_seed(ac::Integer, sub::Integer = 0) = FL_BASE + 50 * ac + sub
const FL_TYPES = (Float64, Double64)      # Float32 is out (D24); BigFloat only for oracle runs
const FL_BITS = 256                       # numerics-defaults.toml key `bigfloat_bits`
# Package suite sizes (D38, 3 minute ceiling) or the full sizes of the spec: verify/fli-ofli.jl defines FL_FULL_RUN = true before including
# runtests.jl. Every reduction is marked "D38" at its use. The testsets are top-level code (not functions), as in the other suites.
const FL_FULL = isdefined(Main, :FL_FULL_RUN) && Main.FL_FULL_RUN
fl_size(small, full) = FL_FULL ? full : small

fl_eps(::Type{T}) where {T} = Float64(effective_eps(T))     # Float64 eps, Double64 2^-100, BigFloat eps at the current precision

const FL_CACHE = Dict{Any, Any}()         # memoized computations shared by several criteria

# ------------------------------------------------------------------ BigFloat conversions

fl_big(x::Double64) = BigFloat(x.hi) + BigFloat(x.lo)       # exact; call inside setprecision(BigFloat, bits)
fl_big(x::AbstractFloat) = BigFloat(x)

# |value - ref| and value - ref at `bits` bits, as Float64 (value exact in T, ref a BigFloat).
fl_err(value, ref; bits = FL_BITS) = setprecision(BigFloat, bits) do
    Float64(abs(fl_big(value) - ref))
end
fl_diff(value, ref; bits = FL_BITS) = setprecision(BigFloat, bits) do
    Float64(fl_big(value) - ref)
end

# ------------------------------------------------------------------ the panel (S7) and the cases

fl_K(::Type{T}, label::Symbol) where {T} =
    label === :two_pi ? exact(T, pi) * 2 : label === :half ? exact(T, 0.5) : label === :ten ? exact(T, 10) :
    label === :five ? exact(T, 5) : error("unknown K $label")

fl_starts(rng, m::Integer) = [SVector(2pi * rand(rng), 2pi * rand(rng)) for _ in 1:m]    # Float64, uniform in [0, 2 pi)^2 (S7)
fl_cast(::Type{T}, z::SVector{2}) where {T} = SVector{2, T}(T(z[1]), T(z[2]))
fl_frame(seed, j, ::Type{T}) where {T} = initial_deviations(seed, j, Val(2), Val(2), T)   # the pair: first two seeded columns

# Regular-class starts at K = 0.5 (L1 S7), the same ids, starts and seed as lb_regular() of the lyap-benettin helpers (copied: that file cannot be
# loaded from here, see the top): candidates c = 1, 2, ... take (x, p) from the c-th draw of a StableRNG(20261012 + 21070) and the orbit id c for the
# seeded deviation vector; a candidate is regular when X(1e5) <= 5e-4 (Float64, tau = 1, Benettin). Returns (ids, starts, seed) of the first 50.
function fl_regular()
    return get!(FL_CACHE, :regular) do
        seed = FL_SEED + 21070
        sys = StandardMap(0.5)
        rng = StableRNG(seed)
        cands = SVector{2, Float64}[]
        ok = Int[]
        while length(ok) < 50
            lo = length(cands) + 1
            append!(cands, fl_starts(rng, 60))
            for c in lo:length(cands)
                W0 = initial_deviations(seed, c, Val(2), Val(1), Float64)         # orbit id = candidate number
                _, _, st = propagate_orbit(sys, cands[c], W0, 100_000, 1, (Benettin(),); orbit = c)
                reducer_final(Benettin(), st[1]).lambda[1] <= 5.0e-4 && push!(ok, c)
            end
        end
        keep = ok[1:50]
        (ids = keep, starts = cands[keep], seed = seed)
    end
end

# m starts of the panel class `label` (:half = K = 0.5 regular class; :two_pi; :ten; :five) widened exactly into T, with their orbit ids, seeded frames.
function fl_case(label::Symbol, ::Type{T}, m::Int) where {T}
    return get!(FL_CACHE, (:case, label, T, m)) do
        K = fl_K(T, label)
        if label === :half
            r = fl_regular()
            ids, starts, fseed = r.ids[1:m], r.starts[1:m], r.seed
        else
            fseed = fl_seed(0, findfirst(==(label), (:two_pi, :ten, :five)))
            ids = collect(1:m)
            starts = fl_starts(StableRNG(fseed), m)
        end
        x0 = [fl_cast(T, s) for s in starts]
        (sys = StandardMap(K), K = K, ids = ids, x0 = x0, fseed = fseed, W0 = [fl_frame(fseed, id, T) for id in ids])
    end
end

# Checkpoints n <= N: multiples of 10 at 10 per decade (D38 reduction of "every handoff that is a multiple of 10"); all multiples with FL_FULL.
fl_checkpoints(N::Int) = FL_FULL ? collect(10:10:N) :
                         sort(unique(filter(n -> n % 10 == 0 && n <= N, vcat((d .* (1:10) for d in (1, 10, 100, 1000))...))))

# Panel of the same-orbit criteria 6, 7, 10: (class, N, starts). 20 starts per class and type and n up to 1e4 for every class, as the spec says (no direct-form floor
# for FLI, and the measured package runtime allows it); D38: checkpoints at 10 per decade, FL_FULL restores every multiple of 10.
fl_panel(::Type{T}) where {T} = ((:half, 10_000, 20), (:two_pi, 10_000, 20), (:ten, 10_000, 20))

# Orbit of the standard map in T: pts[i] = x_{i-1} (the point before step i, where M_i is evaluated).
function fl_pts(sys, x0, n::Integer)
    pts = Vector{typeof(x0)}(undef, n)
    z = x0
    for i in 1:n
        pts[i] = z
        z = forward(sys, z)
    end
    return pts
end

# ------------------------------------------------------------------ BigFloat reference along an orbit (Q9, Q11, Q7)

# Step-i tangent of the standard map at the realized point pts[i]: [[1 + c, 1], [c, 1]], c = K cos x (L1 E1), as (a, b, c, d).
function fl_mats_std(Kt, pts)
    return function (i)
        c = fl_big(Kt) * cos(fl_big(pts[i][1]))
        return (1 + c, one(c), c, one(c))
    end
end
fl_mats_const(a, b, c, d) = i -> (BigFloat(a), BigFloat(b), BigFloat(c), BigFloat(d))   # BigFloat at the precision in force when called

# Squared norms q_k(i) = norm(Y_i w_k(0))^2 of the two columns of the exact widened frame W0, and G(i) = squared Frobenius norm of Y_i, for i = 1..N
# at `bits` bits (the tangent maps come from mats(i)). S_k(i) = (1/2) ln q_k(i) (L1 E2 with the unit start); no log is taken here (the running
# maximum of ln is ln of the running maximum, Q3), so the cost is one cos per step.
function fl_ref_steps(mats, W0, N::Int; bits = FL_BITS)
    return setprecision(BigFloat, bits) do
        w1 = (fl_big(W0[1, 1]), fl_big(W0[2, 1]))
        w2 = (fl_big(W0[1, 2]), fl_big(W0[2, 2]))
        a, b, c, d = one(BigFloat), zero(BigFloat), zero(BigFloat), one(BigFloat)
        q1, q2, G = Vector{BigFloat}(undef, N), Vector{BigFloat}(undef, N), Vector{BigFloat}(undef, N)
        for i in 1:N
            ma, mb, mc, md = mats(i)
            w1 = (ma * w1[1] + mb * w1[2], mc * w1[1] + md * w1[2])
            w2 = (ma * w2[1] + mb * w2[2], mc * w2[1] + md * w2[2])
            a, b, c, d = ma * a + mb * c, ma * b + mb * d, mc * a + md * c, mc * b + md * d
            q1[i] = w1[1]^2 + w1[2]^2
            q2[i] = w2[1]^2 + w2[2]^2
            G[i] = a^2 + b^2 + c^2 + d^2
        end
        (q1 = q1, q2 = q2, G = G)
    end
end

# Reference values at the handoffs of tau, at the steps `ns` (ascending multiples of tau), from fl_ref_steps output R: the running maxima over the
# handoffs tau, 2 tau, ... of the squared norms (the s = 0 term, q = 1, included: Q1, P4), as logs:
# F1 = FLI_1 (column 1), FB = FLI_B (both columns), j1 / (jB, cB) the step (0 = the s = 0 term) and column of the last strict increase (the maximisers of Q11),
# Fstar = F_*(n) = max over handoffs of ln sigma_1(Y_i) (Q7; sigma_1^2 = (G + sqrt(G^2 - 4)) / 2 for determinant one, increasing in G).
function fl_ref_at(R, tau::Int, ns::AbstractVector{<:Integer}; bits = FL_BITS)
    all(n -> n % tau == 0, ns) || error("fl_ref_at: steps must be multiples of tau = $tau")
    return setprecision(BigFloat, bits) do
        m1 = mB = one(BigFloat)
        Gm = BigFloat(2)
        j1 = jB = cB = 0
        out = NamedTuple[]
        pos = 1
        for i in tau:tau:last(ns)
            q1, q2 = R.q1[i], R.q2[i]
            q1 > m1 && (m1 = q1; j1 = i)
            q1 > mB && (mB = q1; jB = i; cB = 1)
            q2 > mB && (mB = q2; jB = i; cB = 2)
            R.G[i] > Gm && (Gm = R.G[i])
            if i == ns[pos]
                sig2 = (Gm + sqrt(max(Gm^2 - 4, zero(Gm)))) / 2
                push!(out, (n = i, F1 = log(m1) / 2, FB = log(mB) / 2, j1 = j1, jB = jB, cB = cB, Fstar = log(sig2) / 2))
                pos += 1
                pos > length(ns) && break
            end
        end
        out
    end
end

# One orbit of a same-orbit reference: squared norms, Frobenius norms, the Float64 tangents (for Psi), the frame, and memo tables.
struct fl_Orbit{M}
    q1::Vector{BigFloat}
    q2::Vector{BigFloat}
    G::Vector{BigFloat}
    Ms::Vector{SMatrix{2, 2, Float64, 4}}
    W0::M
    psi::Dict{Tuple{Int, Int}, Float64}
    at::Dict{Int, Vector{NamedTuple}}
end

# Psi_j^(i) (L1 Amendment D35, lb_Psi) of column i at step j; j = 0 is the s = 0 term, exact, so 0.
function fl_Psi(o::fl_Orbit, i::Int, j::Int)
    j == 0 && return 0.0
    return get!(o.psi, (i, j)) do
        lb_Psi(view(o.Ms, 1:j), o.W0[:, i])
    end
end

# Same-orbit reference for m starts of a panel class (standard map): BigFloat(256) (or `bits`) products at the realized points of the T orbit.
function fl_orbit_ref(label::Symbol, ::Type{T}, m::Int, N::Int; bits = FL_BITS) where {T}
    return get!(FL_CACHE, (:ref, label, T, m, N, bits)) do
        c = fl_case(label, T, m)
        per = map(1:m) do j
            pts = fl_pts(c.sys, c.x0[j], N)
            R = fl_ref_steps(fl_mats_std(c.K, pts), c.W0[j], N; bits = bits)
            fl_Orbit(R.q1, R.q2, R.G, lb_tangents(c.sys, c.x0[j], N), c.W0[j], Dict{Tuple{Int, Int}, Float64}(), Dict{Int, Vector{NamedTuple}}())
        end
        (cps = fl_checkpoints(N), per = per, case = c, bits = bits)
    end
end

# Reference values of orbit o at the checkpoints cps (multiples of tau only), memoized per tau.
fl_ref_tau!(o::fl_Orbit, tau::Int, cps; bits = FL_BITS) = get!(o.at, tau) do
    fl_ref_at(o, tau, filter(n -> n % tau == 0, cps); bits = bits)       # an fl_Orbit has the fields q1, q2, G that fl_ref_at reads
end

# ------------------------------------------------------------------ bounds (Q11, L1 S3 and Amendment D35)

# B_FLI of the two-point rule (Q11): 8 eps max(Psi_{j_a}^(i_a), Psi_{j_b}^(i_b)) for the computed maximiser (i_a, j_a) and the reference one (i_b, j_b).
fl_BFLI(epsT, o::fl_Orbit, ia, ja, ib, jb) = 8 * epsT * max(fl_Psi(o, max(ia, 1), ja), fl_Psi(o, max(ib, 1), jb))

# 8 eps (n + G_first) of L1 S3 for the constant-tangent cases, G_first = squared Frobenius norm of the product over the first window.
fl_Gfirst(M::SMatrix, tau::Int) = sum(abs2, M^tau)

# ------------------------------------------------------------------ reducers used by the tests (protocol of tangent-core)

# Wraps a reducer: forwards init, step, failed, error to `inner`; after the inner step it records reducer_final(inner) at the handoffs in `steps` (or
# at every handoff when `all`). The stream has no emission channel (L1 Q9), so this is how a series is read. CPU only (the state holds a Vector):
# use it with propagate_orbit.
struct fl_Probe{R}
    inner::R
    steps::Vector{Int}
    all::Bool
end
fl_Probe(inner, steps::AbstractVector{<:Integer} = Int[]; all::Bool = false) = fl_Probe(inner, collect(Int, steps), all)

Chaoslab.reducer_init(r::fl_Probe, ::Type{T}, v::Val, k::Val) where {T} = (Chaoslab.reducer_init(r.inner, T, v, k), 1, Any[])
function Chaoslab.reducer_step(r::fl_Probe, st, n::Int, x, W)
    inner, pos, snaps = st
    inner2, W2 = Chaoslab.reducer_step(r.inner, inner, n, x, W)
    if r.all || (pos <= length(r.steps) && n == r.steps[pos])
        push!(snaps, Chaoslab.reducer_final(r.inner, inner2))
        pos += 1
    end
    return (inner2, pos, snaps), W2
end
Chaoslab.reducer_failed(r::fl_Probe, st) = Chaoslab.reducer_failed(r.inner, st[1])
Chaoslab.reducer_error(r::fl_Probe, st, orbit, backend) = Chaoslab.reducer_error(r.inner, st[1], orbit, backend)
Chaoslab.reducer_final(r::fl_Probe, st) = (snaps = st[3], inner = Chaoslab.reducer_final(r.inner, st[1]))

# Read-only reducer recording the orbit point at every handoff.
struct fl_Trace end
Chaoslab.reducer_init(::fl_Trace, ::Type, ::Val, ::Val) = Any[]
function Chaoslab.reducer_update(::fl_Trace, st, n::Int, x, W)
    push!(st, x)
    return st, W
end
Chaoslab.reducer_final(::fl_Trace, st) = (xs = st,)

# Serial stand-in for `propagate` on the CPU: the per-orbit body (propagate_orbit) over orbits 1:length(x0s), the chain `rs` with every reducer wrapped in
# a probe; returns, per orbit, the tuple of reducer_final(probe) (fields snaps, inner), one entry per reducer of the chain. W0s[j] is the d x k frame of
# orbit j; `all` is a Bool or a tuple of Bools (one per reducer).
function fl_run(rs::Tuple, sys, x0s, W0s, N::Integer, tau::Integer; steps = Int[], all = false)
    alls = all isa Tuple ? all : ntuple(_ -> all, length(rs))
    return map(eachindex(x0s)) do j
        probes = map((r, a) -> fl_Probe(r, steps; all = a), rs, alls)
        _, _, st = propagate_orbit(sys, x0s[j], W0s[j], N, tau, probes; orbit = j)
        map(Chaoslab.reducer_final, probes, st)
    end
end

# Memoized FLI runs for fl_orbit_ref's starts: snapshots at the checkpoints that are handoffs of tau.
function fl_orbit_runs(label::Symbol, ::Type{T}, tau::Int, m::Int, N::Int) where {T}
    return get!(FL_CACHE, (:run, label, T, tau, m, N)) do
        c = fl_case(label, T, m)
        cps = filter(n -> n % tau == 0, fl_checkpoints(N))
        [r[1] for r in fl_run((FLI(),), c.sys, c.x0, c.W0, N, tau; steps = cps)]
    end
end

# One row per (start, checkpoint) of a standard-map case: the run, the BigFloat reference, the two-point bounds B_FLI of Q11 for both forms, and the
# ratios error / bound (0 when both are 0). Errors are formed in BigFloat at `bits` bits against the reference.
function fl_rows(label::Symbol, ::Type{T}, tau::Int, m::Int, N::Int; bits = FL_BITS) where {T}
    return get!(FL_CACHE, (:rows, label, T, tau, m, N, bits)) do
        O = fl_orbit_ref(label, T, m, N; bits = bits)
        R = fl_orbit_runs(label, T, tau, m, N)
        epsT = fl_eps(T)
        rows = NamedTuple[]
        for j in 1:m
            o = O.per[j]
            refs = fl_ref_tau!(o, tau, O.cps; bits = bits)
            for (s, r) in zip(R[j].snaps, refs)
                s.n == r.n || error("fl_rows: run and reference disagree on the step ($(s.n) vs $(r.n))")
                B1 = fl_BFLI(epsT, o, 1, s.step1, 1, r.j1)
                BB = fl_BFLI(epsT, o, s.colB, s.stepB, r.cB, r.jB)
                e1, eB = fl_err(s.fli1, r.F1; bits = bits), fl_err(s.fli, r.FB; bits = bits)
                ratio(e, B) = B == 0 ? (e == 0 ? 0.0 : Inf) : e / B
                push!(rows, (j = j, n = s.n, fli1 = s.fli1, fli = s.fli, S1 = s.sum[1] + s.comp[1], S2 = s.sum[2] + s.comp[2],
                             F1 = r.F1, FB = r.FB, Fstar = r.Fstar, e1 = e1, eB = eB, B1 = B1, BB = BB,
                             ratio1 = ratio(e1, B1), ratioB = ratio(eB, BB),
                             maxdiff = (s.stepB, s.colB) != (r.jB, r.cB), s = s, r = r))
            end
        end
        rows
    end
end

# ------------------------------------------------------------------ test maps (never in src/)

# Fault map on the state (n, y, flag), dimension 2: the counter advances by one per step; only flagged orbits misbehave.
# kind :zero: zero tangent (both columns collapse at the first handoff); :col2: tangent diag(1, 0) (column 2 collapses when the frame is I);
# :nan: NaN in the tangent at counter value `at`.
struct fl_Fault <: Chaoslab.AbstractMap
    kind::Symbol
    at::Int
end
Chaoslab.dimension(::fl_Fault) = 2
Chaoslab.forward(::fl_Fault, z) = SVector(z[1] + 1, z[2], z[3])
function Chaoslab.tangent(m::fl_Fault, z::SVector{3, T}) where {T}
    z[3] == 1 || return one(SMatrix{2, 2, T})
    m.kind === :zero && return zero(SMatrix{2, 2, T})
    m.kind === :col2 && return SMatrix{2, 2, T}(1, 0, 0, 0)
    return SMatrix{2, 2, T}(z[1] == m.at ? T(NaN) : one(T), 0, 0, 1)
end

# Map whose forward and tangent throw: proves a rejection happens before any work.
struct fl_Boom <: Chaoslab.AbstractMap end
Chaoslab.dimension(::fl_Boom) = 2
Chaoslab.forward(::fl_Boom, z) = error("fl_Boom: forward was called")
Chaoslab.tangent(::fl_Boom, z) = error("fl_Boom: tangent was called")

# Map with the prescribed diagonal tangent diag(a_m, 1/a_m), a_m = exp(c_m), at step m (criterion 1): the state is (counter, 0), dimension 2.
# kind :a  c_m = -0.3 for m <= 10, then +0.5;  :b  c_m = -0.2 always;  :c  c_m alternating +0.2, -0.2 (m odd: +0.2).
fl_c(kind::Symbol, m::Integer) = kind === :a ? (m <= 10 ? -0.3 : 0.5) : kind === :b ? -0.2 : (isodd(m) ? 0.2 : -0.2)
struct fl_Diag <: Chaoslab.AbstractMap
    kind::Symbol
end
Chaoslab.dimension(::fl_Diag) = 2
Chaoslab.forward(::fl_Diag, z) = SVector(z[1] + 1, z[2])
function Chaoslab.tangent(m::fl_Diag, z::SVector{2, T}) where {T}
    a = exp(T(fl_c(m.kind, round(Int, Float64(z[1])) + 1)))
    return SMatrix{2, 2, T}(a, 0, 0, 1 / a)
end

# ------------------------------------------------------------------ small utilities

# Run f(); return the NumericalHealthError it throws (nothing if it returns), rethrow anything else.
function fl_health_error(f)
    try
        f()
        return nothing
    catch e
        e isa NumericalHealthError ? e : rethrow()
    end
end

# Bytes allocated by one propagate_orbit call after a warm-up call, measured inside a function.
function fl_alloc_bytes(sys, x0, W0, N, tau, reducers)
    propagate_orbit(sys, x0, W0, N, tau, reducers)
    return @allocated propagate_orbit(sys, x0, W0, N, tau, reducers)
end

fl_mean(v) = sum(v) / length(v)
fl_var(v) = (m = fl_mean(v); sum(x -> (x - m)^2, v) / (length(v) - 1))
fl_median(v) = (s = sort(collect(v)); n = length(s); isodd(n) ? s[(n + 1) ÷ 2] : (s[n ÷ 2] + s[n ÷ 2 + 1]) / 2)
# L1 S5: the difference of the ensemble means is at most 4 combined standard errors (plus an allowance).
fl_agree(a, b; allowance = 0.0) = abs(fl_mean(a) - fl_mean(b)) <= 4 * sqrt(fl_var(a) / length(a) + fl_var(b) / length(b)) + allowance

# Least-squares slope of v against u.
fl_slope(u, v) = sum((u .- fl_mean(u)) .* (v .- fl_mean(v))) / sum((u .- fl_mean(u)) .^ 2)

# E14 of L1: sigma(K) = (S + sqrt(S^2 - 4)) / 2, S^2 = 2 K^2 + 2 |K| + 5 (singular values of the standard-map tangent lie in [1/sigma, sigma]).
function fl_sigma(K::T) where {T}
    S2 = 2 * K^2 + 2 * abs(K) + 5
    return (sqrt(S2) + sqrt(S2 - 4)) / 2
end

# Exact widened columns of a frame, squared norms of Y w for a given BigFloat 2 x 2 matrix Y (call inside setprecision).
function fl_qcols(Y, W0)
    w1 = (fl_big(W0[1, 1]), fl_big(W0[2, 1]))
    w2 = (fl_big(W0[1, 2]), fl_big(W0[2, 2]))
    y1 = (Y[1, 1] * w1[1] + Y[1, 2] * w1[2], Y[2, 1] * w1[1] + Y[2, 2] * w1[2])
    y2 = (Y[1, 1] * w2[1] + Y[1, 2] * w2[2], Y[2, 1] * w2[1] + Y[2, 2] * w2[2])
    return (y1[1]^2 + y1[2]^2, y2[1]^2 + y2[2]^2)
end

# Ensemble of FLI_B(n)/n (or the first-column form) in type T on m starts of a class (the same Float64 starts widened exactly), through the per-orbit
# body with the bare reducer.
function fl_ens(::Type{T}, label::Symbol, m::Int, n::Int; first::Bool = false) where {T}
    c = fl_case(label, T, m)
    return map(1:m) do j
        _, _, st = propagate_orbit(c.sys, c.x0[j], c.W0[j], n, 1, (FLI(),); orbit = j)
        f = reducer_final(FLI(), st[1])
        Float64(first ? f.fli1 : f.fli) / n
    end
end
