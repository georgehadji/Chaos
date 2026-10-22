# Helpers for the rli suite (prefix rl_). Definitions only: nothing here calls a PairMap, pair_start, pair_distance, RLI or check_frame name at load time,
# so a missing symbol errors the testsets that use it and nothing else. Spec: workspaces/module-build/stages/01-spec/output/rli-spec.md (D75).
# API proposal: workspaces/module-build/stages/02-tests/output/rli-test-plan.md (top section).
# lb_regular (the K = 0.5 regular-class screening, L1 S7), lb_Psi / lb_tangents (Psi_n, L1 Amendment D35) are the single definitions in the lyap-benettin
# folder, included once behind a guard as the megno suite does: this folder sorts after lyap-benettin, so under Pkg.test they are already loaded, and
# verify/rli.jl includes this suite alone. The shadow start and the separation are built here with `displaced` and `separation` of the verified
# lyap-two-trajectory module (an oracle independent of the new Systems construct), so the construct is tested against it and not against itself.
using Test, Chaoslab, StaticArrays, LinearAlgebra, StableRNGs, DoubleFloats
@isdefined(lb_regular) || include(joinpath(@__DIR__, "..", "lyap-benettin", "lb_helpers.jl"))

const RL_SEED = 20261012                  # base seed, _config/numerics-defaults.toml key `seed`
const RL_BASE = RL_SEED + 36000           # this slug uses 20261012 + 36000 .. 36999; testset AC<n> uses 50 n + sub (sub < 50)
rl_seed(ac::Integer, sub::Integer = 0) = RL_BASE + 50 * ac + sub
const RL_TYPES = (Float64, Double64)      # Float32 is out (D24); BigFloat only for oracle runs
const RL_BITS = 256                       # numerics-defaults.toml key `bigfloat_bits`
# Package suite sizes (D38, 3 minute ceiling) or the full sizes of the spec: verify/rli.jl defines RL_FULL_RUN = true before including runtests.jl.
# Every reduction is marked "D38" at its use. The testsets are top-level code (not functions), as in the other suites.
const RL_FULL = isdefined(Main, :RL_FULL_RUN) && Main.RL_FULL_RUN

rl_eps(::Type{T}) where {T} = Float64(effective_eps(T))     # Float64 eps, Double64 2^-100, BigFloat eps at the current precision
const RL_CACHE = Dict{Any, Any}()         # memoized computations shared by several criteria

# The convention tag of the spec (Algorithm, Representation), recorded in the run record.
const RL_CONVENTION = "same xi_0, shadow along xi_0, mean over the handoffs 1..K"

# ------------------------------------------------------------------ BigFloat conversions

rl_big(x::Double64) = BigFloat(x.hi) + BigFloat(x.lo)       # exact; call inside setprecision(BigFloat, bits)
rl_big(x::AbstractFloat) = BigFloat(x)
rl_err(value, ref; bits = RL_BITS) = setprecision(BigFloat, bits) do
    Float64(abs(rl_big(value) - ref))
end

# ------------------------------------------------------------------ the panel (S7), the pair and the frame

rl_K(::Type{T}, label::Symbol) where {T} =
    label === :two_pi ? exact(T, pi) * 2 : label === :half ? exact(T, 0.5) : label === :ten ? exact(T, 10) :
    label === :five ? exact(T, 5) : label === :zero ? exact(T, 0) : error("unknown K $label")

rl_starts(rng, m::Integer) = [SVector(2pi * rand(rng), 2pi * rand(rng)) for _ in 1:m]    # Float64, uniform in [0, 2 pi)^2 (S7)
rl_cast(::Type{T}, z::SVector{2}) where {T} = SVector{2, T}(T(z[1]), T(z[2]))
# xi_0: the first seeded unit deviation vector of (seed, pair id) (L1 S2); it does not depend on k.
rl_xi0(seed, j, ::Type{T}) where {T} = SVector{2, T}(initial_deviations(seed, j, Val(2), Val(1), T)[:, 1])
rl_xi2(seed, j, ::Type{T}) where {T} = SVector{2, T}(initial_deviations(seed, j, Val(2), Val(2), T)[:, 2])
# The supplied block frame of RL5: column 1 is (xi_0, 0), column 2 is (0, xi_0). Built by the caller, as the spec says.
rl_frame(xi0::SVector{d, T}) where {d, T} = SMatrix{2d, 2, T}(ntuple(l -> l <= d ? xi0[l] : (l <= 3d ? zero(T) : xi0[l - 3d]), Val(4d)))
rl_col(xi0::SVector{d, T}) where {d, T} = SMatrix{d, 1, T}(xi0)

rl_sys(K, wrapped::Bool) = wrapped ? WrappedStandardMap(K) : StandardMap(K)
rl_astate(::WrappedStandardMap, x0) = torus_state(x0)
rl_astate(::Any, x0) = x0
# Pair start from a lifted base point x0, a direction u (unit) and the nominal separation Dx: shadow = displaced(A, Dx u) of lyap-two-trajectory (reduction counts
# of A kept in wrapped mode, so the lifted separation is exactly the minimal-image one), realized separation = norm of its `separation` (L2 T2, T3).
function rl_pair(sys, x0, u, Dx)
    zA = rl_astate(sys, x0)
    zB = displaced(sys, zA, Dx * u)
    return (z = vcat(zA, zB), zA = zA, zB = zB, Dxr = safe_norm(separation(sys, zA, zB)))
end
rl_half(z::SVector{n, T}, i::Int) where {n, T} = SVector{n ÷ 2, T}(ntuple(k -> z[(i - 1) * (n ÷ 2) + k], Val(n ÷ 2)))

# m starts of the panel class `label` (:half = K = 0.5 regular class; :two_pi; :ten; :five; :zero = K = 0) widened exactly into T, with their pair ids and xi_0.
function rl_case(label::Symbol, ::Type{T}, m::Int; wrapped::Bool = true) where {T}
    return get!(RL_CACHE, (:case, label, T, m, wrapped)) do
        K = rl_K(T, label)
        if label === :half
            r = lb_regular()                                  # ids, starts and the frame seed of the L1 screening (one definition)
            ids, starts, fseed = r.ids[1:m], r.starts[1:m], r.seed
        else
            fseed = rl_seed(0, findfirst(==(label), (:two_pi, :ten, :five, :zero)))
            ids = collect(1:m)
            starts = rl_starts(StableRNG(fseed), m)
        end
        x0 = [rl_cast(T, s) for s in starts]
        (sys = rl_sys(K, wrapped), K = K, ids = ids, x0 = x0, fseed = fseed, xi0 = [rl_xi0(fseed, id, T) for id in ids], label = label)
    end
end
rl_pairs(c, Dx) = [rl_pair(c.sys, c.x0[j], c.xi0[j], Dx) for j in eachindex(c.x0)]
rl_frames(c) = [rl_frame(x) for x in c.xi0]

# Checkpoints n <= N: multiples of 10 at 10 per decade (D38 reduction of "every handoff that is a multiple of 10"); all multiples with RL_FULL.
rl_checkpoints(N::Int) = RL_FULL ? collect(10:10:N) :
                         sort(unique(filter(n -> n % 10 == 0 && n <= N, vcat((d .* (1:10) for d in (1, 10, 100, 1000))...))))

# Panel of the same-orbit criterion 3 (class, N, starts): the spec says 20 starts per class to n = 1e4; D38: 6 starts to 1e3 and 2 starts to 1e4.
rl_panel() = RL_FULL ? ((:half, 10_000, 20), (:two_pi, 10_000, 20), (:ten, 10_000, 20)) :
             ((:half, 1000, 6), (:two_pi, 1000, 6), (:ten, 1000, 6), (:half, 10_000, 2), (:two_pi, 10_000, 2), (:ten, 10_000, 2))

# ------------------------------------------------------------------ BigFloat reference: the exact sums s_i(k) with s_i(0) = 0 (L1 E2, E3, E8)

# Orbit of a single-orbit system: pts[i] = z_{i-1} (the point before step i, where M_i is evaluated).
function rl_pts(sys, z0, n::Integer)
    pts = Vector{typeof(z0)}(undef, n)
    z = z0
    for i in 1:n
        pts[i] = z
        z = forward(sys, z)
    end
    return pts
end

# Step-i tangent of the standard map at the realized point pts[i]: [[1 + c, 1], [c, 1]], c = K cos x (L1 E1), as (a, b, c, d). x is the first coordinate,
# lifted or wrapped (cos has period 2 pi; the reduced x is the number the module takes the cosine of).
function rl_mats_std(Kt, pts)
    return function (i)
        c = rl_big(Kt) * cos(rl_big(pts[i][1]))
        return (1 + c, one(c), c, one(c))
    end
end
rl_mats_const(a, b, c, d) = i -> (BigFloat(a), BigFloat(b), BigFloat(c), BigFloat(d))   # BigFloat at the precision in force when called

# s(k) = ln(norm(Y_k w0) / norm(w0)), k = 1..N, from the stored vector w0 (exact widening) and the tangent maps mats(i), at `bits` bits. Dividing by the norm of the stored
# vector makes the convention s(0) = 0 exact for a vector that is unit only to roundoff.
function rl_ref_s(mats, w0, N::Int; bits = RL_BITS)
    return setprecision(BigFloat, bits) do
        w = (rl_big(w0[1]), rl_big(w0[2]))
        l0 = log(w[1]^2 + w[2]^2) / 2
        s = Vector{BigFloat}(undef, N)
        for i in 1:N
            ma, mb, mc, md = mats(i)
            w = (ma * w[1] + mb * w[2], mc * w[1] + md * w[2])
            s[i] = log(w[1]^2 + w[2]^2) / 2 - l0
        end
        s
    end
end

# RL4 on the step-indexed sums sA, sB (BigFloat): DL_k = |sB(k tau) - sA(k tau)| / (k tau), RLI_k = (1/K) sum_{j<=k} DL_j, k = 1..N/tau. Keyword switches build the
# faulty doubles of the negative controls (criterion 1): `noabs` (no absolute value), `divn` (the mean divided by n instead of K), `mean0` (the mean over the
# handoffs 0..K-1 with DL_0 = 0 instead of 1..K), `literal` (the printed prefactor 1/t taken literally: the mean divided by tau).
function rl_ref_rli(sA::AbstractVector{BigFloat}, sB::AbstractVector{BigFloat}, tau::Integer; noabs::Bool = false, divn::Bool = false,
                    mean0::Bool = false, literal::Bool = false, bits = RL_BITS)
    return setprecision(BigFloat, bits) do
        K = length(sA) ÷ tau
        DL, RLI = Vector{BigFloat}(undef, K), Vector{BigFloat}(undef, K)
        Q = zero(BigFloat)
        for k in 1:K
            t = BigFloat(k) * tau
            d = sB[k * tau] - sA[k * tau]
            DL[k] = (noabs ? d : abs(d)) / t
            Qold = Q
            Q += DL[k]
            RLI[k] = (mean0 ? Qold : Q) / (divn ? t : k) / (literal ? tau : 1)
        end
        (DL = DL, RLI = RLI)
    end
end

rl_H(K::Integer) = sum(1 / k for k in 1:K)                  # harmonic number H_K (Float64)

# ------------------------------------------------------------------ the same-orbit reference of a pair (RL11)

# One pair: BigFloat step-indexed sums of both orbits at the realized points, the Float64 tangents (for Psi) of both orbits, the shared xi_0, the Psi memo.
struct rl_PairRef
    sA::Vector{BigFloat}
    sB::Vector{BigFloat}
    MsA::Vector{SMatrix{2, 2, Float64, 4}}
    MsB::Vector{SMatrix{2, 2, Float64, 4}}
    xi0::Any
    psi::Dict{Tuple{Int, Int}, Float64}
end
# Psi_j^(i) (L1 Amendment D35, lb_Psi) of orbit i (1 = A, 2 = B) at step j.
rl_Psi(o::rl_PairRef, i::Int, j::Int) = get!(o.psi, (i, j)) do
    lb_Psi(view(i == 1 ? o.MsA : o.MsB, 1:j), o.xi0)
end
# Psi_max^(i)(n): the maximum of Psi_g^(i) over the grid g <= n (a grid maximum, RL11).
rl_Psimax(o::rl_PairRef, i::Int, n::Int, grid) = maximum(rl_Psi(o, i, g) for g in grid if g <= n)
# B_R(n) of RL11 for the mean at the handoff K = n / tau, and the DL bound at the handoff n (RL11, first display; the last term covers the subtraction and the division).
rl_BR(epsT, PA, PB, K, n, rli_ref) = 8 * epsT * (PA + PB) * rl_H(K) / n + 4 * epsT * rli_ref
rl_BDL(epsT, PA, PB, n, dl_ref) = 8 * epsT * (PA + PB) / n + 4 * epsT * dl_ref

function rl_pair_ref(sys, Kt, zA, zB, xi0, N::Int; bits = RL_BITS)
    sA = rl_ref_s(rl_mats_std(Kt, rl_pts(sys, zA, N)), xi0, N; bits = bits)
    sB = rl_ref_s(rl_mats_std(Kt, rl_pts(sys, zB, N)), xi0, N; bits = bits)
    return rl_PairRef(sA, sB, lb_tangents(sys, zA, N), lb_tangents(sys, zB, N), xi0, Dict{Tuple{Int, Int}, Float64}())
end

# References for m starts of a class at the nominal separation Dx (shared by criteria 3 and 10(d)).
function rl_pair_refs(label::Symbol, ::Type{T}, Dx, m::Int, N::Int; wrapped::Bool = true, bits = RL_BITS) where {T}
    return get!(RL_CACHE, (:ref, label, T, Dx, m, N, wrapped, bits)) do
        c = rl_case(label, T, m; wrapped = wrapped)
        prs = rl_pairs(c, Dx)
        [rl_pair_ref(c.sys, c.K, pr.zA, pr.zB, c.xi0[j], N; bits = bits) for (j, pr) in enumerate(prs)]
    end
end

# ------------------------------------------------------------------ reducers used by the tests (protocol of tangent-core)

# Wraps a reducer: forwards init, step, failed, error to `inner`; after the inner step it records reducer_final(inner) at the handoffs in `steps` (or
# at every handoff when `all`). The stream has no emission channel (L1 Q9), so this is how a series is read. CPU only (the state holds a Vector):
# use it with propagate_orbit.
struct rl_Probe{R}
    inner::R
    steps::Vector{Int}
    all::Bool
end
rl_Probe(inner, steps = Int[]; all::Bool = false) = rl_Probe(inner, collect(Int, steps), all)

Chaoslab.reducer_init(r::rl_Probe, ::Type{T}, v::Val, k::Val) where {T} = (Chaoslab.reducer_init(r.inner, T, v, k), 1, Any[])
function Chaoslab.reducer_step(r::rl_Probe, st, n::Int, x, W)
    inner, pos, snaps = st
    inner2, W2 = Chaoslab.reducer_step(r.inner, inner, n, x, W)
    if r.all || (pos <= length(r.steps) && n == r.steps[pos])
        push!(snaps, Chaoslab.reducer_final(r.inner, inner2))
        pos += 1
    end
    return (inner2, pos, snaps), W2
end
Chaoslab.reducer_failed(r::rl_Probe, st) = Chaoslab.reducer_failed(r.inner, st[1])
Chaoslab.reducer_error(r::rl_Probe, st, orbit, backend) = Chaoslab.reducer_error(r.inner, st[1], orbit, backend)
Chaoslab.reducer_final(r::rl_Probe, st) = (snaps = st[3], inner = Chaoslab.reducer_final(r.inner, st[1]))

# Read-only reducer recording the orbit point at every handoff.
struct rl_Trace end
Chaoslab.reducer_init(::rl_Trace, ::Type, ::Val, ::Val) = Any[]
function Chaoslab.reducer_update(::rl_Trace, st, n::Int, x, W)
    push!(st, x)
    return st, W
end
Chaoslab.reducer_final(::rl_Trace, st) = (xs = st,)

# Serial stand-in for `propagate` on the CPU: the per-orbit body (propagate_orbit) over orbits 1:length(x0s), the chain `rs` with every reducer wrapped in
# a probe; returns, per orbit, the tuple of reducer_final(probe) (fields snaps, inner), one entry per reducer of the chain. W0s[j] is the d x k frame of
# orbit j; `all` is a Bool or a tuple of Bools (one per reducer). `sys` is the pair system for pair runs, the base map for single-orbit runs.
function rl_run(rs::Tuple, sys, x0s, W0s, N::Integer, tau::Integer; steps = Int[], all = false, ids = eachindex(x0s))
    alls = all isa Tuple ? all : ntuple(_ -> all, length(rs))
    return map(eachindex(x0s)) do j
        probes = map((r, a) -> rl_Probe(r, steps; all = a), rs, alls)
        _, _, st = propagate_orbit(sys, x0s[j], W0s[j], N, tau, probes; orbit = ids[j])
        map(Chaoslab.reducer_final, probes, st)
    end
end

# One row per (pair, checkpoint) of a standard-map class: the run, the BigFloat reference, the bounds of RL11 for the mean and for DL, the errors and the ratios error / bound
# (0 when both are 0). Errors are formed in BigFloat at `bits` bits against the reference. grid = the Psi_max grid (the checkpoints).
function rl_rows(label::Symbol, ::Type{T}, Dx, tau::Int, m::Int, N::Int; wrapped::Bool = true, bits = RL_BITS) where {T}
    return get!(RL_CACHE, (:rows, label, T, Dx, tau, m, N, wrapped, bits)) do
        c = rl_case(label, T, m; wrapped = wrapped)
        prs = rl_pairs(c, Dx)
        refs = rl_pair_refs(label, T, Dx, m, N; wrapped = wrapped, bits = bits)
        grid = rl_checkpoints(N)
        cps = filter(n -> n % tau == 0, grid)
        runs = rl_run((RLI(),), PairMap(c.sys), [p.z for p in prs], rl_frames(c), N, tau; steps = cps)
        epsT = rl_eps(T)
        ratio(x, b) = b == 0 ? (x == 0 ? 0.0 : Inf) : x / b
        rows = NamedTuple[]
        for j in 1:m
            o = refs[j]
            ref = rl_ref_rli(o.sA, o.sB, tau; bits = bits)
            for (f, n) in zip(runs[j][1].snaps, cps)
                f.n == n || error("rl_rows: run and reference disagree on the step ($(f.n) vs $n)")
                k = n ÷ tau
                PA, PB = rl_Psimax(o, 1, n, grid), rl_Psimax(o, 2, n, grid)
                BR = rl_BR(epsT, PA, PB, k, n, Float64(ref.RLI[k]))
                BD = rl_BDL(epsT, rl_Psi(o, 1, n), rl_Psi(o, 2, n), n, Float64(ref.DL[k]))
                eR, eD = rl_err(f.RLI, ref.RLI[k]; bits = bits), rl_err(f.DL, ref.DL[k]; bits = bits)
                push!(rows, (j = j, n = n, k = k, RLI = f.RLI, DL = f.DL, eR = eR, eD = eD, BR = BR, BD = BD, rR = ratio(eR, BR), rD = ratio(eD, BD),
                             refRLI = ref.RLI[k], refDL = ref.DL[k], PA = PA, PB = PB))
            end
        end
        rows
    end
end

# ------------------------------------------------------------------ test maps (never in src/)

# Test map T1 of RL9(b): state (x, p), x' = x + p, p' = p, tangent diag(exp(c x), exp(-c x)) at the point before the step. c is built in the element type.
struct rl_T1{T} <: Chaoslab.AbstractMap
    c::T
end
Chaoslab.dimension(::rl_T1) = 2
Chaoslab.forward(::rl_T1, z) = SVector(z[1] + z[2], z[2])
Chaoslab.tangent(m::rl_T1, z::SVector{2, T}) where {T} = SMatrix{2, 2, T}(exp(m.c * z[1]), 0, 0, exp(-(m.c * z[1])))

# Fault map on the state (n, y, flag), dimension 2: the counter advances by one per step; only flagged orbits misbehave. kind :zero: zero tangent;
# :nan: NaN in the tangent at counter value `at`. A pair of it has state (n, y, flag, n, y, flag): the flag of block A or B decides which column fails.
struct rl_Fault <: Chaoslab.AbstractMap
    kind::Symbol
    at::Int
end
Chaoslab.dimension(::rl_Fault) = 2
Chaoslab.forward(::rl_Fault, z) = SVector(z[1] + 1, z[2], z[3])
function Chaoslab.tangent(m::rl_Fault, z::SVector{3, T}) where {T}
    z[3] == 1 || return one(SMatrix{2, 2, T})
    m.kind === :zero && return zero(SMatrix{2, 2, T})
    return SMatrix{2, 2, T}(z[1] == m.at ? T(NaN) : one(T), 0, 0, 1)
end

# Map whose forward and tangent throw: proves a rejection happens before any work (a base map, so pair_start and PairMap accept it).
struct rl_Boom <: Chaoslab.AbstractMap end
Chaoslab.dimension(::rl_Boom) = 2
Chaoslab.forward(::rl_Boom, z) = error("rl_Boom: forward was called")
Chaoslab.tangent(::rl_Boom, z) = error("rl_Boom: tangent was called")

# ------------------------------------------------------------------ small utilities

# Run f(); return the NumericalHealthError it throws (nothing if it returns), rethrow anything else.
function rl_health_error(f)
    try
        f()
        return nothing
    catch e
        e isa NumericalHealthError ? e : rethrow()
    end
end

# Bytes allocated by one propagate_orbit call after a warm-up call, measured inside a function.
function rl_alloc_bytes(sys, x0, W0, N, tau, reducers)
    propagate_orbit(sys, x0, W0, N, tau, reducers)
    return @allocated propagate_orbit(sys, x0, W0, N, tau, reducers)
end

# Number of elements violating a predicate: `@test rl_nbad(f, xs) == 0` keeps a failure report short (no dump of the data).
rl_nbad(f, xs) = count(x -> !f(x), xs)
rl_mean(v) = sum(v) / length(v)
rl_var(v) = (m = rl_mean(v); sum(x -> (x - m)^2, v) / (length(v) - 1))
# L1 S5: the difference of the ensemble means is at most 4 combined standard errors (plus an allowance).
rl_agree(a, b; allowance = 0.0) = abs(rl_mean(a) - rl_mean(b)) <= 4 * sqrt(rl_var(a) / length(a) + rl_var(b) / length(b)) + allowance
# Least-squares slope of v against u.
rl_slope(u, v) = sum((u .- rl_mean(u)) .* (v .- rl_mean(v))) / sum((u .- rl_mean(u)) .^ 2)

# E14 of L1: sigma(K) = (S + sqrt(S^2 - 4)) / 2, S^2 = 2 K^2 + 2 |K| + 5 (singular values of the standard-map tangent lie in [1/sigma, sigma]).
function rl_sigma(K::T) where {T}
    S2 = 2 * K^2 + 2 * abs(K) + 5
    return (sqrt(S2) + sqrt(S2 - 4)) / 2
end

# rho_j (sigma_max(Y_j) / norm(Y_j w), L1 criterion 7) for j = 1..n, and M_state (max(1, |K|, largest |coordinate| on the lifted orbit up to j), S1) for the
# standard map at K = Kf (Float64), start x0 and vector w0, in BigFloat(256): returns the running maximum over j <= n of rho_j and M_state at n.
function rl_rho_max(Kf::Float64, x0, w0, n::Int)
    return setprecision(BigFloat, RL_BITS) do
        K = BigFloat(Kf)
        x, p = rl_big(x0[1]), rl_big(x0[2])
        a, b, c_, d = one(K), zero(K), zero(K), one(K)
        M = max(1.0, abs(Kf))
        w1, w2 = rl_big(w0[1]), rl_big(w0[2])
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

# Phi_n of RL12 (first-order resolution number, Float64): sum_{j=1..n} norm(M_n ... M_j) / norm(w(n)), w the tangent evolution of the unit vector u. Reported, never gated.
function rl_Phi(Ms::AbstractVector{<:SMatrix{2, 2}}, u)
    n = length(Ms)
    w = SVector{2, Float64}(Float64(u[1]), Float64(u[2]))
    lw = log(norm(w))
    w = w / norm(w)
    for k in 1:n
        w = SMatrix{2, 2, Float64}(Ms[k]) * w
        a = norm(w)
        lw += log(a)
        w = w / a
    end
    P = one(SMatrix{2, 2, Float64})
    lP = 0.0
    phi = 0.0
    for j in n:-1:1
        Q = P * SMatrix{2, 2, Float64}(Ms[j])
        f = norm(Q)
        P = Q / f
        lP += log(f)
        phi += exp(lP + log(lb_norm2(P)) - lw)
    end
    return phi
end

# Minimal image of a coordinate difference into (-pi, pi] and the pair distance of a pair state of the standard map computed in BigFloat from the stored numbers (the
# independent reference of pair_distance; for the Linear map the plain difference).
function rl_dist_ref(zA, zB; torus::Bool = true, bits = RL_BITS)
    return setprecision(BigFloat, bits) do
        tp = 2 * BigFloat(pi)
        mi(u) = torus ? u - tp * round(u / tp) : u
        dx, dp = mi(rl_big(zB[1]) - rl_big(zA[1])), mi(rl_big(zB[2]) - rl_big(zA[2]))
        (dx = dx, dp = dp, d = sqrt(dx^2 + dp^2))
    end
end
