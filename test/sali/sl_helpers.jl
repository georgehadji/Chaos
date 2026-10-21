# Helpers for the sali suite (prefix sl_). Definitions only: nothing here calls an Indicators name at load time, so a missing
# Chaoslab.Indicators errors the testsets that use it and nothing else. Spec: workspaces/module-build/stages/01-spec/output/sali-spec.md
# API proposal: workspaces/module-build/stages/02-tests/output/sali-test-plan.md (top section).
using Test, Chaoslab, StaticArrays, LinearAlgebra, StableRNGs, DoubleFloats

# lb_Psi (Psi_n, L1 Amendment D35), lb_tangents, lb_pts, lb_regular (S7 screening), lb_ref_full (rho_n, M_state): one definition, in the
# lyap-benettin folder, included once.
@isdefined(lb_regular) || include(joinpath(@__DIR__, "..", "lyap-benettin", "lb_helpers.jl"))

const SL_SEED = 20261012                  # base seed, _config/numerics-defaults.toml key `seed`
const SL_BASE = SL_SEED + 31000           # this slug uses 20261012 + 31000 .. 31999; testset AC<n> uses 50 n + sub (sub < 50)
sl_seed(ac::Integer, sub::Integer = 0) = SL_BASE + 50 * ac + sub
const SL_TYPES = (Float64, Double64)      # Float32 is out (D24); BigFloat only for oracle runs
const SL_BITS = 256                       # numerics-defaults.toml key `bigfloat_bits`
# Package suite sizes (D38, 3 minute ceiling) or the full sizes of the spec: verify/sali.jl defines SL_FULL_RUN = true before including
# runtests.jl. Every reduction is marked "D38" at its use. The testsets are top-level code (not functions), as in the other suites.
const SL_FULL = isdefined(Main, :SL_FULL_RUN) && Main.SL_FULL_RUN
sl_size(small, full) = SL_FULL ? full : small

sl_eps(::Type{T}) where {T} = Float64(effective_eps(T))     # Float64 eps, Double64 2^-100, BigFloat eps at the current precision
sl_meas(::Type{T}) where {T} = 4 * sl_eps(T)                # B_meas(T) = 4 eps (A14)

# ------------------------------------------------------------------ BigFloat conversions and the exact pair

sl_big(x::Double64) = BigFloat(x.hi) + BigFloat(x.lo)       # exact; call inside setprecision(BigFloat, bits)
sl_big(x::AbstractFloat) = BigFloat(x)

# |value - ref| at `bits` bits, as Float64 (value exact in T, ref a BigFloat from sl_ref or a closed form).
sl_err(value, ref; bits = SL_BITS) = setprecision(BigFloat, bits) do
    Float64(abs(sl_big(value) - ref))
end

# Exact (d_minus, d_plus, SALI) of two vectors held in T, evaluated at `bits` bits (A2, A3); returns BigFloats.
function sl_sali_exact(a, b; bits = SL_BITS)
    return setprecision(BigFloat, bits) do
        ab, bb = map(sl_big, a), map(sl_big, b)
        ua, ub = ab / sqrt(sum(abs2, ab)), bb / sqrt(sum(abs2, bb))
        dm, dp = sqrt(sum(abs2, ua .- ub)), sqrt(sum(abs2, ua .+ ub))
        (dm, dp, min(dm, dp))
    end
end

# ------------------------------------------------------------------ the panel (S7) and the cases

sl_K(::Type{T}, label::Symbol) where {T} =
    label === :two_pi ? exact(T, pi) * 2 : label === :half ? exact(T, 0.5) : label === :ten ? exact(T, 10) :
    label === :five ? exact(T, 5) : error("unknown K $label")

sl_starts(rng, m::Integer) = [SVector(2pi * rand(rng), 2pi * rand(rng)) for _ in 1:m]    # Float64, uniform in [0, 2 pi)^2 (S7)
sl_cast(::Type{T}, z::SVector{2}) where {T} = SVector{2, T}(T(z[1]), T(z[2]))
sl_frame(seed, j, ::Type{T}) where {T} = initial_deviations(seed, j, Val(2), Val(2), T)   # the pair: first two seeded columns

const SL_CACHE = Dict{Any, Any}()          # memoized computations shared by several criteria

# m starts of the panel class `label` (:half = K = 0.5 regular class, screened as in L1 S7 by lb_regular(); :two_pi; :ten; :five) widened
# exactly into T, with their orbit ids and seeded frames.
function sl_case(label::Symbol, ::Type{T}, m::Int) where {T}
    return get!(SL_CACHE, (:case, label, T, m)) do
        K = sl_K(T, label)
        if label === :half
            ids, starts, _, fseed = lb_regular()
            ids, starts = ids[1:m], starts[1:m]
        else
            fseed = sl_seed(0, findfirst(==(label), (:two_pi, :ten, :five)))
            ids = collect(1:m)
            starts = sl_starts(StableRNG(fseed), m)
        end
        x0 = [sl_cast(T, s) for s in starts]
        (sys = StandardMap(K), K = K, ids = ids, x0 = x0, fseed = fseed, W0 = [sl_frame(fseed, id, T) for id in ids])
    end
end

# Checkpoints n <= N: multiples of 10 at 10 per decade (D38 reduction of "every handoff that is a multiple of 10"); all multiples with SL_FULL.
sl_checkpoints(N::Int) = SL_FULL ? collect(10:10:N) :
                         sort(unique(filter(n -> n % 10 == 0 && n <= N, vcat((d .* (1:10) for d in (1, 10, 100, 1000))...))))

# ------------------------------------------------------------------ BigFloat reference along an orbit (A7, A15)

# Step-i tangent of the standard map at the realized point pts[i] = x_{i-1}: [[1 + c, 1], [c, 1]], c = K cos x (E1), as (a, b, c, d).
function sl_mats_std(Kt, pts)
    return function (i)
        c = sl_big(Kt) * cos(sl_big(pts[i][1]))
        return (1 + c, one(c), c, one(c))
    end
end
sl_mats_const(a, b, c, d) = i -> (BigFloat(a), BigFloat(b), BigFloat(c), BigFloat(d))   # BigFloat at the precision in force when called

# Products Y_n W0 for n = 1..N at `bits` bits (the tangent maps come from mats(i)); W0 is the 2 x 2 frame held in T, widened exactly.
# g = abs(det W0) / (norm w1 norm w2) (A7); SALI_ref = 2 sin(atan2(abs(det W0), abs(w1 . w2)) / 2) (A3, A7) at the kept steps and
# g sqrt(2 / (1 + sqrt(1 - g^2))) (A4, no cancellation for small g) as the Float64 scale lnsali at every step.
# Returns Float64 L and lnsali for every n (scales for Lambda_n), and BigFloat L and SALI at the steps in `keep` (ascending).
function sl_ref(mats, W0, N::Int; keep = Int[], bits = SL_BITS)
    return setprecision(BigFloat, bits) do
        w1 = (sl_big(W0[1, 1]), sl_big(W0[2, 1]))
        w2 = (sl_big(W0[1, 2]), sl_big(W0[2, 2]))
        det0 = abs(w1[1] * w2[2] - w1[2] * w2[1])
        L = Vector{Float64}(undef, N)
        lns = Vector{Float64}(undef, N)
        Lb, Sb = BigFloat[], BigFloat[]
        pos = 1
        for i in 1:N
            a, b, c, d = mats(i)
            w1 = (a * w1[1] + b * w1[2], c * w1[1] + d * w1[2])
            w2 = (a * w2[1] + b * w2[2], c * w2[1] + d * w2[2])
            g = min(det0 / (sqrt(w1[1]^2 + w1[2]^2) * sqrt(w2[1]^2 + w2[2]^2)), one(det0))
            L[i] = Float64(log(g))
            gf = exp(L[i])
            lns[i] = L[i] + 0.5 * log(2 / (1 + sqrt(max(1 - gf^2, 0.0))))
            if pos <= length(keep) && keep[pos] == i
                push!(Lb, log(g))
                # SALI = 2 sin(phi/2), phi = atan2(abs(det W0), abs(w1 . w2)): stable also for an almost orthogonal pair (g near 1)
                push!(Sb, 2 * sin(atan(det0, abs(w1[1] * w2[1] + w1[2] * w2[2])) / 2))
                pos += 1
            end
        end
        (L = L, lnsali = lns, Lbig = Lb, Sbig = Sb, keep = collect(Int, keep), det0 = Float64(det0))
    end
end

# SALI_ref(0) of a frame (BigFloat).
sl_sali0(W0) = sl_sali_exact(W0[:, 1], W0[:, 2])[3]

# Value of the reference at step n (BigFloat) from sl_ref(...; keep).
function sl_at(ref, n::Int)
    i = findfirst(==(n), ref.keep)
    i === nothing && error("sl_at: step $n was not kept")
    return (L = ref.Lbig[i], S = ref.Sbig[i])
end

# ------------------------------------------------------------------ bounds of A14 and A15 (first order, constants proposed, P4)

# Lambda_n = sum_{j=1..n} SALI_ref(n) / SALI_ref(j) (sum over map steps); the ratio is clipped at e^700 (a scale, not a value).
sl_Lambda(lnsali, n::Int) = sum(j -> exp(min(lnsali[n] - lnsali[j], 700.0)), 1:n)
sl_BSALI(epsT, Lam, saliref, psi1, psi2) = 8 * epsT * (Lam + saliref * (psi1 + psi2)) + 4 * epsT      # B_SALI(n)
sl_BL(epsT, psi1, psi2) = 8 * epsT * (psi1 + psi2)                                                      # B_L(n)

# Psi_n^(1), Psi_n^(2) (lb_Psi: Float64 scale from the stored tangents) at the steps `ns`, for the frame W0.
function sl_psis(Ms, W0, ns)
    return [(lb_Psi(view(Ms, 1:n), W0[:, 1]), lb_Psi(view(Ms, 1:n), W0[:, 2])) for n in ns]
end

# Reference, tangent products and Psi for m starts of a panel class (standard map): the quantities every same-orbit criterion needs.
function sl_orbit_ref(label::Symbol, ::Type{T}, m::Int, N::Int) where {T}
    return get!(SL_CACHE, (:ref, label, T, m, N)) do
        c = sl_case(label, T, m)
        cps = sl_checkpoints(N)
        per = map(1:m) do j
            pts, _ = lb_pts(c.sys, c.x0[j], N)
            ref = sl_ref(sl_mats_std(c.K, pts), c.W0[j], N; keep = cps)
            Ms = lb_tangents(c.sys, c.x0[j], N)
            (ref = ref, psi = sl_psis(Ms, c.W0[j], cps), lam = [sl_Lambda(ref.lnsali, n) for n in cps])
        end
        (cps = cps, per = per)
    end
end

# ------------------------------------------------------------------ reducers used by the tests (protocol of tangent-core)

# Wraps a reducer: forwards init, step, failed, error to `inner`; after the inner step it records reducer_final(inner) at the handoffs
# in `steps` (or at every handoff when `all`). The stream has no emission channel (L1 Q9), so this is how a series is read. CPU only
# (the state holds a Vector): use it with propagate_orbit.
struct sl_Probe{R}
    inner::R
    steps::Vector{Int}
    all::Bool
end
sl_Probe(inner, steps::AbstractVector{<:Integer} = Int[]; all::Bool = false) = sl_Probe(inner, collect(Int, steps), all)

Chaoslab.reducer_init(r::sl_Probe, ::Type{T}, v::Val, k::Val) where {T} = (Chaoslab.reducer_init(r.inner, T, v, k), 1, Any[])
function Chaoslab.reducer_step(r::sl_Probe, st, n::Int, x, W)
    inner, pos, snaps = st
    inner2, W2 = Chaoslab.reducer_step(r.inner, inner, n, x, W)
    if r.all || (pos <= length(r.steps) && n == r.steps[pos])
        push!(snaps, Chaoslab.reducer_final(r.inner, inner2))
        pos += 1
    end
    return (inner2, pos, snaps), W2
end
Chaoslab.reducer_failed(r::sl_Probe, st) = Chaoslab.reducer_failed(r.inner, st[1])
Chaoslab.reducer_error(r::sl_Probe, st, orbit, backend) = Chaoslab.reducer_error(r.inner, st[1], orbit, backend)
Chaoslab.reducer_final(r::sl_Probe, st) = (snaps = st[3], inner = Chaoslab.reducer_final(r.inner, st[1]))

# Read-only reducer recording the orbit point at every handoff.
struct sl_Trace end
Chaoslab.reducer_init(::sl_Trace, ::Type, ::Val, ::Val) = Any[]
function Chaoslab.reducer_update(::sl_Trace, st, n::Int, x, W)
    push!(st, x)
    return st, W
end
Chaoslab.reducer_final(::sl_Trace, st) = (xs = st,)

# Serial stand-in for `propagate` on the CPU: the per-orbit body (propagate_orbit) over orbits 1:length(x0) with the reducer wrapped in a
# probe; returns reducer_final(probe) per orbit (fields snaps, inner). W0s[j] is the d x k frame of orbit j.
function sl_run(reducer, sys, x0s, W0s, N::Integer, tau::Integer; steps = Int[], all::Bool = false)
    return map(eachindex(x0s)) do j
        probe = sl_Probe(reducer, steps; all = all)
        _, _, st = propagate_orbit(sys, x0s[j], W0s[j], N, tau, (probe,); orbit = j)
        Chaoslab.reducer_final(probe, st[1])
    end
end

# Memoized SALI runs for sl_orbit_ref's starts: snapshots at the checkpoints that are handoffs of tau.
function sl_orbit_runs(label::Symbol, ::Type{T}, tau::Int, m::Int, N::Int) where {T}
    return get!(SL_CACHE, (:run, label, T, tau, m, N)) do
        c = sl_case(label, T, m)
        cps = filter(n -> n % tau == 0, sl_checkpoints(N))
        sl_run(SALI(), c.sys, c.x0, c.W0, N, tau; steps = cps)
    end
end

# One row per (start, checkpoint) of a standard-map case: values of the run, the reference and the bounds of A15 (all Float64 scales
# except the errors, which are formed against the BigFloat reference).
function sl_rows(label::Symbol, ::Type{T}, tau::Int, m::Int, N::Int) where {T}
    O = sl_orbit_ref(label, T, m, N)
    R = sl_orbit_runs(label, T, tau, m, N)
    epsT = sl_eps(T)
    rows = NamedTuple[]
    for j in 1:m, s in R[j].snaps
        i = findfirst(==(s.n), O.cps)
        ref, (p1, p2), lam = O.per[j].ref, O.per[j].psi[i], O.per[j].lam[i]
        at = sl_at(ref, s.n)
        saliref = Float64(at.S)
        push!(rows, (j = j, n = s.n, sali = s.sali, L = s.L, S1 = s.sum[1] + s.comp[1], S2 = s.sum[2] + s.comp[2],
                     errS = sl_err(s.sali, at.S), errL = sl_err(s.L, at.L), saliref = saliref, Lref = Float64(at.L),
                     BS = sl_BSALI(epsT, lam, saliref, p1, p2), BL = sl_BL(epsT, p1, p2), psi1 = p1, psi2 = p2))
    end
    return rows
end

# ------------------------------------------------------------------ test maps (never in src/)

# Fault map on the state (n, y, flag), dimension 2: the counter advances by one per step; only flagged orbits misbehave.
# kind :zero: zero tangent (both columns collapse at the first handoff); :col2: tangent diag(1, 0) (column 2 collapses when the frame is I);
# :nan: NaN in the tangent at counter value `at`.
struct sl_Fault <: Chaoslab.AbstractMap
    kind::Symbol
    at::Int
end
Chaoslab.dimension(::sl_Fault) = 2
Chaoslab.forward(::sl_Fault, z) = SVector(z[1] + 1, z[2], z[3])
function Chaoslab.tangent(m::sl_Fault, z::SVector{3, T}) where {T}
    z[3] == 1 || return one(SMatrix{2, 2, T})
    m.kind === :zero && return zero(SMatrix{2, 2, T})
    m.kind === :col2 && return SMatrix{2, 2, T}(1, 0, 0, 0)
    return SMatrix{2, 2, T}(z[1] == m.at ? T(NaN) : one(T), 0, 0, 1)
end

# Map whose forward and tangent throw: proves a rejection happens before any work.
struct sl_Boom <: Chaoslab.AbstractMap end
Chaoslab.dimension(::sl_Boom) = 2
Chaoslab.forward(::sl_Boom, z) = error("sl_Boom: forward was called")
Chaoslab.tangent(::sl_Boom, z) = error("sl_Boom: tangent was called")

# ------------------------------------------------------------------ small utilities

# Run f(); return the NumericalHealthError it throws (nothing if it returns), rethrow anything else.
function sl_health_error(f)
    try
        f()
        return nothing
    catch e
        e isa NumericalHealthError ? e : rethrow()
    end
end

# Bytes allocated by one propagate_orbit call after a warm-up call, measured inside a function.
function sl_alloc_bytes(sys, x0, W0, N, tau, reducers)
    propagate_orbit(sys, x0, W0, N, tau, reducers)
    return @allocated propagate_orbit(sys, x0, W0, N, tau, reducers)
end

sl_mean(v) = sum(v) / length(v)
sl_var(v) = (m = sl_mean(v); sum(x -> (x - m)^2, v) / (length(v) - 1))
sl_median(v) = (s = sort(collect(v)); n = length(s); isodd(n) ? s[(n + 1) ÷ 2] : (s[n ÷ 2] + s[n ÷ 2 + 1]) / 2)
# L1 S5: the difference of the ensemble means is at most 4 combined standard errors (plus an allowance).
sl_agree(a, b; allowance = 0.0) = abs(sl_mean(a) - sl_mean(b)) <= 4 * sqrt(sl_var(a) / length(a) + sl_var(b) / length(b)) + allowance

# Least-squares slope of v against u.
sl_slope(u, v) = sum((u .- sl_mean(u)) .* (v .- sl_mean(v))) / sum((u .- sl_mean(u)) .^ 2)

# Random vector in T with a full-width mantissa (Double64 gets a nonzero low part), common scale 2^e.
function sl_randvec(rng, ::Type{T}, d::Int, e::Int) where {T}
    return SVector{d, T}(ntuple(_ -> begin
        x = T(randn(rng))
        T === Double64 && (x += ldexp(T(randn(rng)), -55))
        ldexp(x, e)
    end, d))
end
