# Helpers for the stretching suite (prefix sn_). Definitions only: nothing here calls a StretchingNumbers, spectral_distance, spectrum_cell, spectrum_density,
# stretching_range, check_run or check_frame name at load time, so a missing symbol errors the testsets that use it and nothing else.
# Spec: workspaces/module-build/stages/01-spec/output/stretching-numbers-spec.md (D78). API proposal: workspaces/module-build/stages/02-tests/output/stretching-numbers-test-plan.md.
# The reference code of this file (histogram, D, moments, closed forms, envelopes) is plain-array code written from the spec formulas (SN1 to SN17); it never calls a function of the
# module under test to produce an expected value. lb_regular (the K = 0.5 regular-class screening, L1 S7) and lb_Psi / lb_tangents (Psi_n, L1 Amendment D35) are the single
# definitions of the lyap-benettin folder, included once behind a guard as the rli and megno suites do (this folder sorts after lyap-benettin).
using Test, Chaoslab, StaticArrays, LinearAlgebra, StableRNGs, DoubleFloats
@isdefined(lb_regular) || include(joinpath(@__DIR__, "..", "lyap-benettin", "lb_helpers.jl"))

const SN_SEED = 20261012                  # base seed, _config/numerics-defaults.toml key `seed`
const SN_BASE = SN_SEED + 37000           # this slug uses 20261012 + 37000 .. 37999; testset AC<n> uses 50 n + sub (sub < 50)
sn_seed(ac::Integer, sub::Integer = 0) = SN_BASE + 50 * ac + sub
const SN_TYPES = (Float64, Double64)      # Float32 is out (D24); BigFloat only for oracle runs
const SN_BITS = 256                       # numerics-defaults.toml key `bigfloat_bits`
# Package suite sizes (3 minute class, D38) or the full sizes of the spec: verify/stretching-numbers.jl defines SN_FULL_RUN = true before including runtests.jl.
# Every reduction is marked "D38" at its use. The testsets are top-level code (not functions), as in the other suites.
const SN_FULL = isdefined(Main, :SN_FULL_RUN) && Main.SN_FULL_RUN

sn_eps(::Type{T}) where {T} = Float64(effective_eps(T))     # Float64 eps, Double64 2^-100, BigFloat eps at the current precision
const SN_CACHE = Dict{Any, Any}()         # memoized computations shared by several criteria

# The convention tag of the spec (Algorithm, Representation), recorded in the run record.
const SN_CONVENTION = "cumulative spectra, D at every interval, columns renormalised separately, no orthogonalisation, overflow cells of width h"

# The interval of the evaluations of D in map steps must be a multiple of tau (P4): the smallest multiple of tau not below the printed 100.
sn_interval(tau::Integer) = tau * cld(100, tau)

# ------------------------------------------------------------------ BigFloat conversions

sn_big(x::Double64) = BigFloat(x.hi) + BigFloat(x.lo)       # exact; call inside setprecision(BigFloat, bits)
sn_big(x::AbstractFloat) = BigFloat(x)
sn_big(x::Integer) = BigFloat(x)
sn_err(value, ref; bits = SN_BITS) = setprecision(BigFloat, bits) do
    Float64(abs(sn_big(value) - ref))
end

# ------------------------------------------------------------------ the reference histogram and distance (SN3 to SN6, plain arrays, BigFloat)

# Cell of a BigFloat z for the range [lo, hi] and B cells (SN4): 0 below, B + 1 above, 1..B inside, the last cell closed at hi. Call inside setprecision.
function sn_cell_ref(z::BigFloat, lo::BigFloat, hi::BigFloat, B::Integer)
    z < lo && return 0
    z > hi && return B + 1
    h = (hi - lo) / B
    return min(floor(Int, (z - lo) / h) + 1, B)
end

# Counts of the cells 0..B + 1 (entry c + 1) of the BigFloat samples zs.
function sn_hist_ref(zs, lo, hi, B::Integer; bits = SN_BITS)
    return setprecision(BigFloat, bits) do
        l, u = sn_big(lo), sn_big(hi)
        N = zeros(Int, B + 2)
        for z in zs
            N[sn_cell_ref(BigFloat(z), l, u, B) + 1] += 1
        end
        N
    end
end

# Cell indices (0..B + 1) of a series, as a vector.
sn_cells_ref(zs, lo, hi, B::Integer; bits = SN_BITS) = setprecision(BigFloat, bits) do
    l, u = sn_big(lo), sn_big(hi)
    [sn_cell_ref(BigFloat(z), l, u, B) for z in zs]
end

# The exact integer sum of the squared count differences and D of SN5 in BigFloat: D = sqrt(S) / (k sqrt(h)) with h = (hi - lo) / B exact, k the (equal) total.
sn_S(N1, N2) = sum((Int(a) - Int(b))^2 for (a, b) in zip(N1, N2))
function sn_D_ref(N1, N2, lo, hi, B::Integer; bits = SN_BITS)
    k = sum(Int, N1)
    sum(Int, N2) == k || error("sn_D_ref: unequal totals")
    return setprecision(BigFloat, bits) do
        h = (sn_big(hi) - sn_big(lo)) / B
        k == 0 ? zero(BigFloat) : sqrt(BigFloat(sn_S(N1, N2))) / (k * sqrt(h))
    end
end
# D from count vectors and the width h given as a number (the pure-function reference).
sn_D_h(N1, N2, h; bits = SN_BITS) = setprecision(BigFloat, bits) do
    k = sum(Int, N1)
    k == 0 ? zero(BigFloat) : sqrt(BigFloat(sn_S(N1, N2))) / (k * sqrt(sn_big(h)))
end

# Prefix counts: for the BigFloat series zs (one per handoff) the counts after each k in ks (sorted), as a Dict k => counts.
function sn_prefix_counts(zs, lo, hi, B::Integer, ks; bits = SN_BITS)
    return setprecision(BigFloat, bits) do
        l, u = sn_big(lo), sn_big(hi)
        N = zeros(Int, B + 2)
        out = Dict{Int, Vector{Int}}()
        kset = Set(ks)
        for (k, z) in enumerate(zs)
            N[sn_cell_ref(BigFloat(z), l, u, B) + 1] += 1
            k in kset && (out[k] = copy(N))
        end
        out
    end
end

# ------------------------------------------------------------------ edge conditions (spec, Acceptance criteria, "Edge condition")

# Every sample of zs must lie at least `tol` from every interior cell edge and from lo (an overflow-side edge); at the two end edges the samples that sit inside the range must also be at
# least `endslack` from them (the default range of SN12 is padded outward by a few eps, so that rounding of z never crosses it). Samples outside [lo, hi] only need the interior test.
function sn_edge_ok(zs, lo, hi, B::Integer; tol = 1e-9, endslack = 0.0, bits = SN_BITS)
    return setprecision(BigFloat, bits) do
        l, u = sn_big(lo), sn_big(hi)
        h = (u - l) / B
        for z in zs
            zb = BigFloat(z)
            if zb < l || zb > u
                (abs(zb - l) < tol || abs(zb - u) < tol) && return false
                continue
            end
            f = (zb - l) / h
            c = round(Int, f)                              # nearest edge index
            d = abs(f - c) * h
            if c == 0 || c == B
                d < max(endslack, 0) && return false
            else
                d < tol && return false
            end
        end
        true
    end
end

# The smallest B' >= B0 (at most 60 more) for which the edge condition holds; returns (B', bumped).
function sn_pick_B(B0::Integer, lo, hi, zlists...; tol = 1e-9, endslack = 0.0)
    for B in B0:(B0 + 60)
        all(zs -> sn_edge_ok(zs, lo, hi, B; tol = tol, endslack = endslack), zlists) && return (B, B != B0)
    end
    error("sn_pick_B: no B within 60 of $B0 satisfies the edge condition")
end

# ------------------------------------------------------------------ the panel (S7), the pair of columns and the frames

sn_K(::Type{T}, label::Symbol) where {T} =
    label === :two_pi ? exact(T, pi) * 2 : label === :half ? exact(T, 0.5) : label === :ten ? exact(T, 10) :
    label === :five ? exact(T, 5) : label === :zero ? exact(T, 0) : error("unknown K $label")

sn_starts(rng, m::Integer) = [SVector(2pi * rand(rng), 2pi * rand(rng)) for _ in 1:m]    # Float64, uniform in [0, 2 pi)^2 (S7)
sn_cast(::Type{T}, z::SVector{2}) where {T} = SVector{2, T}(T(z[1]), T(z[2]))
# The seeded pair of columns of orbit id j (the first two columns of the seeded frame, independent of the number of columns, L1 S2).
sn_frame(seed, j, ::Type{T}) where {T} = initial_deviations(seed, j, Val(2), Val(2), T)
sn_astate(::WrappedStandardMap, x0) = torus_state(x0)
sn_astate(::Any, x0) = x0
sn_sys(K, wrapped::Bool) = wrapped ? WrappedStandardMap(K) : StandardMap(K)

# m starts of the panel class `label` (:half = K = 0.5 regular class; :two_pi; :ten; :five; :zero = K = 0) widened exactly into T, with their orbit ids, the stream starts (torus state in
# wrapped mode) and the seeded pairs of columns.
function sn_case(label::Symbol, ::Type{T}, m::Int; wrapped::Bool = true) where {T}
    return get!(SN_CACHE, (:case, label, T, m, wrapped)) do
        K = sn_K(T, label)
        if label === :half
            r = lb_regular()                                  # ids, starts and the frame seed of the L1 screening (one definition)
            ids, starts, fseed = r.ids[1:m], r.starts[1:m], r.seed
        else
            fseed = sn_seed(0, findfirst(==(label), (:two_pi, :ten, :five, :zero)))
            ids = collect(1:m)
            starts = sn_starts(StableRNG(fseed), m)
        end
        sys = sn_sys(K, wrapped)
        x0 = [sn_astate(sys, sn_cast(T, s)) for s in starts]
        (sys = sys, K = K, ids = ids, x0 = x0, fseed = fseed, W0 = [sn_frame(fseed, id, T) for id in ids], label = label)
    end
end

# Checkpoints and grids ------------------------------------------------------

# Grid of the Psi maximum (a test-side reading, as rli RL11: Psi_g on 1..9 and 10 per decade, the maximum over g <= n).
sn_grid(N::Int) = sort(unique(filter(<=(N), vcat(collect(1:9), [d * m for d in (10, 100, 1000, 10000, 100000) for m in 1:10]))))

# Panel sizes (class, N, starts): the spec says 20 starts per class to n = 1e4; D38: 6 starts to 1e3 and 2 starts to 1e4.
sn_panel() = SN_FULL ? ((:half, 10_000, 20), (:two_pi, 10_000, 20), (:ten, 10_000, 20)) :
             ((:half, 1000, 6), (:two_pi, 1000, 6), (:ten, 1000, 6), (:half, 10_000, 2), (:two_pi, 10_000, 2), (:ten, 10_000, 2))

# ------------------------------------------------------------------ BigFloat reference: the exact log growth s(m) = ln(norm(Y_m w0) / norm(w0)) (L1 E2, E3, E8)

function sn_pts(sys, z0, n::Integer)
    pts = Vector{typeof(z0)}(undef, n)
    z = z0
    for i in 1:n
        pts[i] = z
        z = forward(sys, z)
    end
    return pts
end

# Step-i tangent of the standard map at the realized point pts[i]: [[1 + c, 1], [c, 1]], c = K cos x (L1 E1), as (a, b, c, d).
function sn_mats_std(Kt, pts)
    return function (i)
        c = sn_big(Kt) * cos(sn_big(pts[i][1]))
        return (1 + c, one(c), c, one(c))
    end
end

function sn_ref_s(mats, w0, N::Int; bits = SN_BITS)
    return setprecision(BigFloat, bits) do
        w = (sn_big(w0[1]), sn_big(w0[2]))
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

# z_k = (s(k tau) - s((k - 1) tau)) / tau, k = 1..N / tau (SN1, with s(0) = 0).
sn_z_from_s(s::AbstractVector{BigFloat}, tau::Integer; bits = SN_BITS) = setprecision(BigFloat, bits) do
    K = length(s) ÷ tau
    [(s[k * tau] - (k == 1 ? zero(BigFloat) : s[(k - 1) * tau])) / tau for k in 1:K]
end

# z series of a stored cumulative total: the module's own compensated totals S_k = sum + comp (BigFloat), differenced.
function sn_z_from_lights(lights, j::Int, tau::Integer; bits = SN_BITS)
    return setprecision(BigFloat, bits) do
        S = [sn_big(f.sum[j]) + sn_big(f.comp[j]) for f in lights]
        [(S[k] - (k == 1 ? zero(BigFloat) : S[k - 1])) / tau for k in eachindex(S)]
    end
end

# One orbit of the standard map: both columns' step-indexed sums at the realized points, the tangents (for Psi) and the Psi memo, per column.
struct sn_Ref
    s1::Vector{BigFloat}
    s2::Vector{BigFloat}
    Ms::Vector{SMatrix{2, 2, Float64, 4}}
    w1::SVector{2, Float64}
    w2::SVector{2, Float64}
    psi::Dict{Tuple{Int, Int}, Float64}
    pmax::Dict{Int, Vector{Float64}}
end
sn_Psi(o::sn_Ref, i::Int, g::Int) = get!(o.psi, (i, g)) do
    lb_Psi(view(o.Ms, 1:g), i == 1 ? o.w1 : o.w2)
end
# pm[n] = max over the grid g <= n of Psi_g^(i) (n = 1..N); the same for every n between two grid points.
function sn_Pmax(o::sn_Ref, i::Int, N::Int)
    return get!(o.pmax, i * 10_000_000 + N) do
        grid = sn_grid(N)
        pm = zeros(N)
        cur, gi = 0.0, 1
        for n in 1:N
            while gi <= length(grid) && grid[gi] <= n
                cur = max(cur, sn_Psi(o, i, grid[gi]))
                gi += 1
            end
            pm[n] = cur
        end
        pm
    end
end

function sn_ref_orbit(sys, Kt, z0, W0, N::Int; bits = SN_BITS)
    pts = sn_pts(sys, z0, N)
    mats = sn_mats_std(Kt, pts)
    return sn_Ref(sn_ref_s(mats, W0[:, 1], N; bits = bits), sn_ref_s(mats, W0[:, 2], N; bits = bits), lb_tangents(sys, z0, N),
                  SVector{2, Float64}(Float64(W0[1, 1]), Float64(W0[2, 1])), SVector{2, Float64}(Float64(W0[1, 2]), Float64(W0[2, 2])),
                  Dict{Tuple{Int, Int}, Float64}(), Dict{Int, Vector{Float64}}())
end

# References for m starts of a class (shared by criteria 5, 6, 8, 9, 13).
function sn_refs(label::Symbol, ::Type{T}, m::Int, N::Int; wrapped::Bool = true, bits = SN_BITS) where {T}
    return get!(SN_CACHE, (:ref, label, T, m, N, wrapped, bits)) do
        c = sn_case(label, T, m; wrapped = wrapped)
        [sn_ref_orbit(c.sys, c.K, c.x0[j], c.W0[j], N; bits = bits) for j in 1:m]
    end
end

# The half-width eps_{j,k} of SN14 for the samples k = 1..K: (B_j(k tau) + B_j((k - 1) tau)) / tau + delta_edge, B_j(n) = 8 eps Psi_max^(j)(n), B_j(0) = 0.
function sn_halfwidths(o::sn_Ref, i::Int, epsT::Float64, tau::Int, K::Int, lo, hi)
    pm = sn_Pmax(o, i, K * tau)
    delta = 4 * epsT * Float64(max(abs(sn_big(lo)), abs(sn_big(hi))))
    return [(8 * epsT * (pm[k * tau] + (k == 1 ? 0.0 : pm[(k - 1) * tau]))) / tau + delta for k in 1:K]
end

# The envelope of SN14 at the prefixes ks of the series zs with half-widths eps: Dict k => (N, Nminus, Nplus, U) (counts of the cells 0..B + 1; U = samples whose interval meets two cells).
function sn_envelope(zs, eps::AbstractVector{Float64}, lo, hi, B::Integer, ks; bits = SN_BITS)
    return setprecision(BigFloat, bits) do
        l, u = sn_big(lo), sn_big(hi)
        N, Nm, Np = zeros(Int, B + 2), zeros(Int, B + 2), zeros(Int, B + 2)
        U = 0
        out = Dict{Int, NamedTuple}()
        kset = Set(ks)
        for (k, z) in enumerate(zs)
            zb = BigFloat(z)
            e = BigFloat(eps[k])
            c0 = sn_cell_ref(zb, l, u, B)
            cl, cu = sn_cell_ref(zb - e, l, u, B), sn_cell_ref(zb + e, l, u, B)
            N[c0 + 1] += 1
            if cl == cu
                Nm[cl + 1] += 1
            else
                U += 1
            end
            for c in cl:cu
                Np[c + 1] += 1
            end
            k in kset && (out[k] = (N = copy(N), Nm = copy(Nm), Np = copy(Np), U = U))
        end
        out
    end
end

# ------------------------------------------------------------------ moments (SN7), BigFloat from a stored series

# Prefix mean and standard deviation (division by the count) of the BigFloat series zs at the prefixes ks: Dict k => (mean, sd).
function sn_prefix_moments(zs, ks; bits = SN_BITS)
    return setprecision(BigFloat, bits) do
        s1, s2 = zero(BigFloat), zero(BigFloat)
        out = Dict{Int, Tuple{BigFloat, BigFloat}}()
        kset = Set(ks)
        for (k, z) in enumerate(zs)
            s1 += z
            s2 += z^2
            if k in kset
                m = s1 / k
                out[k] = (m, sqrt(max(s2 / k - m^2, zero(BigFloat))))
            end
        end
        out
    end
end

# Histogram mean and standard deviation from the counts of the cells 1..B at the cell centres (SN7); Float64.
function sn_hist_moments(counts::AbstractVector{<:Integer}, lo::Float64, h::Float64)
    B = length(counts) - 2
    n = Float64(sum(@view counts[2:(B + 1)]))
    m1 = sum(Float64(counts[c + 1]) * (lo + (c - 0.5) * h) for c in 1:B) / n
    m2 = sum(Float64(counts[c + 1]) * (lo + (c - 0.5) * h - m1)^2 for c in 1:B) / n
    return m1, sqrt(m2)
end

# ------------------------------------------------------------------ closed forms (SN8), BigFloat from the stored (exactly widened) numbers

# SN8(d), test map T1: z_{1,k} = c [x0 + p0 (tau (2k - 1) - 1) / 2], z_{2,k} = -z_{1,k}.
sn_T1_z1(c, x0, p0, tau::Int, K::Int; bits = SN_BITS) = setprecision(BigFloat, bits) do
    cb, xb, pb = sn_big(c), sn_big(x0), sn_big(p0)
    [cb * (xb + pb * (tau * (2k - 1) - 1) / 2) for k in 1:K]
end
# L1 S3 (K_s = 0) for T1: G_first = squared Frobenius norm of the product of the first tau tangents = exp(2 C) + exp(-2 C), C = c sum_{j < tau} x_j.
sn_T1_G(c, x0, p0, tau::Int) = (C = Float64(c) * sum(Float64(x0) + j * Float64(p0) for j in 0:(tau - 1)); exp(2C) + exp(-2C))

# SN8(b), Y5 A = [[2, 1], [1, 1]]: z_k = ln mu + (1 / (2 tau)) ln[(1 + r mu^(-4 k tau)) / (1 + r mu^(-4 (k - 1) tau))], r = (c2 / c1)^2, c_i the components of w on the unit eigenvectors.
function sn_y5_z(w, tau::Int, K::Int; bits = SN_BITS)
    return setprecision(BigFloat, bits) do
        mu = (3 + sqrt(BigFloat(5))) / 2
        t = mu - 2
        nrm = sqrt(1 + t^2)
        e1, e2 = (1 / nrm, t / nrm), (-t / nrm, 1 / nrm)
        wb = (sn_big(w[1]), sn_big(w[2]))
        c1, c2 = wb[1] * e1[1] + wb[2] * e1[2], wb[1] * e2[1] + wb[2] * e2[2]
        r = (c2 / c1)^2
        # the stored w is unit only to roundoff: ln |w| is not zero, so the first term is relative to |w| (convention s(0) = 0 with the stored vector's norm, as in L1)
        [log(mu) + log((1 + r * mu^(-4 * k * tau)) / (1 + r * mu^(-4 * (k - 1) * tau))) / (2 * tau) for k in 1:K]
    end
end
sn_mu(; bits = SN_BITS) = setprecision(BigFloat, bits) do
    (3 + sqrt(BigFloat(5))) / 2
end

# SN8(c), K_std = 0 shear, start (a, b): z_k = (1 / (2 tau)) ln[((a + k tau b)^2 + b^2) / ((a + (k - 1) tau b)^2 + b^2)]; for k = 1 the denominator is |w0|^2 of the stored vector
# (the convention s(0) = 0 of L1), so a vector that is unit only to roundoff needs no special case.
sn_shear_z(a, b, tau::Int, K::Int; bits = 128) = setprecision(BigFloat, bits) do
    ab, bb = sn_big(a), sn_big(b)
    [log(((ab + k * tau * bb)^2 + bb^2) / ((ab + (k - 1) * tau * bb)^2 + bb^2)) / (2 * tau) for k in 1:K]
end

# L1 E14, sigma(K) = (S + sqrt(S^2 - 4)) / 2, S^2 = 2 K^2 + 2 |K| + 5, and L_sigma = ln sigma, in BigFloat from the stored K.
sn_Lsigma(Kt; bits = SN_BITS) = setprecision(BigFloat, bits) do
    Kb = sn_big(Kt)
    S2 = 2 * Kb^2 + 2 * abs(Kb) + 5
    log((sqrt(S2) + sqrt(S2 - 4)) / 2)
end

# ------------------------------------------------------------------ the faulty doubles of the negative controls (plain arrays, never the module)

# Counts of a series with the cell rule of SN4 but a fault: :drop drops a sample outside the range, :clip puts it into the end cell, :abs stores |z|, :window keeps only the last w samples.
function sn_hist_fault(zs, lo, hi, B::Integer, fault::Symbol; window::Int = 100, bits = SN_BITS)
    return setprecision(BigFloat, bits) do
        l, u = sn_big(lo), sn_big(hi)
        N = zeros(Int, B + 2)
        ks = fault === :window ? (max(1, length(zs) - window + 1):length(zs)) : (1:length(zs))
        for k in ks
            z = BigFloat(zs[k])
            fault === :abs && (z = abs(z))
            c = sn_cell_ref(z, l, u, B)
            if fault === :drop && (c == 0 || c == B + 1)
                continue
            elseif fault === :clip
                c = c == 0 ? 1 : (c == B + 1 ? B : c)
            end
            N[c + 1] += 1
        end
        N
    end
end

# ------------------------------------------------------------------ reducers used by the tests (protocol of tangent-core with the workspace surface of the plan)

# Wraps a reducer: forwards init, workspace, step, failed, error to `inner`; after the inner step it records reducer_final(inner, state, ws) (a copy; counts and D included) at the handoffs in
# `steps`, and, when `light`, the scalar part reducer_final(inner, state) at every handoff. The stream has no emission channel (L1 Q9), so this is how a series is read. CPU only (the
# state holds Vectors): use it with propagate_orbit.
struct sn_Probe{R}
    inner::R
    steps::Vector{Int}
    light::Bool
end
sn_Probe(inner, steps = Int[]; light::Bool = false) = sn_Probe(inner, collect(Int, steps), light)

Chaoslab.reducer_init(r::sn_Probe, ::Type{T}, v::Val, k::Val) where {T} = (Chaoslab.reducer_init(r.inner, T, v, k), 1, Any[], Any[])
# reducer_workspace is a new name of the plan (Tangent workspace surface): before stage 03 it does not exist, and a method for it cannot be defined at load time.
if isdefined(Chaoslab, :reducer_workspace)
    Chaoslab.reducer_workspace(r::sn_Probe, ::Type{T}, v::Val, k::Val) where {T} = Chaoslab.reducer_workspace(r.inner, T, v, k)
end
function Chaoslab.reducer_step(r::sn_Probe, st, n::Int, x, W, ws)
    inner, pos, snaps, lights = st
    inner2, W2 = Chaoslab.reducer_step(r.inner, inner, n, x, W, ws)
    if pos <= length(r.steps) && n == r.steps[pos]
        push!(snaps, map(v -> v isa AbstractArray ? copy(v) : v, Chaoslab.reducer_final(r.inner, inner2, ws)))
        pos += 1
    end
    r.light && push!(lights, Chaoslab.reducer_final(r.inner, inner2))
    return (inner2, pos, snaps, lights), W2
end
Chaoslab.reducer_failed(r::sn_Probe, st) = Chaoslab.reducer_failed(r.inner, st[1])
Chaoslab.reducer_error(r::sn_Probe, st, orbit, backend) = Chaoslab.reducer_error(r.inner, st[1], orbit, backend)
Chaoslab.reducer_final(r::sn_Probe, st) = (snaps = st[3], lights = st[4])
Chaoslab.reducer_final(r::sn_Probe, st, ws) = (snaps = st[3], lights = st[4], inner = Chaoslab.reducer_final(r.inner, st[1], ws))

# Read-only reducer recording the orbit point at every handoff.
struct sn_Trace end
Chaoslab.reducer_init(::sn_Trace, ::Type, ::Val, ::Val) = Any[]
function Chaoslab.reducer_update(::sn_Trace, st, n::Int, x, W)
    push!(st, x)
    return st, W
end
Chaoslab.reducer_final(::sn_Trace, st) = (xs = st,)

# The workspace length of a chain (the sum over its reducers; the plan allows at most one reducer with a nonzero request) for the per-orbit body.
sn_wslen(rs::Tuple, ::Type{T}, d::Int) where {T} = sum(r -> Chaoslab.reducer_workspace(r, T, Val(d), Val(2)), rs)

# Serial stand-in for `propagate` on the CPU: the per-orbit body (propagate_orbit, with the caller's workspace) over the orbits, the chain `rs` with every reducer wrapped in a probe; returns, per
# orbit, the tuple of reducer_final(probe, state, ws) (fields snaps, lights, inner), one entry per reducer of the chain. W0s[j] is the d x 2 frame of orbit j.
function sn_run(rs::Tuple, sys, x0s, W0s, N::Integer, tau::Integer; steps = Int[], light = false, ids = eachindex(x0s))
    lights = light isa Tuple ? light : ntuple(_ -> light, length(rs))
    return map(eachindex(x0s)) do j
        probes = map((r, l) -> sn_Probe(r, steps; light = l), rs, lights)
        T = eltype(x0s[j])
        ws = zeros(UInt32, sn_wslen(probes, T, size(W0s[j], 1)))
        _, _, st = propagate_orbit(sys, x0s[j], W0s[j], N, tau, probes; orbit = ids[j], workspace = ws)
        map((p, s) -> Chaoslab.reducer_final(p, s, ws), probes, st)
    end
end

# The reducer of the criteria: the range from the derived bound of SN12 unless given, B cells, interval a multiple of tau.
sn_red(sys, B::Int, tau::Int; replace::Bool = true, range = nothing, min_per_cell::Int = 100, interval::Int = sn_interval(tau)) =
    StretchingNumbers(; range = range === nothing ? stretching_range(sys) : range, bins = B, interval = interval, min_per_cell = min_per_cell, replace = replace)

# Evaluation handoffs (map steps) of D up to N: the multiples of the interval, and the last handoff.
sn_eval_steps(N::Int, interval::Int) = sort(unique(vcat(collect(interval:interval:N), N)))

# ------------------------------------------------------------------ test maps (never in src/)

# Test map T1 of rli RL9(b): state (x, p), x' = x + p, p' = p, tangent diag(exp(c x), exp(-c x)) at the point before the step. c is built in the element type.
struct sn_T1{T} <: Chaoslab.AbstractMap
    c::T
end
Chaoslab.dimension(::sn_T1) = 2
Chaoslab.forward(::sn_T1, z) = SVector(z[1] + z[2], z[2])
Chaoslab.tangent(m::sn_T1, z::SVector{2, T}) where {T} = SMatrix{2, 2, T}(exp(m.c * z[1]), 0, 0, exp(-(m.c * z[1])))

# Constant one-step matrix diag(l1, l2): the state does not move (SN8a).
struct sn_Diag{T} <: Chaoslab.AbstractMap
    l1::T
    l2::T
end
Chaoslab.dimension(::sn_Diag) = 2
Chaoslab.forward(::sn_Diag, z) = z
Chaoslab.tangent(m::sn_Diag, z::SVector{2, T}) where {T} = SMatrix{2, 2, T}(m.l1, 0, 0, m.l2)

# Fault map on the state (n, y, flag), dimension 2: the counter advances by one per step; only flagged orbits misbehave. :zero both columns zero (zero tangent), :zero2 the second column only
# (tangent diag(1, 0)).
struct sn_Fault <: Chaoslab.AbstractMap
    kind::Symbol
end
Chaoslab.dimension(::sn_Fault) = 2
Chaoslab.forward(::sn_Fault, z) = SVector(z[1] + 1, z[2], z[3])
function Chaoslab.tangent(m::sn_Fault, z::SVector{3, T}) where {T}
    z[3] == 1 || return one(SMatrix{2, 2, T})
    m.kind === :zero && return zero(SMatrix{2, 2, T})
    return SMatrix{2, 2, T}(1, 0, 0, 0)
end

# ------------------------------------------------------------------ small utilities

function sn_health_error(f)
    try
        f()
        return nothing
    catch e
        e isa NumericalHealthError ? e : rethrow()
    end
end
function sn_error(f)
    try
        f()
        return nothing
    catch e
        return e
    end
end

# Bytes allocated by one propagate_orbit call after a warm-up call (the caller's workspace), measured inside a function.
function sn_alloc_bytes(sys, x0, W0, N, tau, reducers, ws)
    propagate_orbit(sys, x0, W0, N, tau, reducers; workspace = ws)
    return @allocated propagate_orbit(sys, x0, W0, N, tau, reducers; workspace = ws)
end

# Counts of column j of a snapshot (rows are the cells 0..B + 1) as Int.
sn_cnt(f, j) = Int.(f.counts[:, j])
sn_nbad(f, xs) = count(x -> !f(x), xs)
sn_mean(v) = sum(v) / length(v)
sn_var(v) = (m = sn_mean(v); sum(x -> (x - m)^2, v) / (length(v) - 1))
# L1 S5: the difference of the ensemble means is at most 4 combined standard errors (plus an allowance).
sn_agree(a, b; allowance = 0.0) = abs(sn_mean(a) - sn_mean(b)) <= 4 * sqrt(sn_var(a) / length(a) + sn_var(b) / length(b)) + allowance

# Multinomial sampler: Z categorical draws by inverse cdf (exact), Dirichlet(1, ..., 1) cell probabilities.
function sn_multinomial(rng, Z::Int, cum::Vector{Float64})
    N = zeros(Int, length(cum))
    for _ in 1:Z
        N[min(searchsortedfirst(cum, rand(rng)), length(cum))] += 1
    end
    return N
end
function sn_dirichlet(rng, B::Int)
    e = -log.(rand(rng, B))
    return e ./ sum(e)
end

# Bitwise equality of nested results (NaN equals NaN, -0.0 differs from 0.0).
sn_same(a, b) = isequal(a, b)
