# Stage 04 verification of the diagnostics module (M2 Batch A).
# Spec: workspaces/module-build/stages/01-spec/output/diagnostics-spec.md (no Amendments section).
# Run from the repo root:  julia --project=verify -t 4 verify/diagnostics.jl
# One line per criterion:  AC<n> | <quantity> | <value> | <bound> | PASS/FAIL/STAGE04-GPU/DEFERRED/INFO
# Child mode (--child) is used by AC10 to repeat the 64-orbit run with one thread; it prints a digest only.
using Chaoslab, StaticArrays, LinearAlgebra, StableRNGs, DoubleFloats, Printf, BenchmarkTools

const CHILD = "--child" in ARGS
CHILD || @eval import CairoMakie as CM

const ROOT = dirname(@__DIR__)
const OUT = joinpath(ROOT, "workspaces", "module-build", "stages", "04-verify", "output")
const SEED = 20261012                     # base seed, _config/numerics-defaults.toml key `seed`
const DET_TOL = 4.0                       # det_residual_tol
const PAIR_TOL = 16.0                     # pairing_residual_tol
const ORTH_TOL = 1.0e-6                   # absolute orthogonality tolerance for runs that do not probe it
const TYPES = (Float64, Double64)         # Float32 dropped (D24)
const PANEL = (0.5, 0.9716354, 2.0, 5.0, 10.0)
const N_LONG = parse(Int, get(ENV, "DIAG_NLONG", "1000000"))                  # AC7 steps per orbit (full size, see log for the runtime)
const N_BIG = parse(Int, get(ENV, "DIAG_NBIG", "100000"))                     # AC8 panel (256-bit det recomputation) steps per K
const N_BIGSUM = 10_000                   # AC8 panel (256-bit estimator) steps per K

include(joinpath(ROOT, "test", "diagnostics", "faultmaps.jl"))   # dg_FaultMap, dg_Renorm (test-only helpers, read only)

# ---------------------------------------------------------------- output helpers

fmt(x::AbstractFloat) = @sprintf("%.4e", Float64(x))
fmt(x) = string(x)
verdict(b) = b ? "PASS" : "FAIL"
function row(ac, qty, val, bound, v)
    println("AC", ac, " | ", qty, " | ", fmt(val), " | ", fmt(bound), " | ", v)
    flush(stdout)
end

# Run f; return the NumericalHealthError it throws, or nothing. Other exceptions propagate.
function catchhe(f)
    try
        f()
        return nothing
    catch e
        e isa NumericalHealthError || rethrow()
        return e
    end
end

# Powers of two strictly bracketing a threshold: raise side [2 theta, 4 theta), quiet side (theta/4, theta/2].
above(theta) = [exp2(j) for j in -130:20 if 2theta <= exp2(j) < 4theta]
below(theta) = [exp2(j) for j in -130:20 if theta / 4 < exp2(j) <= theta / 2]

maxabs(M) = maximum(abs, M)

# Reference orthonormalisation (modified Gram-Schmidt) of a 2 x 2 matrix, generic in T.
function mgs(W::SMatrix{2,2,T}) where {T}
    q1 = W[:, 1] / sqrt(sum(abs2, W[:, 1]))
    v = W[:, 2] - dot(q1, W[:, 2]) * q1
    return hcat(q1, v / sqrt(sum(abs2, v)))
end

# ---------------------------------------------------------------- AC1: integer maps are exact

function ac1(::Type{T}) where {T}
    mats = ((1, 0, 1, 1), (2, 1, 1, 1), (0, -1, 1, 0), (1, 1, -1, 0))   # shear, cat, rot4, rot6 (column-major)
    maxd = zero(T); maxs = zero(T); maxfin = zero(T); errs = 0
    for ent in mats
        sys = LinearSymplecticMap(SMatrix{2,2,T}(ent...))
        mon = ResidualMonitor(sys, DET_TOL, ORTH_TOL)
        st = reducer_init(mon, T, Val(2), Val(2))
        W = SMatrix{2,2,T}(1, 0, 0, 1)
        z = SVector{2,T}(1, 2)
        try
            for j in 1:1000
                M = tangent(sys, z)
                maxd = max(maxd, eps_det(M)); maxs = max(maxs, eps_symp(M))
                check_residual("det", eps_det(M), det_threshold(DET_TOL, T, M); step = j, orbit = 1)
                st, W = reducer_update(mon, st, j, z, W)
                z = wrap(forward(sys, z))
            end
        catch e
            e isa NumericalHealthError || rethrow()
            errs += 1
        end
        f = reducer_final(mon, st)
        maxfin = max(maxfin, f.det_max, f.symp_max, f.orth_max)
    end
    row(1, "$T integer maps x1000 steps: max(eps_det, eps_symp, monitor maxima), health errors=$errs",
        Float64(max(maxd, maxs, maxfin)), 0.0, verdict(maxd == 0 && maxs == 0 && maxfin == 0 && errs == 0))
end

# ---------------------------------------------------------------- AC2: rotation by 0.7311

function ac2(::Type{T}, n) where {T}
    sys = LinearSymplecticMap(rotation_block(T, T(0.7311)))
    z = SVector{2,T}(1, 0)
    Q = SMatrix{2,2,T}(1, 0, 0, 1)
    maxdet = zero(T); maxorth = zero(T); viol = 0
    for _ in 1:n
        M = tangent(sys, z)
        d = eps_det(M)
        maxdet = max(maxdet, d)
        d <= det_threshold(DET_TOL, T, M) || (viol += 1)
        z = forward(sys, z)
        Q = mgs(M * Q)
        maxorth = max(maxorth, eps_orth(Q))
    end
    e = eps(T)
    row(2, "$T rotation 0.7311, $n steps: max eps_det/eps (threshold violations=$viol)", Float64(maxdet / e), 8.0,
        verdict(viol == 0 && maxdet <= 8e))
    row(2, "$T rotation 0.7311, $n steps: max eps_orth/eps (reference MGS)", Float64(maxorth / e), 8.0,
        verdict(maxorth <= 8e))
end

# ---------------------------------------------------------------- reference QR estimator (AC3, AC7, AC8c)

# Window = one step (tau = 1). Matrix QR of DT*Q, log|R_ii| in CompensatedSum, lambda_i(n) = sum/n (estimator order).
# G from R (same factorisation): sigma_1^2 = (S + sqrt(S^2 - 4 d^2))/2 with S = |R|_F^2, d = |R11 R22| = prod sigma;
# G = sigma_1/sigma_2 = sigma_1^2/d. Returns the running-maximum-normalised sum residual at every n (log-binned maxima kept).
const BINS_PER_DECADE = 20
function qr_estimator(sys, advance, x0::SVector{2,T}, n; Q0 = Matrix{T}(I, 2, 2), reference = nothing) where {T}
    z = x0
    Q = Q0
    a1 = CompensatedSum{T}(); a2 = CompensatedSum{T}()
    Gmax = zero(T); Mmax = zero(T)
    ratio_max = 0.0; sum_max = zero(T); detratio = 0.0; viol4 = 0; viol16 = 0
    nb = BINS_PER_DECADE * ceil(Int, log10(n)) + 1
    bins = zeros(nb)
    lam = SVector{2,T}(0, 0)
    for k in 1:n
        M = tangent(sys, z)
        z = advance(z)
        F = qr(Matrix(M) * Q)
        Q = Matrix(F.Q)
        r1 = abs(F.R[1, 1]); r2 = abs(F.R[2, 2])
        a1 = add(a1, log(r1)); a2 = add(a2, log(r2))
        lam = SVector{2,T}(total(a1) / k, total(a2) / k)
        S = r1^2 + F.R[1, 2]^2 + r2^2; d = r1 * r2
        Gmax = max(Gmax, (S + sqrt(S^2 - 4 * d^2)) / 2 / d)
        Mmax = max(Mmax, maxabs(M))
        es = eps_sum(lam)
        sum_max = max(sum_max, es)
        r = Float64(es / (eps(T) * (Gmax + Mmax)))
        ratio_max = max(ratio_max, r)
        r > 4 && (viol4 += 1)
        r > PAIR_TOL && (viol16 += 1)
        detratio = max(detratio, Float64(eps_det(M) / (eps(T) * maxabs(M))))
        b = floor(Int, BINS_PER_DECADE * log10(k)) + 1
        bins[b] = max(bins[b], r)
    end
    return (; ratio_max, sum_max, Gmax, Mmax, detratio, viol4, viol16, bins, lam)
end

# ---------------------------------------------------------------- AC3: cat map, reference QR, n = 1e5

function ac3(n)
    sys = LinearSymplecticMap(SMatrix{2,2,Float64}(2, 1, 1, 1))
    Q0 = [cos(0.3) -sin(0.3); sin(0.3) cos(0.3)]
    r = qr_estimator(sys, z -> wrap(forward(sys, z)), SVector(1.0, 2.0), n; Q0)
    Gexact = ((3 + sqrt(5)) / 2)^2
    row(3, "cat map tau=1 n=$n: G_max (vs mu^2 = $(round(Gexact; digits = 4)))", r.Gmax, Gexact, verdict(isapprox(r.Gmax, Gexact; rtol = 1e-12)))
    row(3, "cat map: max_n eps_sum / (eps (G_max/tau + M_max)), violations of theta_sum=$(r.viol16)", r.ratio_max, PAIR_TOL,
        verdict(r.viol16 == 0))
    # N = 1: eps_pair is the same quantity as eps_sum (same expression, gap 0 by construction)
    lam = r.lam
    row(3, "cat map: |eps_pair - eps_sum| at n=$n; |lambda_1 - ln mu|", abs(eps_pair(lam) - eps_sum(lam)),
        0.0, verdict(eps_pair(lam) == eps_sum(lam) && abs(lam[1] - log((3 + sqrt(5)) / 2)) <= 1e-3))
end

# ---------------------------------------------------------------- AC4: pairing index mapping

function ac4(::Type{T}) where {T}
    mu = (3 + sqrt(T(5))) / 2
    L = log(mu); d = T(2)^-20
    lam1 = SVector{2,T}(L + d, -L)
    a = T(96) / 100; b = T(31) / 100; d1 = T(2)^-20; d2 = T(2)^-18
    lam2 = SVector{4,T}(a, b, -b + d2, -a + d1)
    e1 = max(abs(eps_pair(lam1) - d), abs(eps_sum(lam1) - d))
    tol1 = eps(T) * (abs(L) + d)
    e2 = max(abs(eps_pair(lam2) - d2), abs(eps_sum(lam2) - (d1 + d2)))
    tol2 = eps(T) * 2 * maximum(abs, lam2)
    adjacent = maximum(i -> abs(lam2[i] + lam2[i+1]), 1:2:3)   # what (i, i+1) pairing would give
    row(4, "$T N=1: max(|eps_pair-d|,|eps_sum-d|)/tol", Float64(e1 / tol1), 1.0, verdict(e1 <= tol1))
    row(4, "$T N=2: max(|eps_pair-d2|,|eps_sum-(d1+d2)|)/tol, pair=$(Float64(eps_pair(lam2))) (adjacent would be $(Float64(adjacent)))",
        Float64(e2 / tol2), 1.0, verdict(e2 <= tol2 && eps_pair(lam2) < 0.5))
end

# ---------------------------------------------------------------- AC5: injected faults

function ac5a(::Type{T}) where {T}
    A0 = SMatrix{2,2,T}(2, 1, 1, 1)
    theta = Float64(det_threshold(DET_TOL, T, A0))
    x0 = SVector{2,T}(0.5, 0.25)
    W0 = SMatrix{2,1,T}(1, 0)
    run(delta, N) = begin
        sys = dg_FaultMap(SMatrix{2,2,T}(2 + T(delta), 1, 1, 1))
        mon = ResidualMonitor(sys, DET_TOL, ORTH_TOL)
        catchhe(() -> propagate_orbit(sys, x0, W0, N, 1, (dg_Renorm(), mon); orbit = 7))
    end
    ups = [(d, run(d, 1)) for d in above(theta)]
    dns = [(d, run(d, 10_000)) for d in below(theta)]
    okU = count(((d, e),) -> e isa NumericalHealthError && e.what == "det" && e.value == d && e.orbit == 7 && e.step == 1, ups)
    okD = count(((d, e),) -> e === nothing, dns)
    row(5, "det fault $T theta=$(fmt(theta)): raise-side ok (bitwise value) / quiet-side ok", "$okU/$(length(ups)), $okD/$(length(dns))",
        "all", verdict(okU == length(ups) > 0 && okD == length(dns) > 0))
end

function ac5b(::Type{T}, theta_o) where {T}
    sys = dg_FaultMap(SMatrix{2,2,T}(1, 0, 0, 1))
    mon = ResidualMonitor(sys, DET_TOL, theta_o)
    x = SVector{2,T}(0.5, 0.25)
    run(delta) = begin
        Q = SMatrix{2,2,T}(1 + T(delta), 0, 0, 1)
        st = reducer_init(mon, T, Val(2), Val(2))
        (delta, eps_orth(Q), catchhe(() -> reducer_update(mon, st, 5, x, Q)))
    end
    ups = [run(r / 2) for r in above(theta_o)]
    dns = [run(r / 2) for r in below(theta_o)]
    err(c) = abs(c[2] - (2 * c[1] + c[1]^2))
    okU = count(c -> c[3] isa NumericalHealthError && c[3].what == "orthogonality" && c[3].step == 5 &&
                     abs(c[3].value - (2 * c[1] + c[1]^2)) <= eps(T) && err(c) <= eps(T), ups)
    okD = count(c -> c[3] === nothing && err(c) <= eps(T), dns)
    worst = maximum(c -> Float64(err(c)) / Float64(eps(T)), vcat(ups, dns))
    row(5, "orthogonality fault $T tol=$(fmt(theta_o)): raise ok / quiet ok; max |eps_orth-(2d+d^2)|/eps=$(round(worst; digits = 3))",
        "$okU/$(length(ups)), $okD/$(length(dns))", "all", verdict(okU == length(ups) > 0 && okD == length(dns) > 0))
end

function ac5c(::Type{T}) where {T}
    mu = (3 + sqrt(T(5))) / 2
    L = log(mu); G = mu^2; Mmax = T(2)
    theta = sum_threshold(PAIR_TOL, T, G, 1, Mmax)
    mon = ResidualMonitor(StandardMap(T(2)), DET_TOL, ORTH_TOL)
    tol = eps(T) * (abs(L) + 1)
    run(delta) = begin
        lam = SVector{2,T}(L + T(delta), -L)
        st = reducer_init(mon, T, Val(2), Val(2))
        em = catchhe(() -> monitor_exponents(st, lam, G, 1; pairing_tol = PAIR_TOL, tau = 1, M_max = Mmax, step = 3, orbit = 9))
        ep = catchhe(() -> check_residual("pairing", eps_pair(lam), theta; step = 3, orbit = 9))
        (; delta, es = eps_sum(lam), ep = eps_pair(lam), em, epe = ep)
    end
    thf = Float64(theta)
    ups = [run(d) for d in above(thf)]; dns = [run(d) for d in below(thf)]
    val(c) = abs(c.es - c.delta) <= tol && abs(c.ep - c.delta) <= tol
    okU = count(c -> val(c) && c.em isa NumericalHealthError && c.em.what == "sum" && c.em.step == 3 && c.em.orbit == 9 &&
                     c.epe isa NumericalHealthError && c.epe.what == "pairing", ups)
    okD = count(c -> val(c) && c.em === nothing && c.epe === nothing, dns)
    row(5, "sum+pairing fault $T theta_sum=$(fmt(thf)): raise ok / quiet ok (values within eps(|L|+1))",
        "$okU/$(length(ups)), $okD/$(length(dns))", "all", verdict(okU == length(ups) > 0 && okD == length(dns) > 0))
end

function ac5d(::Type{T}) where {T}
    B1 = SMatrix{2,2,T}(2, 1, 1, 1)
    theta = symp_threshold(DET_TOL, T, block_diagonal(B1, B1))
    run(delta) = begin
        A = block_diagonal(B1, SMatrix{2,2,T}(2 + T(delta), 1, 1, 1))
        (; delta, es = eps_symp(A), err = catchhe(() -> check_residual("symplecticity", eps_symp(A),
            symp_threshold(DET_TOL, T, A); step = 1, orbit = 3)))
    end
    thf = Float64(theta)
    ups = [run(d) for d in above(thf)]; dns = [run(d) for d in below(thf)]
    tol = 16 * eps(T)
    okU = count(c -> abs(c.es - c.delta) <= tol && c.err isa NumericalHealthError && c.err.what == "symplecticity", ups)
    okD = count(c -> abs(c.es - c.delta) <= tol && c.err === nothing, dns)
    D = SMatrix{4,4,T}(2, 0, 0, 0, 0, 0.5, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1)   # det 1, not symplectic
    mon = ResidualMonitor(dg_FaultMap(D), DET_TOL, ORTH_TOL)
    st = reducer_init(mon, T, Val(4), Val(1))
    e = catchhe(() -> reducer_update(mon, st, 4, SVector{4,T}(0.1, 0.2, 0.3, 0.4), SMatrix{4,1,T}(1, 0, 0, 0)))
    okM = e isa NumericalHealthError && e.what == "symplecticity" && abs(e.value - 1) <= tol
    row(5, "symplecticity fault $T theta_symp=$(fmt(thf)): raise ok / quiet ok / det-1 non-symplectic named correctly",
        "$okU/$(length(ups)), $okD/$(length(dns)), $okM", "all", verdict(okU == length(ups) > 0 && okD == length(dns) > 0 && okM))
end

function ac5e(::Type{T}) where {T}
    x = SVector{2,T}(0.5, 0.25)
    W1 = SMatrix{2,1,T}(1, 0)
    ok = dg_FaultMap(SMatrix{2,2,T}(1, 0, 0, 1))
    cases = Tuple{String,Any}[]
    for what in ("det", "symplecticity", "orthogonality", "sum", "pairing"), v in (T(NaN), T(Inf))
        push!(cases, (what, catchhe(() -> check_residual(what, v, T(Inf); step = 1, orbit = 1))))
    end
    for bad in (T(NaN), T(Inf))
        mon = ResidualMonitor(dg_FaultMap(SMatrix{2,2,T}(bad, 1, 1, 1)), Inf, Inf)
        st = reducer_init(mon, T, Val(2), Val(1))
        push!(cases, ("det", catchhe(() -> reducer_update(mon, st, 1, x, W1))))
    end
    for bad in (T(NaN), T(Inf))
        mon = ResidualMonitor(ok, Inf, Inf)
        st = reducer_init(mon, T, Val(2), Val(2))
        Q = SMatrix{2,2,T}(bad, 0, 0, 1)
        push!(cases, ("orthogonality", catchhe(() -> reducer_update(mon, st, 1, x, Q))))
    end
    for bad in (T(NaN), T(Inf))
        mon = ResidualMonitor(ok, Inf, Inf)
        st = reducer_init(mon, T, Val(2), Val(2))
        lam = SVector{2,T}(bad, -0.5)
        push!(cases, ("sum", catchhe(() -> monitor_exponents(st, lam, T(6.8), 1; pairing_tol = Inf, tau = 1, M_max = T(2), step = 1, orbit = 1))))
    end
    nok = count(c -> c[2] isa NumericalHealthError && c[2].what == c[1], cases)
    unclipped = all(isfinite, (eps_orth(SMatrix{2,2,T}(T(NaN), 0, 0, 1)), eps_sum(SVector{2,T}(T(NaN), 0)),
        eps_pair(SVector{2,T}(T(Inf), 0)), eps_det(SMatrix{2,2,T}(T(NaN), 1, 1, 1)))) # must be false: NaN stays NaN
    row(5, "non-finite $T (threshold Inf): correctly named errors / never-clipped (residuals stay non-finite)",
        "$nok/$(length(cases)), $(!unclipped)", "all", verdict(nok == length(cases) && !unclipped))
end

# ---------------------------------------------------------------- AC6 and AC7 on the panel (threaded over K x T)

# AC6: 50 seeded starts x 1e4 steps; per-start max eps_det/(eps M). Same starts for both element types.
function ac6_run(::Type{T}, iK, nstarts, nsteps) where {T}
    sys = StandardMap(T(PANEL[iK]))
    rng = StableRNG(SEED + 600 + iK)
    per = zeros(nstarts); viol = 0; errs = 0
    for s in 1:nstarts
        z = SVector{2,T}(T(2pi * rand(rng)), T(2pi * rand(rng)))
        for _ in 1:nsteps
            M = tangent(sys, z)
            e = eps_det(M)
            per[s] = max(per[s], Float64(e / (eps(T) * maxabs(M))))
            e <= det_threshold(DET_TOL, T, M) || (viol += 1)
            z = forward(sys, z)
        end
    end
    return (; per, viol)
end

function ac7_run(::Type{T}, iK, n) where {T}
    sys = StandardMap(T(PANEL[iK]))
    rng = StableRNG(SEED + 700 + iK)
    x0 = SVector{2,T}(T(2pi * rand(rng)), T(2pi * rand(rng)))
    t = @elapsed r = qr_estimator(sys, z -> forward(sys, z), x0, n)
    return (; r..., seconds = t)
end

function panel(f)   # f(T, iK) for every (T, iK), threaded; returns Dict
    jobs = [(T, iK) for T in TYPES for iK in eachindex(PANEL)]
    out = Vector{Any}(undef, length(jobs))
    Threads.@threads :dynamic for j in eachindex(jobs)
        out[j] = f(jobs[j]...)
    end
    return Dict(jobs[j] => out[j] for j in eachindex(jobs))
end

# ---------------------------------------------------------------- AC8: 256-bit recomputation

function ac8a()   # package-suite item re-measured: K = 5, 1e3 steps, stored entries
    sys = StandardMap(5.0)
    rng = StableRNG(SEED + 8)
    z = SVector(2pi * rand(rng), 2pi * rand(rng))
    worst = setprecision(BigFloat, 256) do
        w = 0.0
        for _ in 1:1000
            M = tangent(sys, z)
            Mb = BigFloat.(M)
            ref = abs(Mb[1, 1] * Mb[2, 2] - Mb[1, 2] * Mb[2, 1] - 1)
            w = max(w, abs(Float64(ref) - eps_det(M)))
            z = forward(sys, z)
        end
        w
    end
    row(8, "K=5, 1e3 steps, stored entries: max |eps_det(Float64) - eps_det(256 bit)| / eps", worst / eps(), 1.0, verdict(worst <= eps()))
end

# Panel: same comparison, 1e5 steps per K (stage 04 long run).
function ac8b(iK, n)
    sys = StandardMap(PANEL[iK])
    rng = StableRNG(SEED + 800 + iK)
    z = SVector(2pi * rand(rng), 2pi * rand(rng))
    return setprecision(BigFloat, 256) do
        w = 0.0
        for _ in 1:n
            M = tangent(sys, z)
            Mb = BigFloat.(M)
            w = max(w, abs(Float64(abs(Mb[1, 1] * Mb[2, 2] - Mb[1, 2] * Mb[2, 1] - 1)) - eps_det(M)))
            z = forward(sys, z)
        end
        w
    end
end

# 256-bit estimator on matrices recomputed at the stored binary64 orbit points (K is the binary64 value).
# Reports the arbitrary-precision sum residual and the binary64 estimator's distance to it.
function ac8c(iK, n)
    Kf = PANEL[iK]
    sys = StandardMap(Kf)
    rng = StableRNG(SEED + 900 + iK)
    z = SVector(2pi * rand(rng), 2pi * rand(rng))
    Q64 = Matrix{Float64}(I, 2, 2)
    a1 = CompensatedSum{Float64}(); a2 = CompensatedSum{Float64}()
    return setprecision(BigFloat, 256) do
        Kb = BigFloat(Kf)
        Qb = Matrix{BigFloat}(I, 2, 2)
        b1 = CompensatedSum{BigFloat}(); b2 = CompensatedSum{BigFloat}()
        Gmax = BigFloat(0); Mmax = BigFloat(0); sum_big = BigFloat(0); sum_64 = 0.0; dlam = 0.0; ratio_big = BigFloat(0)
        for k in 1:n
            c = Kb * cos(BigFloat(z[1]))
            Mb = [1+c 1; c 1]
            M64 = tangent(sys, z)
            z = forward(sys, z)
            F = qr(Mb * Qb); Qb = Matrix(F.Q)
            r1 = abs(F.R[1, 1]); r2 = abs(F.R[2, 2])
            b1 = add(b1, log(r1)); b2 = add(b2, log(r2))
            F64 = qr(Matrix(M64) * Q64); Q64 = Matrix(F64.Q)
            a1 = add(a1, log(abs(F64.R[1, 1]))); a2 = add(a2, log(abs(F64.R[2, 2])))
            l64 = SVector(total(a1) / k, total(a2) / k)
            lb = (total(b1) / k, total(b2) / k)
            S = r1^2 + F.R[1, 2]^2 + r2^2; d = r1 * r2
            Gmax = max(Gmax, (S + sqrt(S^2 - 4 * d^2)) / 2 / d)
            Mmax = max(Mmax, maximum(abs, Mb))
            sb = abs(lb[1] + lb[2])
            sum_big = max(sum_big, sb)
            ratio_big = max(ratio_big, sb / (eps(BigFloat) * (Gmax + Mmax)))
            sum_64 = max(sum_64, eps_sum(l64))
            dlam = max(dlam, abs(l64[1] - Float64(lb[1])))
        end
        (; sum_big = Float64(sum_big), ratio_big = Float64(ratio_big), sum_64, dlam)
    end
end

# ---------------------------------------------------------------- AC9: zero allocation

function alloc_orbit(sys, mon, x0, W0)
    propagate_orbit(sys, x0, W0, 100, 1, (mon,))
    return @allocated propagate_orbit(sys, x0, W0, 100, 1, (mon,))
end
function alloc_update(mon, st, x, W)
    reducer_update(mon, st, 1, x, W)
    return @allocated begin
        r = reducer_update(mon, st, 2, x, W)
        r
    end
end
function alloc_exponents(st, lam, G, M)
    monitor_exponents(st, lam, G, 1; pairing_tol = PAIR_TOL, tau = 1, M_max = M, step = 1, orbit = 1)
    return @allocated begin
        r = monitor_exponents(st, lam, G, 2; pairing_tol = PAIR_TOL, tau = 1, M_max = M, step = 2, orbit = 1)
        r
    end
end
function ac9(::Type{T}) where {T}
    sys = StandardMap(T(2))
    mon = ResidualMonitor(sys, DET_TOL, ORTH_TOL)
    x0 = SVector{2,T}(1, 2)
    st2 = reducer_init(mon, T, Val(2), Val(2))
    W2 = SMatrix{2,2,T}(1, 0, 0, 1)
    lam = SVector{2,T}(0.5, -0.5)
    a = (alloc_orbit(sys, mon, x0, SMatrix{2,1,T}(1, 0)), alloc_update(mon, st2, x0, W2), alloc_exponents(st2, lam, T(2), T(2)))
    row(9, "$T bytes allocated (propagate_orbit 100 steps, reducer_update, monitor_exponents)", join(a, "/"), "0/0/0", verdict(all(==(0), a)))
    ns = @belapsed reducer_update($mon, $st2, 2, $x0, $W2) evals = 100 samples = 200
    row(9, "$T reducer_update time per call, ns (loaded machine; informational)", ns * 1e9, "none", "INFO")
end

# ---------------------------------------------------------------- AC10 / AC11 helpers

bits(v) = reinterpret(UInt64, collect(Float64, v))

function ac10_digest()
    rng = StableRNG(SEED + 10)
    sys = StandardMap(2.0)
    mon = ResidualMonitor(sys, DET_TOL, ORTH_TOL)
    x0s = [SVector(2pi * rand(rng), 2pi * rand(rng)) for _ in 1:64]
    W0s = [(a = 2pi * rand(rng); SMatrix{2,1,Float64}(cos(a), sin(a))) for _ in 1:64]
    res = propagate(sys, x0s, 100; k = 1, tau = 1, seed = SEED + 10, reducers = (mon,), W0 = W0s)
    all64 = UInt64[]
    for j in 1:64
        f = res.finals[j][1]
        append!(all64, bits(res.x[j])); append!(all64, bits(res.W[j]))
        append!(all64, bits([f.det_max, f.symp_max, f.orth_max, f.sum_max, f.pair_max]))
    end
    ser = map(1:64) do j
        x, W, st = propagate_orbit(sys, x0s[j], W0s[j], 100, 1, (mon,); orbit = j)
        (x, W, (reducer_final(mon, st[1]),))
    end
    same_serial = all(j -> res.x[j] === ser[j][1] && res.W[j] === ser[j][2] && res.finals[j] === ser[j][3], 1:64)
    return (; digest = hash(all64), same_serial, detmax = maximum(res.finals[j][1].det_max for j in 1:64))
end

# orbit i starts at (i, 0); tangent has det 1 + 2^-20 only at orbits whose first coordinate equals a target
function ac10_fault(targets)
    delta = 2.0^-20
    healthy = SMatrix{2,2,Float64}(1, 0, 0, 1)
    sys = dg_FaultMap(healthy, SMatrix{2,2,Float64}(1 + delta, 0, 0, 1), targets[1], targets[2])
    mon = ResidualMonitor(sys, DET_TOL, ORTH_TOL)
    x0s = [SVector(Float64(i), 0.0) for i in 1:64]
    W0s = [SMatrix{2,1,Float64}(1, 0) for _ in 1:64]
    e = catchhe(() -> propagate(sys, x0s, 10; k = 1, tau = 1, seed = SEED + 10, reducers = (mon,), W0 = W0s))
    return e === nothing ? 0 : e.orbit
end
const FAULT_CASES = ((33.0, 33.0), (17.0, 40.0), (40.0, 17.0))

if CHILD
    d = ac10_digest()
    println("CHILD digest=", d.digest, " threads=", Threads.nthreads(), " faults=", join((ac10_fault(t) for t in FAULT_CASES), ","))
    exit(0)
end

function ac10()
    d = ac10_digest()
    mine = (d.digest, Tuple(ac10_fault(t) for t in FAULT_CASES))
    cmd = `$(Base.julia_cmd()) --project=$(joinpath(ROOT, "verify")) -t 1 $(@__FILE__) --child`
    out = read(cmd, String)
    m = match(r"CHILD digest=(\d+) threads=(\d+) faults=([\d,]+)", out)
    m === nothing && (row(10, "child run (1 thread) failed to report", out, "-", "FAIL"); return)
    theirs = (parse(UInt64, m[1]), Tuple(parse.(Int, split(m[3], ","))))
    nt = Threads.nthreads()
    row(10, "64 orbits: x, W, residual maxima bitwise equal, $(nt) threads vs $(m[2]) thread(s) (max eps_det in run = $(fmt(d.detmax)))",
        mine[1] == theirs[1], true, verdict(mine[1] == theirs[1] && d.same_serial))
    row(10, "64 orbits: propagate bitwise equal to serial propagate_orbit map", d.same_serial, true, verdict(d.same_serial))
    exp = (33, 17, 17)
    row(10, "fault orbit ids reported ($(nt) threads / 1 thread) for targets (33), (17,40), (40,17)",
        "$(join(mine[2], ","))/$(join(theirs[2], ","))", "33,17,17", verdict(mine[2] == theirs[2] == exp))
end

function ac11(starts_by_K)
    # item 1 through the kernel (CPU backend)
    for T in TYPES
        mats = ((1, 0, 1, 1), (2, 1, 1, 1), (0, -1, 1, 0), (1, 1, -1, 0))
        mx = zero(T)
        for ent in mats
            sys = LinearSymplecticMap(SMatrix{2,2,T}(ent...))
            mon = ResidualMonitor(sys, DET_TOL, ORTH_TOL)
            res = propagate(sys, [SVector{2,T}(0, 0)], 1000; k = 1, tau = 1, seed = SEED + 11,
                reducers = (dg_Renorm(), mon), backend = :cpu, W0 = [SMatrix{2,1,T}(1, 0)])
            f = res.finals[1][2]
            mx = max(mx, f.det_max, f.symp_max, f.orth_max)
        end
        row(11, "CPU kernel item 1 $T: max residual over 4 integer maps x1000 steps", Float64(mx), 0.0, verdict(mx == 0))
    end
    # item 5 through the kernel
    A0 = SMatrix{2,2,Float64}(2, 1, 1, 1)
    theta = Float64(det_threshold(DET_TOL, Float64, A0))
    run5(delta, N) = begin
        sys = dg_FaultMap(SMatrix{2,2,Float64}(2 + delta, 1, 1, 1))
        mon = ResidualMonitor(sys, DET_TOL, ORTH_TOL)
        catchhe(() -> propagate(sys, [SVector(0.5, 0.25)], N; k = 1, tau = 1, seed = SEED + 11,
            reducers = (dg_Renorm(), mon), backend = :cpu, W0 = [SMatrix{2,1,Float64}(1, 0)]))
    end
    ups = [run5(d, 1) for d in above(theta)]; dns = [run5(d, 1000) for d in below(theta)]
    okU = count(e -> e isa NumericalHealthError && e.what == "det", ups); okD = count(isnothing, dns)
    row(11, "CPU kernel item 5 (det fault, theta=$(fmt(theta))): raise ok / quiet ok", "$okU/$(length(ups)), $okD/$(length(dns))", "all",
        verdict(okU == length(ups) > 0 && okD == length(dns) > 0))
    # item 6 through the kernel, and bitwise agreement with the direct loop on the same starts
    worst_kernel = 0.0; worst_gap = 0.0; errs = 0
    for (iK, K) in enumerate(PANEL)
        sys = StandardMap(K)
        mon = ResidualMonitor(sys, DET_TOL, ORTH_TOL)
        rng = StableRNG(SEED + 1100 + iK)
        x0s = [SVector(2pi * rand(rng), 2pi * rand(rng)) for _ in 1:50]
        W0s = [(a = 2pi * rand(rng); SMatrix{2,1,Float64}(cos(a), sin(a))) for _ in 1:50]
        res = nothing
        e = catchhe(() -> (res = propagate(sys, x0s, 10_000; k = 1, tau = 1, seed = SEED + 11, reducers = (dg_Renorm(), mon),
                                           backend = :cpu, W0 = W0s)))
        e === nothing || (errs += 1; continue)
        for j in 1:50
            z = x0s[j]; dmax = 0.0
            for _ in 1:10_000
                dmax = max(dmax, eps_det(tangent(sys, z))); z = forward(sys, z)
            end
            kmax = res.finals[j][2].det_max
            worst_kernel = max(worst_kernel, kmax / eps())
            worst_gap = max(worst_gap, abs(kmax - dmax))
        end
    end
    row(11, "CPU kernel item 6 (panel, 50 starts x 1e4 steps): health errors; max det_max/eps; max |kernel - direct loop| det_max",
        "$errs; $(round(worst_kernel; digits = 3)); $(fmt(worst_gap))", "0; -; 0", verdict(errs == 0 && worst_gap == 0))
    # arbitrary precision on a GPU request falls back with the recorded reason
    c = resolve_backend(:gpu, BigFloat)
    row(11, "resolve_backend(:gpu, BigFloat): backend / reason", "$(c.backend) / $(c.reason)", "cpu / " * FALLBACK_REASONS[3],
        verdict(c.backend === :cpu && c.reason == FALLBACK_REASONS[3]))
    row(11, "GPU-vs-CPU verdicts and residual differences on the GTX 1050 (items 1, 5, 6)", "see verify/gpu.jl", "<= threshold in force", "STAGE04-GPU")
end

# ---------------------------------------------------------------- figures

function figures(six, seven)
    mkpath(OUT)
    wong = CM.Makie.wong_colors()
    # Figure 1: per-start max eps_det/(eps M) (AC6) and the 1e6-step orbit (AC7)
    fig = CM.Figure(size = (900, 480))
    ax = CM.Axis(fig[1, 1]; xlabel = "K", ylabel = "max over steps of eps_det / (eps(T) M)", xticks = (1:5, string.(collect(PANEL))),
        title = "Standard map: realised determinant residual against B1 bound 1.5")
    for (it, T) in enumerate(TYPES)
        off = it == 1 ? -0.17 : 0.17
        xs = Float64[]; ys = Float64[]
        for iK in eachindex(PANEL)
            per = six[(T, iK)].per
            append!(xs, iK .+ off .+ range(-0.12, 0.12; length = length(per))); append!(ys, per)
        end
        CM.scatter!(ax, xs, ys; color = wong[it], markersize = 5, label = "$T, per start (50 starts, 1e4 steps)")
        CM.scatter!(ax, [iK + off for iK in eachindex(PANEL)], [seven[(T, iK)].detratio for iK in eachindex(PANEL)];
            color = wong[it], marker = :star5, markersize = 16, strokewidth = 1, label = "$T, one orbit of $(N_LONG) steps")
    end
    CM.hlines!(ax, [1.5]; color = :black, linestyle = :dash, label = "bound 1.5 (B1)")
    CM.ylims!(ax, 0, 1.7)
    CM.Legend(fig[1, 2], ax)
    CM.save(joinpath(OUT, "diagnostics-det-residual.png"), fig; px_per_unit = 2)
    # Figure 2: eps_sum/(eps (G_max/tau + M_max)) against n (AC7), running maximum per log bin
    fig2 = CM.Figure(size = (1000, 520))
    ax2 = CM.Axis(fig2[1, 1]; xlabel = "n (map steps)", ylabel = "max eps_sum(n) / (eps(T) (G_max/tau + M_max)) per bin",
        xscale = log10, yscale = log10, title = "Sum-of-exponents residual, reference QR estimator, tau = 1")
    for (it, T) in enumerate(TYPES), iK in eachindex(PANEL)
        b = seven[(T, iK)].bins
        used = findall(>(0), b)
        xs = [10.0^((i - 1) / BINS_PER_DECADE) for i in used]
        CM.lines!(ax2, xs, max.(b[used], 1e-3); color = wong[iK], linestyle = it == 1 ? :solid : :dash,
            label = "K=$(PANEL[iK]) $T")
    end
    CM.hlines!(ax2, [4.0]; color = :black, linestyle = :dash, label = "bound 4 (AC7)")
    CM.hlines!(ax2, [PAIR_TOL]; color = :gray, linestyle = :dot, label = "key pairing_residual_tol = 16")
    CM.Legend(fig2[1, 2], ax2)
    CM.save(joinpath(OUT, "diagnostics-sum-residual.png"), fig2; px_per_unit = 2)
end

# ---------------------------------------------------------------- main

function main()
    println("# diagnostics verify, threads=$(Threads.nthreads()), Julia $(VERSION), seed=$SEED, started $(time())")
    println("# format: AC<n> | quantity | value | bound | verdict")
    for T in TYPES; ac1(T); end
    for T in TYPES; ac2(T, N_LONG); end
    ac3(100_000)
    for T in TYPES; ac4(T); end
    for T in TYPES; ac5a(T); end
    ac5b(Float64, 2.0^-30); ac5b(Double64, 2.0^-60)
    for T in TYPES; ac5c(T); end
    for T in TYPES; ac5d(T); end
    for T in TYPES; ac5e(T); end
    t6 = @elapsed six = panel((T, iK) -> ac6_run(T, iK, 50, 10_000))
    for T in TYPES
        for iK in eachindex(PANEL)
            r = six[(T, iK)]
            mx = maximum(r.per)
            row(6, "$T K=$(PANEL[iK]) 50 starts x 1e4 steps: max eps_det/(eps M), threshold violations=$(r.viol)", mx, 1.5, verdict(mx <= 1.5 && r.viol == 0))
        end
    end
    println("# AC6 panel wall time (s): ", round(t6; digits = 1))
    seven = panel((T, iK) -> ac7_run(T, iK, N_LONG))
    for T in TYPES, iK in eachindex(PANEL)
        r = seven[(T, iK)]
        row(7, "$T K=$(PANEL[iK]) n=$N_LONG: max_n eps_sum/(eps (G_max/tau + M_max)) [G_max=$(round(Float64(r.Gmax); sigdigits = 4)), M_max=$(round(Float64(r.Mmax); sigdigits = 4)), snapshots>4: $(r.viol4), >16: $(r.viol16), lambda_1=$(round(Float64(r.lam[1]); digits = 4)), $(round(r.seconds; digits = 1)) s]",
            r.ratio_max, 4.0, verdict(r.ratio_max <= 4))
        row(6, "$T K=$(PANEL[iK]) one orbit of n=$N_LONG (during AC7): max eps_det/(eps M)", r.detratio, 1.5, verdict(r.detratio <= 1.5))
    end
    row(7, "ChaosTools full-spectrum sum of exponents next to eps_sum (comparison only)", "n/a", "n/a", "DEFERRED lyap-benettin 04")
    ac8a()
    t8 = @elapsed (b = panel_big(ac8b, N_BIG))
    for iK in eachindex(PANEL)
        row(8, "K=$(PANEL[iK]) $N_BIG steps, stored entries vs 256 bit: max |diff| / eps", b[iK] / eps(), 1.0, verdict(b[iK] <= eps()))
    end
    c = panel_big(ac8c, N_BIGSUM)
    for iK in eachindex(PANEL)
        r = c[iK]
        row(8, "K=$(PANEL[iK]) $N_BIGSUM steps, 256-bit QR estimator on recomputed entries: max eps_sum/(eps_256 (G+M)) [binary64 eps_sum max = $(fmt(r.sum_64)), max|lambda_1(64) - lambda_1(256)| = $(fmt(r.dlam))]",
            r.ratio_big, PAIR_TOL, verdict(r.ratio_big <= PAIR_TOL))
    end
    println("# AC8 panel wall time (s): ", round(t8; digits = 1))
    for T in TYPES; ac9(T); end
    ac10()
    ac11(nothing)
    figures(six, seven)
    println("# figures written to ", OUT)
end

function panel_big(f, n)
    out = Vector{Any}(undef, length(PANEL))
    Threads.@threads :dynamic for iK in eachindex(PANEL)
        out[iK] = f(iK, n)
    end
    return out
end

main()
