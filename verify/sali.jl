# Stage 04 verification of module sali (reducer SALI, alignment_indices, sali_pair, sali_from_log, check_frame; layer Indicators).
# Run from the repo root:  julia -t 2 --project=verify verify/sali.jl
# Spec: workspaces/module-build/stages/01-spec/output/sali-spec.md (criteria 1 to 17, Amendment D67). Criterion 18 (CPU vs GPU) is in verify/gpu.jl, tag "sali".
# Section "suite": the package suite at full size (SL_FULL_RUN = true: 20 starts for every Double64 panel, checkpoints at every multiple of 10).
# Other sections: one printed line per criterion  AC<n> | quantity | value | bound | PASS/FAIL  with the MEASURED number, reusing the helpers of
# test/sali/sl_helpers.jl (one home for the tolerance formulas B_meas, B_SALI, B_L, Lambda_n, Psi_n and the BigFloat reference). Values are
# "worst ratio error / bound" (pass when <= 1) unless stated. The cases use verify seeds (SL_SEED + 56000 + 50 ac + sub), distinct from the package
# tests (SL_SEED + 31000 ..); the K = 0.5 regular class is screened here again (Benettin X(1e5) <= 5e-4, as L1 S7) with a verify seed.
# SL_ONLY="AC16,figure" runs only the named sections (debugging aid) and SL_LOG names the log file; the default runs everything.
# SL_CHILD=1 is the child-process mode of the thread-invariance check of AC14.
using Printf, Test, Random, LinearAlgebra, Statistics
using Chaoslab, StaticArrays, DoubleFloats, StableRNGs

const SL_FULL_RUN = true
const ROOT = dirname(@__DIR__)
const TESTDIR = joinpath(ROOT, "test", "sali")
include(joinpath(TESTDIR, "sl_helpers.jl"))

vseed(ac, sub = 0) = SL_SEED + 56000 + 50 * ac + sub

# Thread-count invariance (AC14b): the same series hashed in this process and in child processes with other thread counts.
function sl_thread_hash()
    m, n = 4096, 1000
    h = UInt(0)
    for tau in (1, 4)
        sys = StandardMap(sl_K(Float64, :two_pi))
        x0 = [sl_cast(Float64, s) for s in sl_starts(StableRNG(vseed(14, 3)), m)]
        res = propagate(sys, x0, n; k = 2, tau = tau, seed = vseed(14, 3), reducers = (SALI(),))
        flat = Float64[]
        foreach(z -> append!(flat, z), res.x)
        foreach(W -> append!(flat, W), res.W)
        for f in res.finals
            g = f[1]
            push!(flat, g.sali, g.L)
            append!(flat, g.sum)
            append!(flat, g.comp)
        end
        h = hash(reinterpret(UInt64, flat), h)
    end
    return h
end
if get(ENV, "SL_CHILD", "") == "1"
    println("HASH ", sl_thread_hash(), " THREADS ", Threads.nthreads())
    exit(0)
end

import ChaosTools as CT
import CairoMakie as CM

const OUT = joinpath(ROOT, "workspaces", "module-build", "stages", "04-verify", "output")
mkpath(OUT)
const ONLY = filter(!isempty, split(get(ENV, "SL_ONLY", ""), ","))
const LOG = open(get(ENV, "SL_LOG", joinpath(OUT, "sali-verify-log.txt")), "w")
const VERDICTS = String[]

fmt(x::AbstractFloat) = @sprintf("%.4g", x)
fmt(x) = string(x)
function say(s)
    println(s); println(LOG, s); flush(stdout); flush(LOG)
end
function emit(ac, quantity, value, bound, verdict)
    push!(VERDICTS, verdict)
    say(string(ac, " | ", quantity, " | ", fmt(value), " | ", bound, " | ", verdict))
end
pf(ok::Bool) = ok ? "PASS" : "FAIL"
ratio_line(ac, quantity, r) = emit(ac, quantity, r, "<= 1", pf(r <= 1))
info(ac, quantity, value) = emit(ac, quantity, value, "reported", "INFO")

function section(f, name)
    isempty(ONLY) || name in ONLY || return
    t0 = time()
    try
        f()
    catch e
        say(string(name, " | ERROR | ", sprint(showerror, e)[1:min(end, 300)], " | - | FAIL"))
        push!(VERDICTS, "FAIL")
    end
    say(@sprintf("# section %s took %.1f s", name, time() - t0))
end

const T0 = time()
const DATA = Dict{Any, Any}()
const EPS64 = sl_eps(Float64)
const EPSDD = sl_eps(Double64)
say("# sali verification, julia $(VERSION), threads $(Threads.nthreads()), base seed $SL_SEED (verify offset 56000), time $(round(time())) (unix)")
say("# timings in this log were taken on a machine that was not otherwise idle")

# ---------------------------------------------------------------------------------------------------
# Cases with verify seeds (the package helpers sl_case / sl_rows use the package seeds)
const VC = Dict{Any, Any}()

# Regular class K = 0.5 (L1 S7): candidates c = 1, 2, .. draw (x, p) from a StableRNG, orbit id c; regular when Benettin X(1e5) <= 5e-4.
function vregular(m::Int)
    return get!(VC, (:regular, m)) do
        seed = vseed(0, 4)
        sys = StandardMap(0.5)
        rng = StableRNG(seed)
        cands = SVector{2, Float64}[]
        ok = Int[]
        while length(ok) < m
            lo = length(cands) + 1
            append!(cands, sl_starts(rng, 60))
            W0 = [initial_deviations(seed, c, Val(2), Val(1), Float64) for c in lo:length(cands)]
            res = propagate(sys, cands[lo:end], 100_000; k = 1, tau = 1, seed = seed, reducers = (Benettin(),), W0 = W0)
            append!(ok, (lo - 1) .+ findall(<=(5.0e-4), [f[1].lambda[1] for f in res.finals]))
        end
        keep = ok[1:m]
        (ids = keep, starts = cands[keep], rejected = last(keep) - m, seed = seed)
    end
end

function vcase(label::Symbol, ::Type{T}, m::Int) where {T}
    return get!(VC, (:case, label, T, m)) do
        K = sl_K(T, label)
        if label === :half
            r = vregular(50)
            ids, starts, fseed = r.ids[1:m], r.starts[1:m], r.seed
        else
            fseed = vseed(0, findfirst(==(label), (:two_pi, :ten, :five)))
            ids = collect(1:m)
            starts = sl_starts(StableRNG(fseed), m)
        end
        (sys = StandardMap(K), K = K, ids = ids, x0 = [sl_cast(T, s) for s in starts], fseed = fseed,
         W0 = [sl_frame(fseed, id, T) for id in ids])
    end
end

# BigFloat reference, Psi_n and Lambda_n along the orbit of every start at the steps `cps` (A15).
function vref(label::Symbol, ::Type{T}, m::Int, N::Int, cps::Vector{Int}) where {T}
    return get!(VC, (:ref, label, T, m, N, cps)) do
        c = vcase(label, T, m)
        map(1:m) do j
            pts, _ = lb_pts(c.sys, c.x0[j], N)
            ref = sl_ref(sl_mats_std(c.K, pts), c.W0[j], N; keep = cps)
            Ms = lb_tangents(c.sys, c.x0[j], N)
            (ref = ref, psi = sl_psis(Ms, c.W0[j], cps), lam = [sl_Lambda(ref.lnsali, n) for n in cps])
        end
    end
end

# One row per (start, checkpoint) with the run values, the reference and the bounds (as sl_rows of the package helpers).
function vrows(label::Symbol, ::Type{T}, tau::Int, m::Int, N::Int; cps::Vector{Int} = collect(10:10:N)) where {T}
    return get!(VC, (:rows, label, T, tau, m, N, cps)) do
        O = vref(label, T, m, N, cps)
        c = vcase(label, T, m)
        hs = filter(n -> n % tau == 0, cps)
        R = sl_run(SALI(), c.sys, c.x0, c.W0, N, tau; steps = hs)
        epsT = sl_eps(T)
        rows = NamedTuple[]
        for j in 1:m, s in R[j].snaps
            i = findfirst(==(s.n), cps)
            ref, (p1, p2), lam = O[j].ref, O[j].psi[i], O[j].lam[i]
            at = sl_at(ref, s.n)
            saliref = Float64(at.S)
            push!(rows, (j = j, n = s.n, sali = s.sali, L = s.L, S1 = s.sum[1] + s.comp[1], S2 = s.sum[2] + s.comp[2],
                         errS = sl_err(s.sali, at.S), errL = sl_err(s.L, at.L), saliref = saliref, Lref = Float64(at.L),
                         BS = sl_BSALI(epsT, lam, saliref, p1, p2), BL = sl_BL(epsT, p1, p2), psi1 = p1, psi2 = p2))
        end
        rows
    end
end

sl_panel_v(::Type{T}) where {T} = ((:half, 10_000, 20), (:two_pi, 1_000, 20), (:ten, 1_000, 20))
mean_se(v) = (mean(v), std(v) / sqrt(length(v)))
zdiff(a, b) = abs(mean(a) - mean(b)) / sqrt(var(a) / length(a) + var(b) / length(b))

# ---------------------------------------------------------------------------------------------------
# The package suite at full size
section("suite") do
    local ts
    t = @elapsed try
        ts = @testset verbose = true "sali verify" begin
            include(joinpath(TESTDIR, "runtests.jl"))
        end
    catch e
        e isa Test.TestSetException || rethrow()
        ts = nothing
        say("# suite: TestSetException: $(e.pass) passed, $(e.fail) failed, $(e.error) errors")
    end
    ok = ts !== nothing
    c = ok ? Test.get_test_counts(ts) : nothing
    npass = ok ? c.passes + c.cumulative_passes : -1
    nbad = ok ? c.fails + c.errors + c.cumulative_fails + c.cumulative_errors : -1
    emit("SUITE", "package suite at full size (SL_FULL_RUN), assertions passed (fail + error = $nbad); wall time $(round(t; digits = 1)) s", npass, "0 fail, 0 error", pf(ok && nbad == 0 && npass > 0))
end

# ---------------------------------------------------------------------------------------------------
# AC1: algebra of the definition
function ac1(T)
    rng = StableRNG(vseed(1))
    epsT, meas = sl_eps(T), sl_meas(T)
    emax = T === Float64 ? 200 : 60
    wa = wb = wc = wd = 0.0
    bad_min = bad_bits = 0
    for i in 1:11_000
        dim = i <= 10_000 ? 2 : 4
        e1, e2 = rand(rng, -emax:emax), rand(rng, -emax:emax)
        a = sl_randvec(rng, T, dim, e1)
        if iseven(i)
            b = sl_randvec(rng, T, dim, e2)
        else
            b = a + T(10.0^(-1 - 13 * rand(rng))) * sl_randvec(rng, T, dim, e1)
        end
        rand(rng, Bool) && (b = -b)
        dm, dp = alignment_indices(a, b)
        s = sali_pair(a, b)
        bad_min += s == min(dm, dp) ? 0 : 1
        dmx, dpx, sx = sl_sali_exact(a, b)
        wa = max(wa, sl_err(dm, dmx) / meas, sl_err(dp, dpx) / meas)
        wb = max(wb, Float64(abs(dm^2 + dp^2 - 4)) / (16 * epsT))
        wc = max(wc, Float64(abs(dm * dp / 2 - s * max(dm, dp) / 2)) / (8 * epsT))
        wd = max(wd, sl_err(s, sx) / meas)
        bad_bits += (s === sali_pair(-a, b) && s === sali_pair(a, -b) && s === sali_pair(b, a) && s === sali_pair(-b, -a)) ? 0 : 1
    end
    wcon = wconb = 0.0
    nfar = 0
    for m in 1:14, sgn in (1, -1)
        ang = T(10.0)^(-m)
        u1 = SVector{2, T}(cos(T(0.3)), sin(T(0.3)))
        u2 = sgn * SVector{2, T}(cos(T(0.3) + ang), sin(T(0.3) + ang))
        dm, dp = alignment_indices(u1, u2)
        s = sali_pair(u1, u2)
        wcon = max(wcon, sl_err(s, sl_sali_exact(u1, u2)[3]) / meas)
        wconb = max(wconb, Float64(abs(dm^2 + dp^2 - 4)) / (16 * epsT))
        nfar += s < 3 * 10.0^(-m) ? 0 : 1
    end
    nzero = count(1:50) do _
        a = sl_randvec(rng, T, 2, rand(rng, -emax:emax))
        sali_pair(a, a) == 0 && sali_pair(a, -a) == 0 && sali_pair(-a, a) == 0
    end
    wg = 0.0
    for (p, q) in ((SVector{2, T}(1, 0), SVector{2, T}(0, 1)), (SVector{2, T}(3, 4), SVector{2, T}(-4, 3)),
                   (SVector{2, T}(5, 12), SVector{2, T}(-12, 5)), (SVector{2, T}(0, 7), SVector{2, T}(7, 0))), e in (-100, 0, 40)
        wg = max(wg, Float64(abs(sali_pair(ldexp.(p, e), ldexp.(q, e)) - sqrt(T(2)))) / (4 * epsT))
    end
    u = SVector{2, T}(cos(T(1)), sin(T(1)))
    neg = Float64(alignment_indices(u, -u)[1])
    return (; wa, wb, wc, wd, bad_min, bad_bits, wcon, wconb, nfar, nzero, wg, neg, sneg = sali_pair(u, -u))
end

section("AC1") do
    for T in SL_TYPES
        r = ac1(T)
        ratio_line("AC1", "(a) max |d-, d+ - BigFloat(256)| / B_meas, 11000 random and near-parallel pairs (d = 2 and 4), $T", r.wa)
        emit("AC1", "(a) SALI == min(d-, d+) failures over the 11000 pairs, $T", r.bad_min, "0", pf(r.bad_min == 0))
        ratio_line("AC1", "(b) max |d-^2 + d+^2 - 4| / (16 eps) (Amendment D67), $T", r.wb)
        info("AC1", "(b) same ratio against the original 8 eps of the spec text (above 1 means the amendment was needed), $T", 2 * r.wb)
        ratio_line("AC1", "(b) same for the 28 constructed pairs (rotation by 10^-m, parallel and antiparallel), $T", r.wconb)
        ratio_line("AC1", "(c) max |d- d+/2 - SALI max(d-, d+)/2| / (8 eps), $T", r.wc)
        ratio_line("AC1", "(d) max |SALI - SALI_exact| / B_meas, random pairs, $T", r.wd)
        ratio_line("AC1", "(d) same for the 28 constructed pairs, $T", r.wcon)
        emit("AC1", "(e) pairs whose SALI changes bitwise under a sign change of either vector or an exchange, $T", r.bad_bits, "0", pf(r.bad_bits == 0))
        emit("AC1", "(f) identical and exactly opposite vectors giving exactly 0 (of 50), $T", r.nzero, "50", pf(r.nzero == 50))
        ratio_line("AC1", "(g) exactly orthogonal pairs: max |SALI - sqrt 2| / (4 eps), $T", r.wg)
        emit("AC1", "constructed pairs with SALI >= 3 10^-m (the pair is not that close), $T", r.nfar, "0", pf(r.nfar == 0))
        emit("AC1", "negative control: d- alone for an antiparallel pair (sali_pair gives $(fmt(Float64(r.sneg)))), $T", r.neg, "> 1.9", pf(r.neg > 1.9 && r.sneg == 0))
    end
end

# ---------------------------------------------------------------------------------------------------
# AC2: tau invariance
section("AC2") do
    for T in SL_TYPES, label in (:half, :two_pi)
        base = Dict((r.j, r.n) => r for r in vrows(label, T, 1, 20, 100))
        wS = wL = 0.0
        for tau in (2, 5, 10), r in vrows(label, T, tau, 20, 100)
            r.n in (10, 100) || continue
            b = base[(r.j, r.n)]
            wS = max(wS, abs(Float64(r.sali - b.sali)) / (r.BS + b.BS))
            wL = max(wL, abs(Float64(r.L - b.L)) / (r.BL + b.BL))
        end
        ratio_line("AC2", "K = $(label === :half ? "0.5 regular" : "2 pi"), $T: max |SALI_tau - SALI_1| / (B_SALI^tau + B_SALI^1), tau in {2,5,10}, n in {10,100}, 20 starts", wS)
        ratio_line("AC2", "K = $(label === :half ? "0.5 regular" : "2 pi"), $T: max |L_tau - L_1| / (B_L^tau + B_L^1), same", wL)
    end
end

# ---------------------------------------------------------------------------------------------------
# AC3: K = 0 shear (A10)
function v_shear_cf(W0, n::Integer)
    return setprecision(BigFloat, SL_BITS) do
        a1, b1, a2, b2 = sl_big(W0[1, 1]), sl_big(W0[2, 1]), sl_big(W0[1, 2]), sl_big(W0[2, 2])
        det0 = abs(a1 * b2 - a2 * b1)
        w1, w2 = (a1 + n * b1, b1), (a2 + n * b2, b2)
        (S = 2 * sin(atan(det0, abs(w1[1] * w2[1] + w1[2] * w2[2])) / 2),
         L = log(det0) - log(sqrt(w1[1]^2 + w1[2]^2)) - log(sqrt(w2[1]^2 + w2[2]^2)))
    end
end
function v_shear_lnsali(W0, N::Integer)
    a1, b1, a2, b2 = Float64(W0[1, 1]), Float64(W0[2, 1]), Float64(W0[1, 2]), Float64(W0[2, 2])
    det0 = abs(a1 * b2 - a2 * b1)
    return [log(2 * sin(atan(det0, abs((a1 + j * b1) * (a2 + j * b2) + b1 * b2)) / 2)) for j in 1:N]
end

function ac3(T)
    rng = StableRNG(vseed(3))
    epsT = sl_eps(T)
    sys = StandardMap(exact(T, 0))
    M = SMatrix{2, 2, Float64}(1, 0, 1, 1)
    wS = wL = wI = 0.0
    rep = NaN
    for tau in (1, 10), (case, nid) in (("seeded", 20), ("identity", 1)), N in (10_000, 1_000_000)
        m = case == "identity" ? 1 : (N == 1_000_000 ? 3 : (T === Double64 ? 8 : nid))
        x0 = [sl_cast(T, SVector(2pi * rand(rng), 2pi * rand(rng))) for _ in 1:m]
        W0 = case == "identity" ? [SMatrix{2, 2, T}(1, 0, 0, 1)] : [sl_frame(vseed(3), j, T) for j in 1:m]
        ns = N == 1_000_000 ? [1_000_000] : [10, 100, 10_000]
        runs = sl_run(SALI(), sys, x0, W0, N, tau; steps = ns)
        Ms = fill(M, N)
        for j in 1:m
            lns = v_shear_lnsali(W0[j], N)
            for s in runs[j].snaps
                n = s.n
                cf = v_shear_cf(W0[j], n)
                p1, p2 = lb_Psi(view(Ms, 1:n), W0[j][:, 1]), lb_Psi(view(Ms, 1:n), W0[j][:, 2])
                BS = sl_BSALI(epsT, sl_Lambda(lns, n), Float64(cf.S), p1, p2)
                wS = max(wS, sl_err(s.sali, cf.S) / BS)
                wL = max(wL, sl_err(s.L, cf.L) / sl_BL(epsT, p1, p2))
                case == "identity" && (wI = max(wI, Float64(abs(s.sali * n - 1)) / (3 / (8 * n^2) + n * BS)))     # in T: 3/(8 n^2) is below the Float64 resolution of SALI n at n = 1e6 for Double64
                if case == "seeded" && N == 10_000 && n == 10_000 && j == 1 && tau == 1
                    rep = Float64(s.sali * n^2 * abs(W0[j][2, 1] * W0[j][2, 2]) / abs(det(W0[j])))
                end
            end
        end
    end
    return (; wS, wL, wI, rep)
end

section("AC3") do
    for T in SL_TYPES
        r = ac3(T)
        ratio_line("AC3", "max |SALI - SALI_cf| / B_SALI, K = 0, seeded frame and W0 = I, tau in {1,10}, n in {10,100,1e4,1e6}, $T", r.wS)
        ratio_line("AC3", "max |L - L_cf| / B_L, same, $T", r.wL)
        ratio_line("AC3", "W0 = I: max |SALI n - 1| / (3/(8 n^2) + n B_SALI), $T", r.wI)
        info("AC3", "SALI n^2 |b1 b2| / |det W0| at n = 1e4, first seeded start (A10: 1 - (a1/b1 + a2/b2)/n), $T", r.rep)
    end
end

# ---------------------------------------------------------------------------------------------------
# AC4 to AC6: analytic cases
function v_y5_eigen()
    mu = (3 + sqrt(BigFloat(5))) / 2
    phi = (1 + sqrt(BigFloat(5))) / 2
    nv = sqrt(phi^2 + 1)
    return (mu = mu, vp = (phi / nv, 1 / nv), vm = (-1 / nv, phi / nv))
end

function ac4(T, nid)
    seed = vseed(4)
    keep = [1, 2, 4, 5, 8, 10, 20, 100, 200]
    epsT, meas = sl_eps(T), sl_meas(T)
    sys = LinearSymplecticMap(SMatrix{2, 2, T}(2, 1, 1, 1))
    x0 = [SVector{2, T}(0, 0) for _ in 1:nid]
    W0 = [sl_frame(seed, j, T) for j in 1:nid]
    Ms = lb_tangents(sys, x0[1], 200)
    refs = [sl_ref(sl_mats_const(2, 1, 1, 1), W0[j], 200; keep = keep) for j in 1:nid]
    psis = [sl_psis(Ms, W0[j], keep) for j in 1:nid]
    c = setprecision(BigFloat, SL_BITS) do
        e = v_y5_eigen()
        map(1:nid) do j
            w1 = (sl_big(W0[j][1, 1]), sl_big(W0[j][2, 1]))
            w2 = (sl_big(W0[j][1, 2]), sl_big(W0[j][2, 2]))
            c1 = (w1[1] * e.vp[1] + w1[2] * e.vp[2], w2[1] * e.vp[1] + w2[2] * e.vp[2])
            c2 = (w1[1] * e.vm[1] + w1[2] * e.vm[2], w2[1] * e.vm[1] + w2[2] * e.vm[2])
            (c1 = c1, c2 = c2, D = c1[1] * c2[2] - c2[1] * c1[2])
        end
    end
    excluded = count(cj -> abs(cj.c1[1] * cj.c1[2]) < 0.01, c)
    helper = setprecision(BigFloat, SL_BITS) do
        e = v_y5_eigen()
        maximum(1:5) do j
            n2(i) = e.mu^16 * c[j].c1[i]^2 + e.mu^(-16) * c[j].c2[i]^2
            Float64(abs(log(abs(c[j].D)) - log(n2(1)) / 2 - log(n2(2)) / 2 - sl_at(refs[j], 8).L))
        end
    end
    wa = wb = wc = 0.0
    nres = 0
    for tau in (1, 5)
        run = sl_run(SALI(), sys, x0, W0, 200, tau; steps = filter(n -> n % tau == 0, keep))
        for j in 1:nid
            cj = c[j]
            abs(cj.c1[1] * cj.c1[2]) < 0.01 && continue
            for s in run[j].snaps
                n = s.n
                at = sl_at(refs[j], n)
                p1, p2 = psis[j][findfirst(==(n), keep)]
                saliref = Float64(at.S)
                if n in (tau == 1 ? (1, 2, 4, 8) : (5, 10)) && saliref >= 10 * meas
                    nres += 1
                    wa = max(wa, sl_err(s.sali, at.S) / sl_BSALI(epsT, sl_Lambda(refs[j].lnsali, n), saliref, p1, p2))
                end
                n in (20, 100, 200) && (wb = max(wb, sl_err(s.L, at.L) / sl_BL(epsT, p1, p2)))
                if n in (20, 100)
                    asym, mu4 = setprecision(BigFloat, SL_BITS) do
                        mu = v_y5_eigen().mu
                        (2 * n * log(mu) - log(abs(cj.D)) + log(abs(cj.c1[1] * cj.c1[2])), Float64(mu^(-4n)))
                    end
                    wc = max(wc, sl_err(-s.L, asym) / (sl_BL(epsT, p1, p2) + mu4))
                end
            end
        end
    end
    return (; wa, wb, wc, nres, excluded, helper)
end

function ac5(T, K, nid)
    seed = vseed(5)
    keep = [1, 2, 4, 5, 8, 10, 20, 100, 200]
    epsT, meas = sl_eps(T), sl_meas(T)
    sys = StandardMap(exact(T, K))
    x0 = [SVector{2, T}(0, 0) for _ in 1:nid]
    W0 = [sl_frame(seed, j, T) for j in 1:nid]
    Ms = lb_tangents(sys, x0[1], 200)
    refs = [sl_ref(sl_mats_const(1 + K, 1, K, 1), W0[j], 200; keep = keep) for j in 1:nid]
    psis = [sl_psis(Ms, W0[j], keep) for j in 1:nid]
    d = setprecision(BigFloat, SL_BITS) do
        Kb = BigFloat(K)
        s = sqrt(Kb * (Kb + 4))
        mu = (2 + Kb + s) / 2
        vu, vs = (one(Kb), (s - Kb) / 2), (one(Kb), -(Kb + s) / 2)
        ai(w) = (w[2] + w[1] * (Kb + s) / 2) / s
        per = map(1:nid) do j
            w = ((sl_big(W0[j][1, 1]), sl_big(W0[j][2, 1])), (sl_big(W0[j][1, 2]), sl_big(W0[j][2, 2])))
            a = ntuple(i -> ai(w[i]), 2)
            b = ntuple(i -> w[i][1] - a[i], 2)
            (a = a, b = b, ldet = log(abs(w[1][1] * w[2][2] - w[1][2] * w[2][1])))
        end
        (mu = mu, vu = vu, vs = vs, per = per)
    end
    helper = setprecision(BigFloat, SL_BITS) do
        vv = d.vu[1] * d.vs[1] + d.vu[2] * d.vs[2]
        maximum(1:nid) do j
            n = 8
            n2(i) = d.mu^(2n) * d.per[j].a[i]^2 * (d.vu[1]^2 + d.vu[2]^2) + 2 * d.per[j].a[i] * d.per[j].b[i] * vv + d.mu^(-2n) * d.per[j].b[i]^2 * (d.vs[1]^2 + d.vs[2]^2)
            Float64(abs(d.per[j].ldet - log(n2(1)) / 2 - log(n2(2)) / 2 - sl_at(refs[j], n).L))
        end
    end
    wa = wb = 0.0
    nres = 0
    maxcoef = 0.0
    for tau in (1, 5)
        run = sl_run(SALI(), sys, x0, W0, 200, tau; steps = filter(n -> n % tau == 0, keep))
        for j in 1:nid, s in run[j].snaps
            n = s.n
            at = sl_at(refs[j], n)
            p1, p2 = psis[j][findfirst(==(n), keep)]
            saliref = Float64(at.S)
            if n in (tau == 1 ? (1, 2, 4, 8) : (5, 10)) && saliref >= 10 * meas
                nres += 1
                wa = max(wa, sl_err(s.sali, at.S) / sl_BSALI(epsT, sl_Lambda(refs[j].lnsali, n), saliref, p1, p2))
            end
            if n in (20, 100, 200)
                asym, tail, coef = setprecision(BigFloat, SL_BITS) do
                    nvu = sqrt(d.vu[1]^2 + d.vu[2]^2)
                    coef = sum(i -> abs(d.per[j].b[i] / d.per[j].a[i]), 1:2) * abs(d.vu[1] * d.vs[1] + d.vu[2] * d.vs[2]) / nvu^2
                    (2 * n * log(d.mu) - d.per[j].ldet + sum(i -> log(abs(d.per[j].a[i]) * nvu), 1:2), Float64(1.01 * (coef + 1) * d.mu^(-2n)), Float64(coef))
                end
                maxcoef = max(maxcoef, coef)
                wb = max(wb, sl_err(-s.L, asym) / (sl_BL(epsT, p1, p2) + tail))
            end
        end
    end
    return (; wa, wb, nres, helper, maxcoef)
end

function ac6(T, nid)
    seed = vseed(6)
    epsT = sl_eps(T)
    rot = LinearSymplecticMap(SMatrix{2, 2, T}(0, -1, 1, 0))
    Mrot = SMatrix{2, 2, Float64}(0, -1, 1, 0)
    W0 = [sl_frame(seed, 1, T)]
    S0 = sl_sali0(W0[1])
    ns = (12, 120, 1200, 9996)
    Ms = fill(Mrot, 9999)
    ref = sl_ref(sl_mats_const(0, 1, -1, 0), W0[1], 9999; keep = collect(ns))
    wrot = 0.0
    refangle = 0.0
    for tau in (1, 3)
        run = sl_run(SALI(), rot, [SVector{2, T}(0, 0)], W0, tau == 1 ? 10_000 : 9_999, tau; steps = collect(ns))
        for s in run[1].snaps
            n = s.n
            at = sl_at(ref, n)
            p1, p2 = lb_Psi(view(Ms, 1:n), W0[1][:, 1]), lb_Psi(view(Ms, 1:n), W0[1][:, 2])
            wrot = max(wrot, sl_err(s.sali, S0) / sl_BSALI(epsT, sl_Lambda(ref.lnsali, n), Float64(S0), p1, p2))
            refangle = max(refangle, Float64(abs(at.S - S0)))
        end
    end
    sys = StandardMap(exact(T, -2))
    x0 = [SVector{2, T}(0, 0) for _ in 1:nid]
    W0 = [sl_frame(seed, j, T) for j in 1:nid]
    nsamp = collect(1:1000)
    Ms = lb_tangents(sys, x0[1], 1000)
    refs = [sl_ref(sl_mats_const(-1, 1, -2, 1), W0[j], 1000; keep = nsamp) for j in 1:nid]
    psis = [sl_psis(Ms, W0[j], nsamp) for j in 1:nid]
    S0s = [sl_sali0(W0[j]) for j in 1:nid]
    wref = weven = 0.0
    min_gap = Inf
    min_sali = Inf
    for tau in (1, 2)
        run = sl_run(SALI(), sys, x0, W0, 1000, tau; steps = filter(n -> n % tau == 0, nsamp))
        for j in 1:nid, s in run[j].snaps
            n = s.n
            at = sl_at(refs[j], n)
            p1, p2 = psis[j][n]
            saliref = Float64(at.S)
            B = sl_BSALI(epsT, sl_Lambda(refs[j].lnsali, n), saliref, p1, p2)
            wref = max(wref, sl_err(s.sali, at.S) / B)
            iseven(n) && (weven = max(weven, sl_err(s.sali, S0s[j]) / B))
            min_gap = min(min_gap, Float64(s.sali) - (0.146 - B))
            min_sali = min(min_sali, Float64(s.sali))
        end
    end
    return (; wrot, refangle, wref, weven, min_gap, min_sali)
end

section("AC4-6") do
    for T in SL_TYPES
        r = ac4(T, 100)
        emit("AC4", "ids with |c1(1) c1(2)| < 0.01 excluded (of 100), $T", r.excluded, "<= 10", pf(r.excluded <= 10))
        emit("AC4", "BigFloat product against the full eigen-expansion of A11 at n = 8, max |L difference|, $T", r.helper, "< 1e-60", pf(r.helper < 1e-60))
        ratio_line("AC4", "(a) max |SALI - SALI_cf| / B_SALI at resolved n in {1,2,4,8} (tau = 5: {5,10}), $(r.nres) points, 100 ids, tau in {1,5}, $T", r.wa)
        ratio_line("AC4", "(b) max |L - L_cf| / B_L, n in {20,100,200}, $T", r.wb)
        ratio_line("AC4", "(c) max |(-L) - (2 n ln mu - ln|D| + ln|c1 c1|)| / (B_L + mu^(-4n)), n in {20,100}, $T", r.wc)
    end
    for T in SL_TYPES, K in (1, 2, 6)
        r = ac5(T, K, 100)
        emit("AC5", "BigFloat product against the eigen-expansion of A12 at n = 8, K = $K, max |L difference|, $T", r.helper, "< 1e-60", pf(r.helper < 1e-60))
        ratio_line("AC5", "(a) K = $K: max |SALI - SALI_closed| / B_SALI at resolved n ($(r.nres) points), 100 ids, tau in {1,5}, $T", r.wa)
        ratio_line("AC5", "(b) K = $K: max |(-L) - [2 n ln mu+ - ln|det W0| + sum ln(|a_i| ||vu||)]| / (B_L + 1.01 (c + 1) mu+^(-2n)), n in {20,100,200}, $T (largest frame coefficient c = $(fmt(r.maxcoef)))", r.wb)
    end
    for T in SL_TYPES
        r = ac6(T, 100)
        ratio_line("AC6", "(a) order-4 rotation: max |SALI(n) - SALI(0)| / B_SALI, n in {12,120,1200,9996}, tau in {1,3}, $T", r.wrot)
        info("AC6", "(a) the BigFloat reference itself keeps the angle: max |SALI_ref - SALI(0)|, $T", r.refangle)
        ratio_line("AC6", "(b) K = -2 at (0,0): max |SALI - SALI_ref| / B_SALI over all n in 1..1000, tau in {1,2}, 100 ids (odd n: the closed form, even n: SALI(0) through the reference), $T", r.wref)
        ratio_line("AC6", "(b) even n: max |SALI(n) - SALI(0)| / B_SALI, $T", r.weven)
        emit("AC6", "(b) min over n and ids of SALI(n) - (0.146 - B_SALI); smallest SALI observed $(fmt(r.min_sali)) (A13 lower bound 0.1459), $T", r.min_gap, ">= 0", pf(r.min_gap >= 0))
    end
end

# ---------------------------------------------------------------------------------------------------
# AC7 and AC8: same-orbit reference, log form
section("AC7-8") do
    for T in SL_TYPES, (label, N, m) in sl_panel_v(T)
        meas = sl_meas(T)
        K = label === :half ? "0.5 regular" : label === :two_pi ? "2 pi" : "10"
        rows = vrows(label, T, 1, m, N)
        emit("AC7", "K = $K, $T: rows (starts x checkpoints at every multiple of 10 up to $N)", length(rows), "$(m * (N ÷ 10))", pf(length(rows) == m * (N ÷ 10)))
        ratio_line("AC7", "K = $K, $T: max |SALI_T - SALI_ref| / B_SALI (realized ratio R4), $m starts, n up to $N", maximum(r.errS / r.BS for r in rows))
        info("AC7", "K = $K, $T: max |error| in units of B_meas / median ratio error/B_SALI", string(fmt(maximum(r.errS for r in rows) / meas), " / ", fmt(sl_median([r.errS / r.BS for r in rows]))))
        if label !== :half
            below = [r for r in rows if r.saliref < meas]
            emit("AC7", "K = $K, $T: floor control, handoffs with SALI_ref < B_meas ($(length(below))) where SALI_T > 2 B_meas; resolved handoffs exist", string(count(r -> r.sali > 2 * meas, below), ", ", any(r -> r.saliref >= 10 * meas, rows)), "0, true",
                 pf(!isempty(below) && all(r -> r.sali <= 2 * meas, below) && any(r -> r.saliref >= 10 * meas, rows)))
        end
        for tau in (1, 4)
            rows = vrows(label, T, tau, m, N)
            epsT = sl_eps(T)
            ratio_line("AC8", "(a) K = $K, $T, tau = $tau: max |L - L_ref| / B_L (realized ratio R4), $(length(rows)) points", maximum(r.errL / r.BL for r in rows))
            res = [r for r in rows if r.saliref >= 10 * meas]
            if !isempty(res)
                ratio_line("AC8", "(c) K = $K, $T, tau = $tau: max |SALI_direct - sali_from_log(L)| / (B_SALI + 8 eps (Psi1 + Psi2) SALI_ref) over the $(length(res)) resolved points",
                           maximum(abs(Float64(r.sali - sali_from_log(r.L))) / (r.BS + 8 * epsT * (r.psi1 + r.psi2) * r.saliref) for r in res))
            else
                emit("AC8", "(c) K = $K, $T, tau = $tau: no resolved point at the handoffs (chaotic, first handoff below the floor)", 0, "-", "INFO")
            end
            if label !== :half
                late = [r for r in rows if r.n >= 50]
                dmin = minimum(abs(Float64(-r.S1) - r.Lref) / (10 * r.BL) for r in late)
                emit("AC8", "(d) K = $K, $T, tau = $tau: negative control S_1 only, min |(-S_1) - L_ref| / (10 B_L) over n >= 50 ($(length(late)) points)", dmin, "> 1", pf(dmin > 1))
            end
        end
    end
    # (b) L finite at every handoff for K = 2 pi, n = 1e3, Float64, where exp(L) underflows
    c = vcase(:two_pi, Float64, 20)
    run = sl_run(SALI(), c.sys, c.x0, c.W0, 1_000, 1; all = true)
    Ls = [s.L for r in run for s in r.snaps]
    emit("AC8", "(b) K = 2 pi, Float64, 20 starts x 1000 handoffs: all L finite; min L (exp(min L) = $(exp(minimum(Ls))))", minimum(Ls), "finite, < -745, exp = 0", pf(length(Ls) == 20 * 1000 && all(isfinite, Ls) && minimum(Ls) < -745 && exp(minimum(Ls)) == 0))
end

# ---------------------------------------------------------------------------------------------------
# AC9: cross-check with L1
section("AC9") do
    T = Float64
    ns = [100, 1_000, 10_000]
    for label in (:half, :two_pi, :ten)
        m = 20
        c = vcase(label, T, m)
        run = sl_run(SALI(), c.sys, c.x0, c.W0, 10_000, 1; steps = ns)
        l1a = sl_run(Benettin(), c.sys, c.x0, [SMatrix{2, 1, T}(W[:, 1]) for W in c.W0], 10_000, 1; steps = ns)
        l1b = sl_run(Benettin(), c.sys, c.x0, [SMatrix{2, 1, T}(W[:, 2]) for W in c.W0], 10_000, 1; steps = ns)
        wa = wb = wc = 0.0
        gaps = Float64[]
        for j in 1:m
            Ms = lb_tangents(c.sys, c.x0[j], 10_000)
            for (i, n) in enumerate(ns)
                s = run[j].snaps[i]
                nX1 = l1a[j].snaps[i].sum[1] + l1a[j].snaps[i].comp[1]
                nX2 = l1b[j].snaps[i].sum[1] + l1b[j].snaps[i].comp[1]
                p1, p2 = lb_Psi(view(Ms, 1:n), c.W0[j][:, 1]), lb_Psi(view(Ms, 1:n), c.W0[j][:, 2])
                wa = max(wa, abs((s.sum[1] + s.comp[1]) - nX1) / (2 * 8 * EPS64 * p1))
                wb = max(wb, abs((s.sum[2] + s.comp[2]) - nX2) / (2 * 8 * EPS64 * p2))
                wc = max(wc, abs(s.L + nX1 + nX2) / (2 * sl_BL(EPS64, p1, p2)))
                n == 10_000 && push!(gaps, abs(nX2 - nX1) / 2)
            end
        end
        K = label === :half ? "0.5 regular" : label === :two_pi ? "2 pi" : "10"
        ratio_line("AC9", "(a) K = $K: max |S_1 - n X_L1(first column)| / (2 B_S^Psi,1), 20 starts, n in {100,1e3,1e4}", wa)
        ratio_line("AC9", "(b) K = $K: max |S_2 - n X_L1(second column)| / (2 B_S^Psi,2)", wb)
        ratio_line("AC9", "(c) K = $K: max |L + n (X1 + X2)| / (2 B_L)", wc)
        info("AC9", "K = $K: |ln(norm w2/norm w1)|/2 at n = 1e4 (A8, of order 1): median; max", string(fmt(sl_median(gaps)), "; ", fmt(maximum(gaps))))
    end
end

# ---------------------------------------------------------------------------------------------------
# AC10 and AC11: chaotic rate and regular decay law (ensemble level)
section("AC10-11") do
    m = 1000
    for (ci, label) in enumerate((:two_pi, :five, :ten))
        seed = vseed(10, ci)
        sys = StandardMap(sl_K(Float64, label))
        x0 = [sl_cast(Float64, s) for s in sl_starts(StableRNG(seed), m)]
        W0 = [sl_frame(seed, j, Float64) for j in 1:m]
        run = sl_run(SALI(), sys, x0, W0, 200, 1; steps = 20:200)
        l1 = sl_run(Benettin(), sys, x0, [SMatrix{2, 1, Float64}(W[:, 1]) for W in W0], 10_000, 1; steps = [10_000])
        hs = [sl_slope(Float64[s.n for s in r.snaps], Float64[-s.L for s in r.snaps]) / 2 for r in run]
        X = [(r.snaps[1].sum[1] + r.snaps[1].comp[1]) / 10_000 for r in l1]
        se = sqrt(var(hs) / m + var(X) / m)
        d = abs(mean(hs) - mean(X))
        K = label === :two_pi ? "2 pi" : string(sl_K(Float64, label))
        ratio_line("AC10", "K = $K, 1000 starts: |mean(m/2) - mean(X(1e4))| / (4 combined SE + 2/180); mean m/2 = $(fmt(mean(hs))), mean X = $(fmt(mean(X))), sd of m/2 = $(fmt(std(hs)))", d / (4 * se + 2 / 180))
        info("AC10", "K = $K: realized |difference| of the means; combined SE; |difference| / SE (finite-window allowance of the spec: $(fmt(2 / 180)))", string(fmt(d), "; ", fmt(se), "; ", fmt(d / se)))
        info("AC10", "K = $K: Chirikov sanity rate 2 ln(K/2) (loose, K18) against mean(m)", string(fmt(2 * log(sl_K(Float64, label) / 2)), " vs ", fmt(2 * mean(hs))))
    end
    # AC11
    mr = 20
    c = vcase(:half, Float64, mr)
    ns = unique(round.(Int, 10 .^ range(4, 6; length = 21)))
    run = sl_run(SALI(), c.sys, c.x0, c.W0, 1_000_000, 1; steps = ns)
    slopes = Float64[sl_slope(log.(Float64.(ns)), Float64[s.L for s in r.snaps]) for r in run]
    transient = count(r -> r.snaps[1].sali > 1e-2, run)
    passed = count(i -> run[i].snaps[1].sali <= 1e-2 && -2.4 <= slopes[i] <= -1.6, 1:mr)
    emit("AC11", "starts with a LS slope of L against ln n over [1e4,1e6] in [-2.4,-1.6] (of $mr; slopes: min $(fmt(minimum(slopes))), median $(fmt(sl_median(slopes))), max $(fmt(maximum(slopes))))", passed, ">= $(ceil(Int, 0.8 * mr))", pf(passed >= ceil(Int, 0.8 * mr)))
    emit("AC11", "starts in the transient at n = 1e4 (SALI > 1e-2)", transient, "<= $(floor(Int, 0.2 * mr))", pf(transient <= 0.2 * mr))
    n2a = sort([Float64(r.snaps[1].sali * r.snaps[1].n^2) for r in run])
    n2b = sort([Float64(r.snaps[end].sali * r.snaps[end].n^2) for r in run])
    info("AC11", "R2: n^2 SALI at n = 1e4 (min; median; max)", string(fmt(n2a[1]), "; ", fmt(sl_median(n2a)), "; ", fmt(n2a[end])))
    info("AC11", "R2: n^2 SALI at n = 1e6 (min; median; max)", string(fmt(n2b[1]), "; ", fmt(sl_median(n2b)), "; ", fmt(n2b[end])))
    DATA[:slopes] = slopes
end

# ---------------------------------------------------------------------------------------------------
# AC12: precision ladder
function v_ens(::Type{T}, label::Symbol, m::Int, n::Int) where {T}
    c = vcase(label, T, m)
    return map(1:m) do j
        _, _, st = propagate_orbit(c.sys, c.x0[j], c.W0[j], n, 1, (SALI(),); orbit = j)
        -Float64(reducer_final(SALI(), st[1]).L) / (2n)
    end
end

section("AC12") do
    # (a) floor scaling and R1
    nres = Dict(:five => (18, 36, 94), :two_pi => (14, 28, 74), :ten => (10, 20, 54))
    for label in (:five, :two_pi, :ten)
        nfl = Dict{DataType, Vector{Int}}()
        for T in SL_TYPES
            c = vcase(label, T, 100)
            run = sl_run(SALI(), c.sys, c.x0, c.W0, 150, 1; all = true)
            nfl[T] = map(run) do r
                i = findfirst(s -> s.sali < sl_meas(T), r.snaps)
                i === nothing ? 0 : r.snaps[i].n
            end
        end
        both = findall(j -> nfl[Float64][j] > 0 && nfl[Double64][j] > 0, 1:100)     # starts that reach the floor in both types within 150 steps (the others stay regular)
        med = sl_median(nfl[Double64][both] ./ nfl[Float64][both])
        K = label === :two_pi ? "2 pi" : string(sl_K(Float64, label))
        if label !== :ten
            emit("AC12", "(a) K = $K, 100 starts, n = 150: median n_floor(Double64) / n_floor(Float64) over the $(length(both)) starts that reach the floor in both types (never reaching it within 150 steps: Float64 $(count(==(0), nfl[Float64])), Double64 $(count(==(0), nfl[Double64])), regular islands)", med, "[1.5, 2.4], >= 80 starts", pf(length(both) >= 80 && 1.5 <= med <= 2.4))
        end
        nf, nd = nfl[Float64][both], nfl[Double64][both]
        info("AC12", "R1 K = $K: n_floor over the $(length(both)) starts reaching the floor: Float64 (min; median; max) | Double64 (min; median; max) | A14 loose n_res Float64, Double64, BigFloat",
             string(minimum(nf), "; ", sl_median(nf), "; ", maximum(nf), " | ", minimum(nd), "; ", sl_median(nd), "; ", maximum(nd), " | ", nres[label]))
    end
    # (b) BigFloat(256) run against the 512-bit reference
    setprecision(BigFloat, SL_BITS) do
        m = 50
        ns = [10, 20, 40, 60]
        c = vcase(:two_pi, BigFloat, m)
        epsT = sl_eps(BigFloat)
        run = sl_run(SALI(), c.sys, c.x0, c.W0, 60, 1; steps = ns)
        w = 0.0
        nr = 0
        for j in 1:m
            pts, _ = lb_pts(c.sys, c.x0[j], 60)
            ref = sl_ref(sl_mats_std(c.K, pts), c.W0[j], 60; keep = ns, bits = 512)
            Ms = lb_tangents(c.sys, c.x0[j], 60)
            for s in run[j].snaps
                at = sl_at(ref, s.n)
                p1, p2 = lb_Psi(view(Ms, 1:s.n), c.W0[j][:, 1]), lb_Psi(view(Ms, 1:s.n), c.W0[j][:, 2])
                saliref = Float64(at.S)
                nr += saliref >= 10 * 4 * epsT
                w = max(w, sl_err(s.sali, at.S; bits = 512) / sl_BSALI(epsT, sl_Lambda(ref.lnsali, s.n), saliref, p1, p2))
            end
        end
        ratio_line("AC12", "(b) BigFloat(256), K = 2 pi, 50 starts, n in {10,20,40,60}: max |SALI - SALI_ref(512 bits)| / B_SALI (own eps $(fmt(epsT))); $nr points resolved", w)
        emit("AC12", "(b) resolved points exist", nr, "> 0", pf(nr > 0))
    end
    # (c) ensembles
    for label in (:two_pi, :five)
        a, b = v_ens(Float64, label, 1000, 1_000), v_ens(Double64, label, 1000, 1_000)
        K = label === :two_pi ? "2 pi" : "5"
        ratio_line("AC12", "(c) K = $K, 1000 starts, n = 1e3: |mean(-L/(2n)) Float64 - Double64| / (4 combined SE); means $(fmt(mean(a))) and $(fmt(mean(b))), z = $(fmt(zdiff(a, b)))", zdiff(a, b) / 4)
    end
    setprecision(BigFloat, SL_BITS) do
        for label in (:two_pi, :five)
            a, b = v_ens(BigFloat, label, 100, 300), v_ens(Double64, label, 100, 300)
            K = label === :two_pi ? "2 pi" : "5"
            ratio_line("AC12", "(c) K = $K, 100 starts, n = 300: BigFloat(256) vs Double64 |difference of means| / (4 combined SE); means $(fmt(mean(a))) and $(fmt(mean(b))), z = $(fmt(zdiff(a, b)))", zdiff(a, b) / 4)
        end
    end
end

# ---------------------------------------------------------------------------------------------------
# AC13: health and rejection
section("AC13") do
    seed = vseed(13, 1)
    x0s(flagged) = [SVector{3, Float64}(0, 0, j in flagged ? 1 : 0) for j in 1:6]
    zero_map, col2_map = sl_Fault(:zero, 0), sl_Fault(:col2, 0)
    ok = propagate(zero_map, x0s(()), 12; k = 2, tau = 3, seed = seed, reducers = (SALI(),), chunk = 2)
    res = nothing
    err = sl_health_error(() -> (res = propagate(zero_map, x0s((4,)), 12; k = 2, tau = 3, seed = seed, reducers = (SALI(),), chunk = 2)))
    res2 = nothing
    err2 = sl_health_error(() -> (res2 = propagate(zero_map, x0s((5, 2)), 12; k = 2, tau = 3, seed = seed, reducers = (SALI(),), chunk = 2)))
    eye = [SMatrix{2, 2, Float64}(1, 0, 0, 1) for _ in 1:6]
    err3 = sl_health_error(() -> propagate(col2_map, x0s((5,)), 12; k = 2, tau = 3, seed = seed, reducers = (SALI(),), chunk = 2, W0 = eye))
    x = SVector{3, Double64}(0, 0, 1)
    err4 = sl_health_error(() -> propagate_orbit(zero_map, x, sl_frame(seed, 1, Double64), 12, 3, (SALI(),); orbit = 4))
    err5 = sl_health_error(() -> propagate_orbit(col2_map, x, SMatrix{2, 2, Double64}(1, 0, 0, 1), 12, 3, (SALI(),); orbit = 4))
    e1 = try
        propagate(sl_Boom(), [SVector(0.0, 0.0)], 4; k = 1, tau = 2, seed = seed, reducers = (SALI(),)); nothing
    catch ex
        ex
    end
    oka = length(ok.finals) == 6 && err isa NumericalHealthError && res === nothing && err.what == "deviation matrix" && err.orbit == 4 && err.step == 3 &&
          err.backend === :cpu && err.column == 1 && err.value == 0.0 && err2 isa NumericalHealthError && err2.orbit == 2 && res2 === nothing &&
          err3 isa NumericalHealthError && err3.orbit == 5 && err3.column == 2 && err3.step == 3 && err4 isa NumericalHealthError && err4.orbit == 4 && err4.column == 1 &&
          err5 isa NumericalHealthError && err5.column == 2 && e1 isa ArgumentError && occursin("SALI needs k = 2, got k = 1", e1.msg)
    emit("AC13", "(a) zero tangent on a flagged orbit (propagate chunk 2 and Double64 per-orbit): kind deviation matrix, orbit, step, backend cpu, column 1 (2 for diag(1,0)), value 0, nothing stored, lowest id first; k = 1 ArgumentError", oka, "true", pf(oka))
    seed = vseed(13, 2)
    resn = nothing
    errn = sl_health_error(() -> (resn = propagate(sl_Fault(:nan, 3), x0s((3,)), 20; k = 2, tau = 5, seed = seed, reducers = (SALI(),), chunk = 2)))
    errn2 = sl_health_error(() -> propagate_orbit(sl_Fault(:nan, 3), x, sl_frame(seed, 1, Double64), 20, 5, (SALI(),); orbit = 3))
    okb = errn isa NumericalHealthError && resn === nothing && errn.what == "deviation matrix" && errn.orbit == 3 && errn.step == 5 && errn.backend === :cpu &&
          errn.column == 0 && errn.threshold == Inf && !isfinite(errn.value) && errn2 isa NumericalHealthError && errn2.column == 0 && errn2.step == 5
    emit("AC13", "(b) NaN in the tangent: kind deviation matrix, column 0, threshold Inf, orbit 3, step 5 (not masked), propagate and Double64 per-orbit", okb, "true", pf(okb))
    seed = vseed(13, 3)
    r = SALI()
    for T in SL_TYPES
        meas = sl_meas(T)
        frame(S) = (phi = 2 * asin(T(S) / 2); SMatrix{2, 2, T}(1, 0, cos(phi), sin(phi)))
        dep = [sl_health_error(() -> check_frame(r, W; orbit = 7)) for W in (SMatrix{2, 2, T}(1, 0, 1, 0), SMatrix{2, 2, T}(1, 0, -1, 0), SMatrix{2, 2, T}(2, 3, 2, 3))]
        W5, W20 = frame(5 * meas), frame(20 * meas)
        S5, S20 = sl_sali_exact(W5[:, 1], W5[:, 2])[3], sl_sali_exact(W20[:, 1], W20[:, 2])[3]
        e5 = sl_health_error(() -> check_frame(r, W5; orbit = 3))
        s0 = nothing
        e20 = sl_health_error(() -> (s0 = check_frame(r, W20; orbit = 3)))
        nid = T === Float64 ? 10_000 : 1_000
        minS0 = minimum(j -> Float64(check_frame(r, sl_frame(seed, j, T); orbit = j)), 1:nid)
        okc = all(e -> e isa NumericalHealthError && e.what == "initial frame" && e.orbit == 7 && e.step == 0, dep) && 4 * meas < S5 < 6 * meas && 19 * meas < S20 < 21 * meas &&
              e5 isa NumericalHealthError && e5.what == "initial frame" && e20 === nothing && sl_err(s0, S20) <= meas
        emit("AC13", "(c) $T: dependent frames and SALI(0) = $(fmt(Float64(S5) / meas)) B_meas rejected (kind initial frame, step 0, orbit); SALI(0) = $(fmt(Float64(S20) / meas)) B_meas accepted and returned", okc, "true", pf(okc))
        emit("AC13", "(c) $T: seeded frame accepted for all ids 1..$nid; min SALI(0) / (10 B_meas)", minS0 / (10 * meas), ">= 1", pf(minS0 >= 10 * meas))
    end
    W0 = [sl_frame(vseed(13, 4), 1, Float64), SMatrix{2, 2, Float64}(1, 0, 1, 0)]
    erb = sl_health_error(() -> begin
        foreach(j -> check_frame(SALI(), W0[j]; orbit = j), eachindex(W0))
        propagate_orbit(sl_Boom(), SVector(0.0, 0.0), W0[2], 4, 2, (SALI(),))
    end)
    boom = try
        propagate_orbit(sl_Boom(), SVector(0.0, 0.0), W0[1], 4, 2, (SALI(),)); false
    catch ex
        ex isa ErrorException
    end
    okr = erb isa NumericalHealthError && erb.what == "initial frame" && erb.orbit == 2 && erb.step == 0 && boom
    emit("AC13", "(c) caller-side frame check rejects the dependent frame before any map call (sl_Boom throws if called); an accepted frame reaches the map (Amendment D67 item 1)", okr, "true", pf(okr))
    mu = exp(3.0)
    sys = LinearSymplecticMap(SMatrix{2, 2, Float64}(mu, 0, 0, 1 / mu))
    Wd = sl_frame(vseed(13, 5), 1, Float64)
    e200 = sl_health_error(() -> propagate_orbit(sys, SVector(0.0, 0.0), Wd, 200, 200, (SALI(),)))
    e300 = sl_health_error(() -> propagate_orbit(sys, SVector(0.0, 0.0), Wd, 300, 300, (SALI(),)))
    okd = e200 === nothing && e300 isa NumericalHealthError && e300.what == "deviation matrix" && e300.step == 300 && e300.backend === :cpu
    emit("AC13", "(d) A = diag(e^3, e^-3), Float64: tau = 200 no error at n = 200; tau = 300 error at handoff 300", okd, "true", pf(okd))
    res = propagate(StandardMap(2.0), [SVector(1.0, 2.0)], 8; k = 2, tau = 2, seed = vseed(13, 6), reducers = (SALI(),))
    f = res.finals[1][1]
    oke = reducer_provenance(SALI()).indicator == "sali" && res.provenance.indicator == "sali" && res.provenance.k == 2 && res.provenance.tau == 2 && res.provenance.precision == "Float64" &&
          f.n == 8 && all(T -> begin
              _, _, st = propagate_orbit(StandardMap(exact(T, 2)), sl_cast(T, SVector(1.0, 2.0)), sl_frame(vseed(13, 6), 1, T), 8, 2, (SALI(),))
              reducer_final(SALI(), st[1]).floor == 4 * effective_eps(T) && sali_floor(T) == 4 * effective_eps(T)
          end, SL_TYPES)
    emit("AC13", "(e) run record: indicator, k, tau, precision; floor B_meas(T) = 4 eps in the final (Float64 $(fmt(4 * EPS64)), Double64 $(fmt(4 * EPSDD)))", oke, "true", pf(oke))
end

# ---------------------------------------------------------------------------------------------------
# AC14: orbit and execution invariance
section("AC14") do
    seed = vseed(14)
    T = Float64
    bad = 0
    ncmp = 0
    for tau in (1, 4)
        sys = StandardMap(sl_K(T, :two_pi))
        for (j, s) in enumerate(sl_starts(StableRNG(seed + 1), 5))
            x0 = sl_cast(T, s)
            W2 = sl_frame(seed, j, T)
            W1 = SMatrix{2, 1, T}(W2[:, 1])
            xa, _, sta = propagate_orbit(sys, x0, W1, 400, tau, (sl_Trace(),); orbit = j)
            xb, _, stb = propagate_orbit(sys, x0, W1, 400, tau, (Benettin(), sl_Trace()); orbit = j)
            xc, _, stc = propagate_orbit(sys, x0, W2, 400, tau, (SALI(), sl_Trace()); orbit = j)
            xd, _, std = propagate_orbit(sys, x0, W2, 400, tau, (sl_Trace(),); orbit = j)
            tra, trb = reducer_final(sl_Trace(), sta[1]).xs, reducer_final(sl_Trace(), stb[2]).xs
            trc, trd = reducer_final(sl_Trace(), stc[2]).xs, reducer_final(sl_Trace(), std[1]).xs
            bad += !(isequal(tra, trb) && isequal(tra, trc) && isequal(tra, trd) && xa === xb && xa === xc && xa === xd)
            ncmp += 1
        end
    end
    emit("AC14", "(a) orbit points at every handoff differ across the chains {}, {L1} (k=1), {SALI} (k=2), {} (k=2); cases = (5 starts x tau in {1,4})", bad, "0 of $ncmp", pf(bad == 0))
    for tau in (1, 4)
        sys = StandardMap(sl_K(T, :two_pi))
        reducers = (SALI(),)
        m, N = 4096, 1000
        x0 = [sl_cast(T, s) for s in sl_starts(StableRNG(seed + 2), m)]
        W0 = [sl_frame(seed, j, T) for j in 1:m]
        serial = [propagate_orbit(sys, x0[j], W0[j], N, tau, reducers; orbit = j) for j in 1:m]
        xs, Ws = [s[1] for s in serial], [s[2] for s in serial]
        fs = [map(reducer_final, reducers, s[3]) for s in serial]
        for chunk in (64, 1000)
            res = propagate(sys, x0, N; k = 2, tau = tau, seed = seed, reducers = reducers, chunk = chunk, W0 = W0)
            okb = isequal(res.x, xs) && isequal(res.W, Ws) && isequal(res.finals, fs) && res.provenance.indicator == "sali"
            emit("AC14", "(b) propagate equals the serial map of propagate_orbit bitwise (x, W, finals), 4096 orbits, K = 2 pi, n = 1e3, tau = $tau, chunk = $chunk", okb, "true", pf(okb))
        end
    end
    hashes = Dict{Int, String}(Threads.nthreads() => string(sl_thread_hash()))
    for nt in (1, 4)
        nt == Threads.nthreads() && continue
        cmd = addenv(`$(Base.julia_cmd()) -t $nt --project=$(dirname(Base.active_project())) $(@__FILE__)`, "SL_CHILD" => "1")
        txt = read(pipeline(cmd; stderr = devnull), String)
        mt = match(r"HASH (\d+) THREADS (\d+)", txt)
        hashes[nt] = mt === nothing ? "child failed: $(first(txt, 100))" : mt.captures[1]
    end
    okt = length(unique(values(hashes))) == 1
    emit("AC14", "(b) x, W, SALI, L, sums and compensations bitwise identical for thread counts $(sort(collect(keys(hashes)))) (hashes of all outputs, tau in {1,4}, separate processes)", okt, "true", pf(okt))
    for T in SL_TYPES, tau in (1, 4)
        a = sl_alloc_bytes(StandardMap(sl_K(T, :two_pi)), sl_cast(T, SVector(1.0, 2.0)), sl_frame(seed, 1, T), 10_000, tau, (SALI(),))
        emit("AC14", "(c) bytes allocated by propagate_orbit after warm-up, n = 1e4, tau = $tau, $T", a, "0", pf(a == 0))
    end
    m = 1000
    c = vcase(:two_pi, Float64, m)
    run = sl_run(SALI(), c.sys, c.x0, c.W0, 1_000, 1; all = true)
    reached = violations = 0
    for r in run
        sal = [s.sali for s in r.snaps]
        i = findfirst(==(0), sal)
        i === nothing && continue
        reached += 1
        violations += all(==(0), sal[i:end]) ? 0 : 1
    end
    emit("AC14", "(d) absorbing zero, K = 2 pi, Float64, 1000 starts, n = 1e3: starts that reach exactly 0 ($reached); starts with a nonzero value afterwards", violations, "0", pf(violations == 0 && reached >= m ÷ 2))
end

# ---------------------------------------------------------------------------------------------------
# AC15: wrapped versus lifted
section("AC15") do
    seed = vseed(15)
    rng = StableRNG(seed)
    dy() = rand(rng, 0:2048) / 1024
    T = Float64
    nbad = 0
    for tau in (1, 3)
        zs = [SVector{2, T}(T(dy()), T(dy())) for _ in 1:50]
        W0 = [sl_frame(seed, j, T) for j in 1:50]
        rl = sl_run(SALI(), StandardMap(exact(T, 0)), zs, W0, 300, tau; all = true)
        rw = sl_run(SALI(), WrappedStandardMap(exact(T, 0)), torus_state.(zs), W0, 300, tau; all = true)
        nbad += count(j -> !isequal(rl[j].snaps, rw[j].snaps), 1:50)
    end
    emit("AC15", "(a) K = 0, dyadic starts: starts (of 2 x 50, tau in {1,3}, all handoffs to 300) whose SALI, L, sums, compensations differ between wrapped and lifted", nbad, "0", pf(nbad == 0))
    m = 50
    c = vcase(:half, T, m)
    for n in (10, 100)
        rl = sl_run(SALI(), c.sys, c.x0, c.W0, n, 1; steps = [n])
        rw = sl_run(SALI(), WrappedStandardMap(0.5), torus_state.(c.x0), c.W0, n, 1; steps = [n])
        excluded = 0
        w = 0.0
        diffmax = 0.0
        for j in 1:m
            _, rho1, M = lb_ref_full(0.5, c.x0[j], c.W0[j][:, 1], n, 0.0, 0.0)
            _, rho2, _ = lb_ref_full(0.5, c.x0[j], c.W0[j][:, 2], n, 0.0, 0.0)
            if max(rho1, rho2) > 50
                excluded += 1
                continue
            end
            d = abs(rw[j].snaps[1].L - rl[j].snaps[1].L)
            diffmax = max(diffmax, d)
            w = max(w, d / (64 * 0.5 * M * n^2 * EPS64 * (rho1 + rho2)))
        end
        ratio_line("AC15", "(b) K = 0.5 regular class, n = $n, 50 starts: max |L_wrapped - L_lifted| / (64 K M_state n^2 eps (rho1 + rho2)); max |difference| $(fmt(diffmax)); excluded (rho > 50) $excluded", w)
        emit("AC15", "(b) n = $n: excluded starts (of 50)", excluded, "<= 5", pf(excluded <= 5))
    end
end

# ---------------------------------------------------------------------------------------------------
# AC16: ChaosTools GALI_2 oracle (P8 confirmations first)
ct_rule(u, p, n) = (pp = u[2] + p[1] * sin(u[1]); SVector(u[1] + pp, pp))
ct_jac(u, p, n) = (c = p[1] * cos(u[1]); SMatrix{2, 2}(1 + c, c, 1, 1))
const CT_TDS = Dict{Float64, Any}()
ct_tds(K) = get!(CT_TDS, K) do
    CT.TangentDynamicalSystem(CT.DeterministicIteratedMap(ct_rule, SVector(0.1, 0.1), [K]); J = ct_jac, k = 2)
end
function ct_gali(K, z0, W0, n; threshold = 0, Q0 = W0)
    tds = ct_tds(K)
    CT.reinit!(tds, SVector{2, Float64}(z0); Q0 = Q0)
    return CT.gali(tds, n; threshold = threshold)[1]
end
# Number of steps (of n) after which the oracle orbit differs bitwise from forward(StandardMap), and the largest coordinate difference.
function ct_orbit_check(K, z0, n)
    tds = ct_tds(K)
    CT.reinit!(tds, SVector{2, Float64}(z0); Q0 = SMatrix{2, 2, Float64}(1, 0, 0, 1))
    sys = StandardMap(K)
    z = SVector{2, Float64}(z0)
    mism, maxd = 0, 0.0
    for _ in 1:n
        CT.step!(tds)
        z = forward(sys, z)
        zo = SVector{2, Float64}(CT.current_state(tds))
        mism += zo == z ? 0 : 1
        maxd = max(maxd, maximum(abs, zo - z))
    end
    return mism, maxd
end
g_from_sali(s) = s * sqrt(4 - s^2) / 2        # A4: GALI_2 = SALI max(d-, d+) / 2, max(d-, d+) = sqrt(4 - SALI^2)

section("AC16") do
    meas = sl_meas(Float64)
    DATA[:ct] = NamedTuple[]
    # (P8-1) threshold: the default 1e-12 stops the run early, threshold 0 runs the full horizon
    c5 = vcase(:five, Float64, 1000)
    lens_def = [length(ct_gali(5.0, c5.x0[j], c5.W0[j], 100; threshold = 1e-12)) for j in 1:200]
    lens_0 = [length(ct_gali(5.0, c5.x0[j], c5.W0[j], 100; threshold = 0)) for j in 1:200]
    nstop = count(<(101), lens_def)
    ok1 = nstop >= 180 && all(==(101), lens_0)
    emit("AC16", "(P8-1) K = 5, 200 starts, 100 steps requested: gali with the default threshold 1e-12 stops early for $nstop of 200 starts (steps at stop: min; median of those; the others are regular islands that run the full 100); threshold = 0 gives the full length for all", string(minimum(lens_def) - 1, "; ", sl_median(filter(<(101), lens_def)) - 1, " | ", unique(lens_0) .- 1), ">= 180 stop early | 100", pf(ok1))
    # (P8-2) prescribed initial frame is used (Q0 of reinit!): the stored deviations equal W0, g(1) agrees with ours, and the default identity frame gives a different value
    n_eq = 0
    wfirst = 0.0
    ndiff = 0
    for j in 1:200
        tds = ct_tds(5.0)
        CT.reinit!(tds, SVector{2, Float64}(c5.x0[j]); Q0 = c5.W0[j])
        n_eq += CT.current_deviations(tds) == c5.W0[j] ? 1 : 0
        g = ct_gali(5.0, c5.x0[j], c5.W0[j], 1)
        gid = ct_gali(5.0, c5.x0[j], c5.W0[j], 1; Q0 = SMatrix{2, 2, Float64}(1, 0, 0, 1))
        _, _, st = propagate_orbit(c5.sys, c5.x0[j], c5.W0[j], 1, 1, (SALI(),); orbit = j)
        gc = exp(reducer_final(SALI(), st[1]).L)
        wfirst = max(wfirst, abs(g[2] - gc))
        ndiff += abs(g[2] - gid[2]) > 1e-3 ? 1 : 0
    end
    ok2 = n_eq == 200 && wfirst <= 4 * EPS64 && ndiff > 100
    emit("AC16", "(P8-2) prescribed frame: stored deviations == W0 after reinit! (of 200: $n_eq); max |g_oracle(1) - exp(L(1))| = $(fmt(wfirst)); starts where the identity frame gives g(1) differing by > 1e-3 from the prescribed frame: $ndiff", wfirst, "<= 4 eps, 200 stored, > 100 differ", pf(ok2))
    # (P8-3 = criterion 16(c)) the oracle orbit equals the stream's lifted orbit bitwise
    for (label, m, n) in ((:half, 20, 10_000), (:two_pi, 100, 12), (:five, 100, 12))
        starts = label === :half ? vcase(:half, Float64, 20).x0 : vcase(label, Float64, 1000).x0[1:m]
        K = Float64(sl_K(Float64, label))
        tot = 0
        mism = 0
        maxd = 0.0
        for z0 in starts
            mm, dd = ct_orbit_check(K, z0, n)
            mism += mm; maxd = max(maxd, dd); tot += n
        end
        emit("AC16", "(c) K = $K: oracle orbit vs forward(StandardMap) over $(length(starts)) starts x $n steps: steps with a bitwise difference; largest coordinate difference $(fmt(maxd))", mism, "0 of $tot (else add the orbit divergence term)", pf(mism == 0))
    end
    # (a) K = 0.5 regular class, 20 starts, n in {10, 100, 1e3, 1e4}
    rows = vrows(:half, Float64, 1, 20, 10_000)
    c = vcase(:half, Float64, 20)
    gs = [ct_gali(0.5, c.x0[j], c.W0[j], 10_000) for j in 1:20]
    wS = wL = 0.0
    maxabs = 0.0
    for r in rows
        r.n in (10, 100, 1000, 10_000) || continue
        go = gs[r.j][r.n + 1]
        gS, gL = g_from_sali(Float64(r.sali)), exp(Float64(r.L))
        wS = max(wS, abs(gS - go) / (2 * r.BS))
        wL = max(wL, abs(gL - go) / (2 * r.BS))
        maxabs = max(maxabs, abs(gS - go))
        push!(DATA[:ct], (K = 0.5, n = r.n, gC = gL, go = go, bound = 2 * r.BS, resolved = r.saliref >= 10 * meas))
    end
    ratio_line("AC16", "(a) K = 0.5 regular, 20 starts, n in {10,100,1e3,1e4}: max |g(SALI by A4) - g_oracle| / (B_SALI + B_oracle); max absolute difference $(fmt(maxabs))", wS)
    ratio_line("AC16", "(a) same with g = exp(L) (the floor-free form)", wL)
    # (b) K = 2 pi and 5, 1000 starts, n in {2, 5, 8, 12}
    for label in (:two_pi, :five)
        K = Float64(sl_K(Float64, label))
        c = vcase(label, Float64, 1000)
        rows = vrows(label, Float64, 1, 1000, 12; cps = [2, 5, 8, 12])
        tot = excl = 0
        wS = wL = 0.0
        for r in rows
            go = ct_gali(K, c.x0[r.j], c.W0[r.j], 12)[r.n + 1]
            tot += 1
            if !(r.saliref >= 10 * meas && go >= 10 * meas)
                excl += 1
                continue
            end
            gS, gL = g_from_sali(Float64(r.sali)), exp(Float64(r.L))
            wS = max(wS, abs(gS - go) / (2 * r.BS))
            wL = max(wL, abs(gL - go) / (2 * r.BS))
            push!(DATA[:ct], (K = K, n = r.n, gC = gL, go = go, bound = 2 * r.BS, resolved = true))
        end
        ratio_line("AC16", "(b) K = $(fmt(K)), 1000 starts, resolved handoffs n in {2,5,8,12}: max |g(SALI by A4) - g_oracle| / (B_SALI + B_oracle)", wS)
        ratio_line("AC16", "(b) K = $(fmt(K)): same with g = exp(L)", wL)
        emit("AC16", "(b) K = $(fmt(K)): pairs excluded as unresolved ($excl of $tot)", excl / tot, "<= 0.10", pf(excl / tot <= 0.10))
    end
    # information: the oracle is blind below about 1e-15 (it returns 0 or noise where our log form is still exact)
    c = vcase(:five, Float64, 1000)
    g30 = [ct_gali(5.0, c.x0[j], c.W0[j], 30)[end] for j in 1:100]
    info("AC16", "K = 5, n = 30, 100 starts: oracle value (min; median; max) while our exp(L) is ~1e-25: the oracle floor (information)", string(fmt(minimum(g30)), "; ", fmt(sl_median(g30)), "; ", fmt(maximum(g30))))
end

# ---------------------------------------------------------------------------------------------------
# AC17: strong chaos and precision at long horizon
section("AC17") do
    for (ci, label) in enumerate((:five, :ten))
        m, n = 200, 10_000
        seed = vseed(17, ci)
        starts = sl_starts(StableRNG(seed), m)
        res = Dict{DataType, Vector{Float64}}()
        for T in SL_TYPES
            sys = StandardMap(sl_K(T, label))
            res[T] = map(1:m) do j
                _, _, st = propagate_orbit(sys, sl_cast(T, starts[j]), sl_frame(seed, j, T), n, 1, (SALI(),); orbit = j)
                -Float64(reducer_final(SALI(), st[1]).L) / (2n)
            end
        end
        a, b = res[Float64], res[Double64]
        ratio_line("AC17", "K = $(Float64(sl_K(Float64, label))), 200 starts, n = 1e4: |mean(-L/(2n)) Float64 - Double64| / (4 combined SE); means $(fmt(mean(a))) (SE $(fmt(std(a) / sqrt(m)))) and $(fmt(mean(b))) (SE $(fmt(std(b) / sqrt(m)))); z = $(fmt(zdiff(a, b)))", zdiff(a, b) / 4)
        info("AC17", "K = $(Float64(sl_K(Float64, label))): Chirikov rate ln(K/2) (loose) against mean(-L/(2n))", string(fmt(log(sl_K(Float64, label) / 2)), " vs ", fmt(mean(b))))
    end
end

# ---------------------------------------------------------------------------------------------------
# Figures
section("figure") do
    # decay laws
    fig = CM.Figure(size = (1500, 480))
    seed = vseed(20)
    ax = CM.Axis(fig[1, 1]; xlabel = "n", ylabel = "log10 of SALI or GALI_2", title = "chaotic orbit, K = 2 pi: direct form hits the floor, log form does not")
    z0 = sl_starts(StableRNG(seed), 1)[1]
    X = 0.0
    runs = Dict{DataType, Any}()
    for T in SL_TYPES
        sys = StandardMap(sl_K(T, :two_pi))
        probe = sl_Probe(SALI(); all = true)
        _, _, st = propagate_orbit(sys, sl_cast(T, z0), sl_frame(seed, 1, T), 150, 1, (probe,); orbit = 1)
        runs[T] = reducer_final(probe, st[1]).snaps
    end
    _, _, stl = propagate_orbit(StandardMap(sl_K(Float64, :two_pi)), z0, SMatrix{2, 1, Float64}(sl_frame(seed, 1, Float64)[:, 1]), 10_000, 1, (Benettin(),); orbit = 1)
    X = (reducer_final(Benettin(), stl[1]).sum[1] + reducer_final(Benettin(), stl[1]).comp[1]) / 10_000
    ns = 1:150
    lg(x, fl) = log10(max(Float64(x), fl))
    CM.scatterlines!(ax, ns, [lg(s.sali, 1e-40) for s in runs[Float64]]; color = :crimson, markersize = 4, label = "SALI direct, Float64 (floor 4 eps = 8.9e-16)")
    CM.scatterlines!(ax, ns, [lg(s.sali, 1e-40) for s in runs[Double64]]; color = :navy, markersize = 4, label = "SALI direct, Double64 (floor 3.2e-30)")
    CM.lines!(ax, ns, [Float64(s.L) / log(10) for s in runs[Float64]]; color = :black, linestyle = :dash, label = "log form L/ln 10, Float64")
    CM.lines!(ax, ns, [Float64(s.L) / log(10) for s in runs[Double64]]; color = :gray50, linestyle = :dot, label = "log form L/ln 10, Double64")
    CM.lines!(ax, ns, [-2 * X * n / log(10) for n in ns]; color = :seagreen, label = "exp(-2 lambda_1 n), lambda_1 = X(1e4) = $(fmt(X)) of this orbit")
    CM.ylims!(ax, -80, 1)
    CM.axislegend(ax; position = :lb, labelsize = 9)
    # regular decay law
    ax2 = CM.Axis(fig[1, 2]; xscale = log10, yscale = log10, xlabel = "n", ylabel = "SALI (from the log form)", title = "regular orbits, K = 0.5: power law n^-2 of A9")
    c = vcase(:half, Float64, 20)
    nn = unique(round.(Int, 10 .^ range(1, 6; length = 40)))
    run = sl_run(SALI(), c.sys, c.x0[1:6], c.W0[1:6], 1_000_000, 1; steps = nn)
    for r in run
        CM.lines!(ax2, Float64.(nn), [Float64(sali_from_log(s.L)) for s in r.snaps]; linewidth = 1)
    end
    CM.lines!(ax2, [1e3, 1e6], [1e3^-2 * 30, 1e6^-2 * 30]; color = :black, linestyle = :dash, label = "30 n^-2 (slope -2)")
    CM.axislegend(ax2; position = :lb, labelsize = 9)
    # ensemble convergence of -L/(2n) to lambda_1
    ax3 = CM.Axis(fig[1, 3]; xscale = log10, xlabel = "n", ylabel = "ensemble mean of -L(n)/(2n)", title = "log form against lambda_1 (200 starts per K)")
    for (K, col) in ((:five, :seagreen), (:two_pi, :darkorange), (:ten, :purple))
        sd = vseed(21, findfirst(==(K), (:five, :two_pi, :ten)))
        sys = StandardMap(sl_K(Float64, K))
        x0 = [sl_cast(Float64, s) for s in sl_starts(StableRNG(sd), 200)]
        W0 = [sl_frame(sd, j, Float64) for j in 1:200]
        nlist = unique(round.(Int, 10 .^ range(1, 3.7; length = 30)))
        r = sl_run(SALI(), sys, x0, W0, maximum(nlist), 1; steps = nlist)
        l1 = sl_run(Benettin(), sys, x0, [SMatrix{2, 1, Float64}(W[:, 1]) for W in W0], 10_000, 1; steps = [10_000])
        lam = mean((q.snaps[1].sum[1] + q.snaps[1].comp[1]) / 10_000 for q in l1)
        CM.lines!(ax3, Float64.(nlist), [mean(-Float64(q.snaps[i].L) / (2 * nlist[i]) for q in r) for i in eachindex(nlist)]; color = col, label = "K = $(Float64(sl_K(Float64, K)))")
        CM.hlines!(ax3, [lam]; color = col, linestyle = :dash)
    end
    CM.axislegend(ax3; position = :rt, labelsize = 9)
    CM.save(joinpath(OUT, "sali-decay-laws.png"), fig)
    say("# figure saved: sali-decay-laws.png")
    # ChaosTools comparison
    haskey(DATA, :ct) || return
    fig2 = CM.Figure(size = (1200, 480))
    pts = [p for p in DATA[:ct] if p.resolved]
    ax4 = CM.Axis(fig2[1, 1]; xscale = log10, yscale = log10, xlabel = "g_oracle = ChaosTools gali(k = 2)", ylabel = "g = exp(L), Chaoslab", title = "AC16: GALI_2, resolved points")
    for (K, col) in ((0.5, :seagreen), (Float64(sl_K(Float64, :two_pi)), :darkorange), (5.0, :purple))
        sel = [p for p in pts if p.K == K && p.go > 0 && p.gC > 0]
        CM.scatter!(ax4, [p.go for p in sel], [p.gC for p in sel]; color = (col, 0.5), markersize = 4, label = "K = $(fmt(K))")
    end
    CM.lines!(ax4, [1e-12, 1], [1e-12, 1]; color = :black)
    CM.axislegend(ax4; position = :lt, labelsize = 9)
    ax5 = CM.Axis(fig2[1, 2]; xscale = log10, yscale = log10, xlabel = "g_oracle", ylabel = "|g_Chaoslab - g_oracle| and the bound B_SALI + B_oracle", title = "AC16: difference against its bound")
    for (K, col) in ((0.5, :seagreen), (Float64(sl_K(Float64, :two_pi)), :darkorange), (5.0, :purple))
        sel = [p for p in pts if p.K == K && p.go > 0]
        CM.scatter!(ax5, [p.go for p in sel], [max(abs(p.gC - p.go), 1e-20) for p in sel]; color = (col, 0.5), markersize = 4, label = "K = $(fmt(K)): difference")
        CM.scatter!(ax5, [p.go for p in sel], [p.bound for p in sel]; color = col, marker = :xcross, markersize = 5, label = "K = $(fmt(K)): bound")
    end
    CM.axislegend(ax5; position = :lt, labelsize = 8)
    CM.save(joinpath(OUT, "sali-chaostools.png"), fig2)
    say("# figure saved: sali-chaostools.png")
end

say(@sprintf("# verdicts: %d PASS, %d FAIL, %d other", count(==("PASS"), VERDICTS), count(==("FAIL"), VERDICTS), count(v -> v != "PASS" && v != "FAIL", VERDICTS)))
say(@sprintf("# total wall time %.1f s", time() - T0))
close(LOG)
