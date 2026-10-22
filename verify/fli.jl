# Stage 04 verification of module fli-ofli (reducer FLI, check_frame, fli_unlogged; layer Indicators; catalog row I5, OFLI is not in this module, D71).
# Run from the repo root:  julia -t 2 --project=verify verify/fli.jl
# Spec: workspaces/module-build/stages/01-spec/output/fli-ofli-spec.md (criteria 1 to 16, Amendment D72). Criterion 16 (CPU vs GPU) is in verify/gpu.jl, tag "fli".
# Section "suite": the package suite at full size (FL_FULL_RUN = true: every multiple of 10 as a checkpoint, n = 1e6 Float64 shear, every n of AC5).
# Other sections: one printed line per criterion  AC<n> | quantity | value | bound | PASS/FAIL  with the MEASURED number, reusing the helpers of
# test/fli-ofli/fl_helpers.jl (one home for the two-point bound B_FLI, B_S, the BigFloat reference, Psi_n) and test/sali/sl_helpers.jl (B_SALI of the
# oracle line only). Values are "worst ratio error / bound" (pass when <= 1) unless stated. The cases use verify seeds (FL_SEED + 58000 + 50 ac + sub),
# distinct from the package tests (FL_SEED + 33000 ..); the K = 0.5 regular class is screened again here (Benettin X(1e5) <= 5e-4, as L1 S7) with a verify seed.
# FL_ONLY="AC15,figure" runs only the named sections (debugging aid) and FL_LOG names the log file; the default runs everything.
# FL_CHILD=1 is the child-process mode of the thread-invariance check of AC13(c).
using Printf, Test, Random, LinearAlgebra, Statistics
using Chaoslab, StaticArrays, DoubleFloats, StableRNGs

# ---------------------------------------------------------------------------------------------------
# Thread-count invariance (AC13c): the same outputs hashed in this process and in child processes with other thread counts.
function v_flatten!(v, x)
    if x isa Number
        push!(v, Float64(x))
    elseif x isa Union{Tuple, NamedTuple, AbstractArray}
        foreach(y -> v_flatten!(v, y), x)
    else
        error("v_flatten!: $(typeof(x))")
    end
    return v
end
function thread_hashes()
    out = Dict{String, String}()
    seed = 20261012 + 58000 + 50 * 13 + 3
    for (cname, reducers) in (("FLI", (FLI(),)), ("FLI,SALI", (FLI(replace = false), SALI()))), tau in (1, 4), chunk in (64, 1000)
        rng = StableRNG(seed)
        x0 = [SVector(2pi * rand(rng), 2pi * rand(rng)) for _ in 1:4096]
        res = propagate(StandardMap(2pi), x0, 1000; k = 2, tau = tau, seed = seed, reducers = reducers, chunk = chunk, backend = :cpu)
        v = Float64[]
        v_flatten!(v, res.x); v_flatten!(v, res.W); v_flatten!(v, res.finals)
        out["$cname tau=$tau chunk=$chunk"] = string(hash(reinterpret(UInt64, v)))
    end
    return out
end
if get(ENV, "FL_CHILD", "") == "1"
    for (k, h) in sort!(collect(thread_hashes())); println("HASH ", k, " = ", h); end
    println("THREADS ", Threads.nthreads())
    exit(0)
end

const FL_FULL_RUN = true
const ROOT = dirname(@__DIR__)
const TESTDIR = joinpath(ROOT, "test", "fli-ofli")
include(joinpath(ROOT, "test", "sali", "sl_helpers.jl"))     # B_SALI of the oracle line (AC15) and lb_* (loaded first: lb_helpers.jl includes lb_psi.jl unguarded)
include(joinpath(TESTDIR, "fl_helpers.jl"))

vseed(ac, sub = 0) = FL_SEED + 58000 + 50 * ac + sub

import ChaosTools as CT
import CairoMakie as CM

const OUT = joinpath(ROOT, "workspaces", "module-build", "stages", "04-verify", "output")
mkpath(OUT)
const ONLY = filter(!isempty, split(get(ENV, "FL_ONLY", ""), ","))
const LOG = open(get(ENV, "FL_LOG", joinpath(OUT, "fli-verify-log.txt")), "w")
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
rat(x, B) = x <= 0 ? 0.0 : (B == 0 ? Inf : x / B)          # ratio of a (signed) excess to its bound; no excess is 0

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
const EPS64 = fl_eps(Float64)
say("# fli-ofli verification, julia $(VERSION), threads $(Threads.nthreads()), base seed $FL_SEED (verify offset 58000), time $(round(time())) (unix)")
say("# timings in this log were taken on a machine that was not otherwise idle")

# ---------------------------------------------------------------------------------------------------
# Cases with verify seeds (the package helpers fl_case / fl_rows use the package seeds)
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
            append!(cands, fl_starts(rng, 60))
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
        K = fl_K(T, label)
        if label === :half
            r = vregular(50)
            ids, starts, fseed = r.ids[1:m], r.starts[1:m], r.seed
        else
            fseed = vseed(0, findfirst(==(label), (:two_pi, :ten, :five)))
            ids = collect(1:m)
            starts = fl_starts(StableRNG(fseed), m)
        end
        (sys = StandardMap(K), K = K, ids = ids, x0 = [fl_cast(T, s) for s in starts], fseed = fseed, W0 = [fl_frame(fseed, id, T) for id in ids])
    end
end

# BigFloat reference products, Float64 tangents (for Psi) and frame of every start of a case (as fl_orbit_ref, verify seeds).
function vorbit(label::Symbol, ::Type{T}, m::Int, N::Int; bits = FL_BITS) where {T}
    return get!(VC, (:orb, label, T, m, N, bits)) do
        c = vcase(label, T, m)
        per = Vector{fl_Orbit}(undef, m)
        Threads.@threads :dynamic for j in 1:m
            pts = fl_pts(c.sys, c.x0[j], N)
            R = fl_ref_steps(fl_mats_std(c.K, pts), c.W0[j], N; bits = bits)
            per[j] = fl_Orbit(R.q1, R.q2, R.G, lb_tangents(c.sys, c.x0[j], N), c.W0[j], Dict{Tuple{Int, Int}, Float64}(), Dict{Int, Vector{NamedTuple}}())
        end
        (per = per, case = c, bits = bits)
    end
end

# One row per (start, checkpoint): the run, the BigFloat reference, the two-point bounds B_FLI of Q11 for both forms and the ratios error / bound.
function vrows(label::Symbol, ::Type{T}, tau::Int, m::Int, N::Int; bits = FL_BITS, cps::Vector{Int} = fl_checkpoints(N)) where {T}
    return get!(VC, (:rows, label, T, tau, m, N, bits, cps)) do
        O = vorbit(label, T, m, N; bits = bits)
        c = O.case
        hs = filter(n -> n % tau == 0, cps)
        R = [r[1] for r in fl_run((FLI(),), c.sys, c.x0, c.W0, N, tau; steps = hs)]
        epsT = fl_eps(T)
        rows = NamedTuple[]
        for j in 1:m
            o = O.per[j]
            refs = fl_ref_at(o, tau, hs; bits = bits)
            for (s, r) in zip(R[j].snaps, refs)
                s.n == r.n || error("vrows: run and reference disagree on the step ($(s.n) vs $(r.n))")
                B1 = fl_BFLI(epsT, o, 1, s.step1, 1, r.j1)
                BB = fl_BFLI(epsT, o, s.colB, s.stepB, r.cB, r.jB)
                e1, eB = fl_err(s.fli1, r.F1; bits = bits), fl_err(s.fli, r.FB; bits = bits)
                push!(rows, (j = j, n = s.n, fli1 = s.fli1, fli = s.fli, S1 = s.sum[1] + s.comp[1], S2 = s.sum[2] + s.comp[2], F1 = r.F1, FB = r.FB, Fstar = r.Fstar,
                             e1 = e1, eB = eB, B1 = B1, BB = BB, ratio1 = rat(e1, B1), ratioB = rat(eB, BB), maxdiff = (s.stepB, s.colB) != (r.jB, r.cB)))
            end
        end
        rows
    end
end

vlabel(label) = label === :half ? "0.5 regular class" : label === :two_pi ? "2 pi" : label === :ten ? "10" : "5"
const VPANEL = ((:half, 10_000, 20), (:two_pi, 10_000, 20), (:ten, 10_000, 20))
mean_se(v) = (mean(v), std(v) / sqrt(length(v)))
zdiff(a, b) = abs(mean(a) - mean(b)) / sqrt(var(a) / length(a) + var(b) / length(b))

# Ensemble of FLI_B(n)/n (or the first-column form) in type T on m verify starts (the same Float64 starts widened exactly).
function v_ens(::Type{T}, label::Symbol, m::Int, n::Int; first::Bool = false) where {T}
    c = vcase(label, T, m)
    return map(1:m) do j
        _, _, st = propagate_orbit(c.sys, c.x0[j], c.W0[j], n, 1, (FLI(),); orbit = j)
        f = reducer_final(FLI(), st[1])
        Float64(first ? f.fli1 : f.fli) / n
    end
end

# ---------------------------------------------------------------------------------------------------
# The package suite at full size
section("suite") do
    local ts
    t = @elapsed try
        ts = @testset verbose = true "fli-ofli verify" begin
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
    emit("SUITE", "package suite at full size (FL_FULL_RUN), assertions passed (fail + error = $nbad); wall time $(round(t; digits = 1)) s", npass, "0 fail, 0 error", pf(ok && nbad == 0 && npass > 0))
end

# ---------------------------------------------------------------------------------------------------
# AC1: algebra of the definition
function v_alg_exact(kind::Symbol, N::Int)
    return setprecision(BigFloat, FL_BITS) do
        P = Vector{BigFloat}(undef, N)
        F1, FB = similar(P), similar(P)
        acc = BigFloat(0)
        r1 = rB = BigFloat(0)
        for m in 1:N
            acc += BigFloat(fl_c(kind, m))
            r1 = max(r1, acc)
            rB = max(rB, acc, -acc)
            P[m], F1[m], FB[m] = acc, r1, rB
        end
        (P = P, F1 = F1, FB = FB)
    end
end
v_alg_ratio(f1, fB, ex, B) = maximum(max(rat(fl_err(f1[j], ex.F1[j]), B[j]), rat(fl_err(fB[j], ex.FB[j]), B[j])) for j in eachindex(f1))

function ac1(T)
    epsT = fl_eps(T)
    N = 40
    worst = 0.0
    bad_bits = bad_sorted = bad_exch = bad_sign = 0
    spec_ok = true
    ctl_pass = 0        # negative controls that wrongly satisfy the predicate (must be 0)
    nctl = 0
    wspec = 0.0         # kind (ii): FLI_B - 0.2 j and FLI_1 excess
    for kind in (:a, :b, :c)
        ex = v_alg_exact(kind, N)
        B = [8 * epsT * (j + 2 + exp(2 * abs(fl_c(kind, 1)))) for j in 1:N]
        sys = fl_Diag(kind)
        x0 = SVector{2, T}(0, 0)
        I2, P2 = SMatrix{2, 2, T}(1, 0, 0, 1), SMatrix{2, 2, T}(0, 1, 1, 0)
        snaps(W) = fl_run((FLI(),), sys, [x0], [W], N, 1; all = true)[1][1].snaps
        sn = snaps(I2)
        f1, fB = [s.fli1 for s in sn], [s.fli for s in sn]
        S1, S2 = [s.sum[1] + s.comp[1] for s in sn], [s.sum[2] + s.comp[2] for s in sn]
        bad_bits += !(isequal(f1, accumulate(max, S1; init = zero(T))) && isequal(fB, accumulate(max, max.(S1, S2); init = zero(T))))
        bad_sorted += !(issorted(f1) && issorted(fB))
        worst = max(worst, v_alg_ratio(f1, fB, ex, B))
        if kind === :a
            spec_ok &= abs(ex.FB[10] - 3) < 1e-14 && abs(ex.FB[22] - 3) < 1e-14 && ex.F1[16] < 1e-14 && abs(ex.F1[17] - 0.5) < 1e-14
            spec_ok &= all(j -> f1[j] <= B[j], 1:15) && f1[17] > 0.4                                  # D72 item 4: FLI_1 = 0 asserted up to j = 15
            wspec = max(wspec, Float64(maximum(f1[j] / B[j] for j in 1:15)))
        elseif kind === :b
            spec_ok &= all(j -> 0 <= f1[j] <= B[j], 1:N)
            wspec = max(wspec, Float64(maximum(f1[j] / B[j] for j in 1:N)), maximum(rat(fl_err(fB[j], BigFloat(0.2) * j), B[j]) for j in 1:N))
        end
        sw = snaps(P2)
        bad_exch += !(isequal([s.fli for s in sw], fB) && isequal([s.fli1 for s in sw], accumulate(max, S2; init = zero(T))))
        kind === :b && (spec_ok &= [s.fli1 for s in sw][end] > 0.1 && f1[end] <= B[end])
        sg = snaps(SMatrix{2, 2, T}(-1, 0, 0, 1))
        bad_sign += !(isequal([s.fli for s in sg], fB) && isequal([s.fli1 for s in sg], f1))
        # negative controls: each double must violate the predicate (ratio above 1)
        ctls = kind !== :b ? [v_alg_ratio(S1, max.(S1, S2), ex, B)] : [v_alg_ratio(accumulate(max, S1), fB, ex, B)]
        push!(ctls, v_alg_ratio(f1, exp.(fB), ex, B))
        ctl_pass += count(<=(1), ctls)
        nctl += length(ctls)
    end
    return (; worst, bad_bits, bad_sorted, bad_exch, bad_sign, spec_ok, ctl_pass, nctl, wspec)
end

section("AC1") do
    for T in FL_TYPES
        r = ac1(T)
        ratio_line("AC1", "max |FLI_1, FLI_B - exact running maxima (BigFloat of the stored c_m)| / B_S(j), maps (i), (ii), (iii), 40 handoffs, $T", r.worst)
        emit("AC1", "FLI_1, FLI_B equal the running maxima of the reducer's own S_1, max(S_1, S_2) with the 0, bitwise: failures over 3 maps, $T", r.bad_bits, "0", pf(r.bad_bits == 0))
        emit("AC1", "nondecreasing at every handoff: failures over 3 maps, $T", r.bad_sorted, "0", pf(r.bad_sorted == 0))
        emit("AC1", "column exchange: FLI_B bitwise unchanged, FLI_1 = running maximum of the other sum: failures, $T", r.bad_exch, "0", pf(r.bad_exch == 0))
        emit("AC1", "sign change of a column: both forms bitwise unchanged: failures, $T", r.bad_sign, "0", pf(r.bad_sign == 0))
        emit("AC1", "printed values of the spec (FLI_B = 3 from j = 10 to 22, FLI_1 = 0 to j = 15 and 0.5 at j = 17, FLI_1 = 0 and FLI_B = 0.2 j for (ii), FLI_1 asymmetric), $T; worst FLI_1 / B_S there $(fmt(r.wspec))", r.spec_ok, "true", pf(r.spec_ok && r.wspec <= 1))
        emit("AC1", "negative controls (no memory, no j = 0 term, no log) that wrongly satisfy the predicate (of $(r.nctl)), $T", r.ctl_pass, "0", pf(r.ctl_pass == 0))
    end
end

# ---------------------------------------------------------------------------------------------------
# AC2: K = 0 shear (Q9a)
function v_shear_q(W0, n::Integer)
    return setprecision(BigFloat, FL_BITS) do
        a1, b1, a2, b2 = fl_big(W0[1, 1]), fl_big(W0[2, 1]), fl_big(W0[1, 2]), fl_big(W0[2, 2])
        ((a1 + n * b1)^2 + b1^2, (a2 + n * b2)^2 + b2^2)
    end
end
function v_cf(q1, q2)
    return setprecision(BigFloat, FL_BITS) do
        one_ = one(BigFloat)
        (F1 = log(max(one_, q1)) / 2, FB = log(max(one_, q1, q2)) / 2)
    end
end

function ac2(T)
    seed = vseed(2)
    rng = StableRNG(seed)
    epsT = fl_eps(T)
    sys = StandardMap(exact(T, 0))
    M = SMatrix{2, 2, Float64}(1, 0, 1, 1)
    w1 = wB = wI = 0.0
    nI = badI = 0
    for tau in (1, 10), (case, nid) in (("seeded", 20), ("identity", 1)), N in (10_000, 1_000_000)
        (N == 1_000_000 && T !== Float64) && continue
        m = case == "identity" ? 1 : (N == 1_000_000 ? 3 : nid)
        x0 = [fl_cast(T, SVector(2pi * rand(rng), 2pi * rand(rng))) for _ in 1:m]
        W0 = case == "identity" ? [SMatrix{2, 2, T}(1, 0, 0, 1)] : [fl_frame(seed, j, T) for j in 1:m]
        ns = N == 1_000_000 ? [1_000_000] : [10, 100, 10_000]
        runs = fl_run((FLI(),), sys, x0, W0, N, tau; steps = ns)
        B(n) = 8 * epsT * (n + fl_Gfirst(M, tau))
        for j in 1:m, s in runs[j][1].snaps
            n = s.n
            cf = v_cf(v_shear_q(W0[j], n)...)
            w1 = max(w1, fl_err(s.fli1, cf.F1) / B(n))
            wB = max(wB, fl_err(s.fli, cf.FB) / B(n))
            if case == "identity"
                nI += 1
                badI += s.fli1 == 0 ? 0 : 1
                wI = max(wI, fl_err(s.fli, setprecision(() -> log(BigFloat(n)^2 + 1) / 2, BigFloat, FL_BITS)) / B(n))
            end
        end
    end
    return (; w1, wB, wI, nI, badI)
end

section("AC2") do
    for T in FL_TYPES
        r = ac2(T)
        ratio_line("AC2", "K = 0 shear, seeded frame and W0 = I, tau in {1,10}, n in {10,100,1e4$(T === Float64 ? ",1e6" : "")}: max |FLI_1 - FLI_cf| / B_S(n), $T", r.w1)
        ratio_line("AC2", "same, max |FLI_B - FLI_cf| / B_S(n), $T", r.wB)
        emit("AC2", "W0 = I: FLI_1 != 0 (tolerance 0) over $(r.nI) handoff values, $T", r.badI, "0", pf(r.badI == 0))
        ratio_line("AC2", "W0 = I: max |FLI_B - (1/2) ln(n^2 + 1)| / B_S(n), $T", r.wI)
    end
end

# ---------------------------------------------------------------------------------------------------
# AC3 to AC5: analytic cases
function v_y5_eigen()
    mu = (3 + sqrt(BigFloat(5))) / 2
    phi = (1 + sqrt(BigFloat(5))) / 2
    nv = sqrt(phi^2 + 1)
    return (mu = mu, vp = (phi / nv, 1 / nv), vm = (-1 / nv, phi / nv))
end
function v_y5_coeffs(W0)
    return setprecision(BigFloat, FL_BITS) do
        e = v_y5_eigen()
        w = ((fl_big(W0[1, 1]), fl_big(W0[2, 1])), (fl_big(W0[1, 2]), fl_big(W0[2, 2])))
        (c1 = ntuple(i -> w[i][1] * e.vp[1] + w[i][2] * e.vp[2], 2), c2 = ntuple(i -> w[i][1] * e.vm[1] + w[i][2] * e.vm[2], 2))
    end
end
function v_y5_frame(c1::Real, ::Type{T}) where {T}
    return setprecision(BigFloat, FL_BITS) do
        e = v_y5_eigen()
        c = BigFloat(c1)
        s = sqrt(1 - c^2)
        w1 = (c * e.vp[1] + s * e.vm[1], c * e.vp[2] + s * e.vm[2])
        w2 = (-w1[2], w1[1])
        SMatrix{2, 2, T}(T(w1[1]), T(w1[2]), T(w2[1]), T(w2[2]))
    end
end
# B_S(n) = 8 eps (n + G_first); with `psi`, max(B_S, 8 eps Psi_n^(i)) (L1 Amendment D35 per vector, D72 item 2).
function v_lin_B(epsT, M::SMatrix, Ms, W0, n::Int, tau::Int, memo, psi::Bool = false)
    psi || return 8 * epsT * (n + fl_Gfirst(M, tau))
    p1, p2 = get!(memo, n) do
        (lb_Psi(view(Ms, 1:n), W0[:, 1]), lb_Psi(view(Ms, 1:n), W0[:, 2]))
    end
    return 8 * epsT * max(n + fl_Gfirst(M, tau), p1, p2)
end

function ac3(T)
    seed = vseed(3)
    epsT = fl_eps(T)
    nid = 100
    sys = LinearSymplecticMap(SMatrix{2, 2, T}(2, 1, 1, 1))
    W0 = vcat([fl_frame(seed, j, T) for j in 1:nid], [v_y5_frame(1e-3, T), v_y5_frame(1e-6, T)])
    nfr = length(W0)
    x0 = [SVector{2, T}(0, 0) for _ in 1:nfr]
    M64 = SMatrix{2, 2, Float64}(2, 1, 1, 1)
    Ms = lb_tangents(sys, x0[1], 200)
    refs = [fl_ref_steps(fl_mats_const(2, 1, 1, 1), W0[j], 200) for j in 1:nfr]
    memo = [Dict{Int, Tuple{Float64, Float64}}() for _ in 1:nfr]
    cs = [v_y5_coeffs(W0[j]) for j in 1:nfr]
    excluded = [j for j in 1:nid if abs(cs[j].c1[1]) < 0.05]
    helper = setprecision(BigFloat, FL_BITS) do
        Y8 = [BigFloat(2) 1; 1 1]^8
        maximum(1:3) do j
            q = fl_qcols(Y8, W0[j])
            Float64(max(abs(q[1] - refs[j].q1[8]) / q[1], abs(q[2] - refs[j].q2[8]) / q[2]))
        end
    end
    frames_ok = abs(abs(cs[nid + 1].c1[1]) - 1e-3) < 1e-12 && abs(abs(cs[nid + 2].c1[1]) - 1e-6) < 1e-12 && abs(cs[nid + 1].c1[2]) > 0.999
    wa = [0.0, 0.0, 0.0, 0.0]        # seeded B_S, c1 = 1e-3 (Psi form), c1 = 1e-6 (Psi form), c1 = 1e-6 plain B_S (reported)
    wc = 0.0
    bok = true
    for tau in (1, 5)
        ns = tau == 1 ? [1, 2, 4, 8, 20, 100, 200] : [5, 10, 20, 100, 200]
        run = fl_run((FLI(),), sys, x0, W0, 200, tau; steps = ns)
        for j in 1:nfr
            j in excluded && continue
            at = fl_ref_at(refs[j], tau, ns)
            for (s, r) in zip(run[j][1].snaps, at)
                n = s.n
                e = max(fl_err(s.fli1, r.F1), fl_err(s.fli, r.FB))
                B = v_lin_B(epsT, M64, Ms, W0[j], n, tau, memo[j], j > nid)
                k = j <= nid ? 1 : j == nid + 1 ? 2 : 3
                wa[k] = max(wa[k], e / B)
                j == nid + 2 && (wa[4] = max(wa[4], e / v_lin_B(epsT, M64, Ms, W0[j], n, tau, memo[j], false)))
                if n in (20, 100)
                    asym, tail = setprecision(BigFloat, FL_BITS) do
                        mu = v_y5_eigen().mu
                        cmax = max(abs(cs[j].c1[1]), abs(cs[j].c1[2]))
                        (n * log(mu) + log(cmax), Float64(mu^(-4n) * (1 + 1 / cmax^2)))
                    end
                    wc = max(wc, fl_err(s.fli, asym) / (B + tail))
                end
                if j == nid + 1 && tau == 1 && n == 4
                    bok &= r.F1 == 0 && s.fli1 == 0 && s.fli > 0.5
                end
            end
        end
    end
    return (; wa, wc, bok, excluded = length(excluded), helper, frames_ok)
end

function ac4(T, K)
    seed = vseed(4)
    epsT = fl_eps(T)
    nid = 100
    sys = StandardMap(exact(T, K))
    W0 = [fl_frame(seed, j, T) for j in 1:nid]
    x0 = [SVector{2, T}(0, 0) for _ in 1:nid]
    M64 = SMatrix{2, 2, Float64}(1 + K, K, 1, 1)
    Ms = lb_tangents(sys, x0[1], 200)
    refs = [fl_ref_steps(fl_mats_const(1 + K, 1, K, 1), W0[j], 200) for j in 1:nid]
    memo = [Dict{Int, Tuple{Float64, Float64}}() for _ in 1:nid]
    helper = setprecision(BigFloat, FL_BITS) do
        Y8 = [BigFloat(1 + K) 1; BigFloat(K) 1]^8
        maximum(1:3) do j
            q = fl_qcols(Y8, W0[j])
            Float64(max(abs(q[1] - refs[j].q1[8]) / q[1], abs(q[2] - refs[j].q2[8]) / q[2]))
        end
    end
    w = 0.0
    for tau in (1, 5)
        ns = tau == 1 ? [1, 2, 4, 8, 20, 100, 200] : [5, 10, 20, 100, 200]
        run = fl_run((FLI(),), sys, x0, W0, 200, tau; steps = ns)
        for j in 1:nid
            at = fl_ref_at(refs[j], tau, ns)
            for (s, r) in zip(run[j][1].snaps, at)
                w = max(w, max(fl_err(s.fli1, r.F1), fl_err(s.fli, r.FB)) / v_lin_B(epsT, M64, Ms, W0[j], s.n, tau, memo[j]))
            end
        end
    end
    return (; w, helper)
end

function ac5(T)
    seed = vseed(5)
    epsT = fl_eps(T)
    rot = LinearSymplecticMap(SMatrix{2, 2, T}(0, -1, 1, 0))
    nrot = nbad = 0
    for tau in (1, 3)
        N = tau == 1 ? 10_000 : 9_999
        snaps = fl_run((FLI(),), rot, [SVector{2, T}(0, 0)], [SMatrix{2, 2, T}(1, 0, 0, 1)], N, tau; all = true)[1][1].snaps
        nrot += length(snaps)
        nbad += count(s -> !(s.fli == 0 && s.fli1 == 0), snaps)
    end
    sys = StandardMap(exact(T, -2))
    nid = 100
    x0 = [SVector{2, T}(0, 0) for _ in 1:nid]
    W0 = [fl_frame(seed, j, T) for j in 1:nid]
    nsamp = collect(1:1000)
    M64 = SMatrix{2, 2, Float64}(-1, -2, 1, 1)
    Ms = lb_tangents(sys, x0[1], 1000)
    refs = [fl_ref_steps(fl_mats_const(-1, 1, -2, 1), W0[j], 1000) for j in 1:nid]
    memo = [Dict{Int, Tuple{Float64, Float64}}() for _ in 1:nid]
    ln0 = setprecision(() -> log((3 + sqrt(BigFloat(5))) / 2), BigFloat, FL_BITS)
    w1 = wB = wln = w2 = 0.0
    nonneg = true
    helper = 0.0
    for tau in (1, 2)
        ns = filter(n -> n % tau == 0, nsamp)
        run = fl_run((FLI(),), sys, x0, W0, 1000, tau; steps = ns)
        for j in 1:nid
            at = fl_ref_at(refs[j], tau, ns)
            cf = setprecision(BigFloat, FL_BITS) do
                l1, l2 = log(refs[j].q1[1]) / 2, log(refs[j].q2[1]) / 2
                e1 = log(sum(abs2, (fl_big(W0[j][1, 1]), fl_big(W0[j][2, 1])))) / 2
                e2 = log(sum(abs2, (fl_big(W0[j][1, 2]), fl_big(W0[j][2, 2])))) / 2
                (F1 = max(zero(l1), l1, e1), FB = max(zero(l1), l1, l2, e1, e2))
            end
            tau == 1 && (helper = max(helper, Float64(max(abs(at[end].F1 - cf.F1), abs(at[end].FB - cf.FB)))))
            for (s, r) in zip(run[j][1].snaps, at)
                B = v_lin_B(epsT, M64, Ms, W0[j], s.n, tau, memo[j])
                if tau == 1
                    w1 = max(w1, fl_err(s.fli1, cf.F1) / B)
                    wB = max(wB, fl_err(s.fli, cf.FB) / B)
                    wln = max(wln, rat(fl_diff(s.fli, ln0), B))
                else
                    nonneg &= s.fli1 >= 0 && s.fli >= 0
                    w2 = max(w2, Float64(max(s.fli1, s.fli)) / B)
                end
            end
        end
    end
    return (; nrot, nbad, w1, wB, wln, w2, nonneg, helper)
end

section("AC3-5") do
    for T in FL_TYPES
        r = ac3(T)
        emit("AC3", "linear hyperbolic map Y5: ids with |c_1| < 0.05 excluded (of 100), supplied frames |c_1| = 1e-3, 1e-6 as intended, BigFloat product against A^8 (max relative difference $(fmt(r.helper))), $T", r.excluded, "<= 10", pf(r.excluded <= 10 && r.frames_ok && r.helper < 1e-60))
        ratio_line("AC3", "(a) seeded ids, tau in {1,5}, n in {1,..,200}: max |FLI - FLI_cf| / B_S(n), both forms, $T", r.wa[1])
        ratio_line("AC3", "(a) supplied frame |c_1| = 1e-3: same over max(B_S, 8 eps Psi_n) (D72 item 2), $T", r.wa[2])
        ratio_line("AC3", "(a) supplied frame |c_1| = 1e-6: same over max(B_S, 8 eps Psi_n) (D72 item 2), $T", r.wa[3])
        info("AC3", "(a) supplied frame |c_1| = 1e-6: the same errors over the plain B_S (above 1 is why D72 item 2 was needed), $T", r.wa[4])
        emit("AC3", "(b) |c_1| = 1e-3, n = 4, tau = 1: FLI_1 == 0 exactly (exact reference 0) and FLI_B > 0.5, $T", r.bok, "true", pf(r.bok))
        ratio_line("AC3", "(c) rate at n in {20,100}: max |FLI_B - (n ln mu + ln |c_max|)| / (B + mu^(-4n) (1 + 1/c_max^2)), $T", r.wc)
    end
    for T in FL_TYPES, K in (1, 2, 6)
        r = ac4(T, K)
        ratio_line("AC4", "hyperbolic fixed point K = $K, 100 ids, tau in {1,5}, n in {1,..,200}: max |FLI - FLI_cf| / B_S(n), both forms (BigFloat product against the matrix power: $(fmt(r.helper))), $T", r.w)
    end
    for T in FL_TYPES
        r = ac5(T)
        emit("AC5", "(a) order-4 rotation, W0 = I, n to 1e4, tau in {1,3}: handoffs with FLI != 0 (of $(r.nrot)), $T", r.nbad, "0", pf(r.nbad == 0))
        ratio_line("AC5", "(b) K = -2 at (0,0), 100 ids, tau = 1, every n in 1..1000: max |FLI_1 - max(0, ln norm(M w_1))| / B_S(n), $T (helper against closed form: $(fmt(r.helper)))", r.w1)
        ratio_line("AC5", "(b) max |FLI_B - closed form| / B_S(n), $T", r.wB)
        ratio_line("AC5", "(b) max (FLI_B - ln((3 + sqrt 5)/2))_+ / B_S(n), $T", r.wln)
        ratio_line("AC5", "(c) tau = 2: max FLI / B_S(n) (0 <= FLI <= B_S, nonnegative: $(r.nonneg)), $T", r.nonneg ? r.w2 : Inf)
    end
end

# ---------------------------------------------------------------------------------------------------
# AC6 to AC8: same-orbit reference, basis-free sandwich, cross-check with L1 and sali
function v_single_fli(o::fl_Orbit, n::Int)
    return setprecision(BigFloat, FL_BITS) do
        acc1 = acc2 = mx = Float32(0)
        p1 = p2 = BigFloat(0)
        for i in 1:n
            s1, s2 = log(o.q1[i]) / 2, log(o.q2[i]) / 2
            acc1 += Float32(s1 - p1)
            acc2 += Float32(s2 - p2)
            p1, p2 = s1, s2
            mx = max(mx, acc1, acc2)
        end
        Float64(mx)
    end
end

section("AC6-7") do
    lnsq2 = setprecision(() -> log(BigFloat(2)) / 2, BigFloat, FL_BITS)
    for T in FL_TYPES, (label, N, m) in VPANEL
        K = vlabel(label)
        for tau in (1, 4)
            rows = vrows(label, T, tau, m, N)
            nexp = m * count(n -> n % tau == 0, fl_checkpoints(N))
            emit("AC6", "K = $K, $T, tau = $tau: rows (starts x checkpoints at every multiple of 10 up to $N)", length(rows), "$nexp", pf(length(rows) == nexp))
            ratio_line("AC6", "K = $K, $T, tau = $tau: max |FLI_1 - FLI_ref,1| / B_FLI (realized ratio R2), $m starts, n up to $N", maximum(r.ratio1 for r in rows))
            ratio_line("AC6", "K = $K, $T, tau = $tau: max |FLI_B - FLI_ref,B| / B_FLI (realized ratio R2)", maximum(r.ratioB for r in rows))
            lastn = maximum(r.n for r in rows)
            info("AC6", "K = $K, $T, tau = $tau: pairs with a differing maximiser (of $(length(rows))); at the final checkpoint (of $m)", string(count(r -> r.maxdiff, rows), "; ", count(r -> r.maxdiff && r.n == lastn, rows)))
            info("AC6", "K = $K, $T, tau = $tau: max |error| in units of eps / median ratio error/B_FLI (FLI_B)", string(fmt(maximum(r.eB for r in rows) / fl_eps(T)), " / ", fmt(fl_median([r.ratioB for r in rows]))))
        end
        if T === Float64
            rows = vrows(label, T, 1, m, N)
            if label !== :half
                O = vorbit(label, T, m, N)
                dmin = minimum(j -> (r = only(r for r in rows if r.j == j && r.n == 1000); abs(v_single_fli(O.per[j], 1000) - Float64(r.FB)) / r.BB), 1:3)
                emit("AC6", "K = $K: negative control, Float32 accumulator at n = 1e3, min |FLI_B(Float32) - ref| / B_FLI over 3 starts", dmin, "> 1", pf(dmin > 1))
            end
        end
        rows = vrows(label, T, 1, m, N)
        ratio_line("AC7", "K = $K, $T: max (FLI_B - F_*)_+ / B_FLI (sandwich upper)", maximum(rat(fl_diff(r.fli, r.Fstar), r.BB) for r in rows))
        ratio_line("AC7", "K = $K, $T: max (F_* - ln sqrt 2 - FLI_B)_+ / B_FLI (sandwich lower)", maximum(rat(fl_diff(r.fli, r.Fstar - lnsq2) * -1, r.BB) for r in rows))
        ratio_line("AC7", "K = $K, $T: max (FLI_1 - F_*)_+ / B_FLI", maximum(rat(fl_diff(r.fli1, r.Fstar), r.B1) for r in rows))
        u = [-fl_diff(r.fli, r.Fstar) / Float64(lnsq2) for r in rows]
        v = [-fl_diff(r.fli1, r.Fstar) for r in rows]
        info("AC7", "K = $K, $T: (F_* - FLI_B) / ln sqrt 2 (min; median; max) and F_* - FLI_1 (min; median; max)", string(fmt(minimum(u)), "; ", fmt(fl_median(u)), "; ", fmt(maximum(u)), " and ", fmt(minimum(v)), "; ", fmt(fl_median(v)), "; ", fmt(maximum(v))))
    end
end

section("AC8") do
    T = Float64
    epsT = EPS64
    N = 10_000
    cps = fl_checkpoints(N)
    for label in (:half, :two_pi, :ten)
        m = 20
        c = vcase(label, T, m)
        run = fl_run((FLI(replace = false), SALI()), c.sys, c.x0, c.W0, N, 1; steps = cps, all = (false, true))
        l1 = [fl_run((Benettin(),), c.sys, c.x0, [SMatrix{2, 1, T}(W[:, i]) for W in c.W0], N, 1; all = true) for i in 1:2]
        drawdown = Dict(n => Float64[] for n in (100, 1_000, 10_000))
        basis_gap = Dict(n => Float64[] for n in (100, 1_000, 10_000))
        bad_sali = bad_l1 = 0
        w_b = w_c1 = w_c2 = 0.0
        npts = 0
        for j in 1:m
            Ms = lb_tangents(c.sys, c.x0[j], N)
            fsn, ssn = run[j][1].snaps, run[j][2].snaps
            X1, X2 = l1[1][j][1].snaps, l1[2][j][1].snaps
            rm = accumulate(max, [i * X1[i].lambda[1] for i in 1:N]; init = 0.0)
            runL = accumulate(max, [-s.L / 2 for s in ssn])
            for s in fsn
                n = s.n
                sa = ssn[n]
                S1, S2 = s.sum[1] + s.comp[1], s.sum[2] + s.comp[2]
                bad_sali += !(S1 == sa.sum[1] + sa.comp[1] && S2 == sa.sum[2] + sa.comp[2] && s.sum == sa.sum && s.comp == sa.comp)
                bad_l1 += !(S1 == X1[n].sum[1] + X1[n].comp[1] && S2 == X2[n].sum[1] + X2[n].comp[1])
                w_b = max(w_b, abs(s.fli1 - rm[n]) / (4 * epsT * max(1, s.fli1)))
                p1, p2 = lb_Psi(view(Ms, 1:n), c.W0[j][:, 1]), lb_Psi(view(Ms, 1:n), c.W0[j][:, 2])
                BL = 8 * epsT * (p1 + p2)
                w_c1 = max(w_c1, rat(runL[n] - s.fli, BL / 2))
                w_c2 = max(w_c2, abs(max(S1, S2) - (-sa.L / 2 + abs(S1 - S2) / 2)) / BL)
                npts += 1
                if haskey(drawdown, n)
                    push!(drawdown[n], s.fli1 - n * X1[n].lambda[1])
                    push!(basis_gap[n], s.fli - s.fli1)
                end
            end
        end
        K = vlabel(label)
        emit("AC8", "(a) K = $K, Float64, $m starts x $(length(cps)) checkpoints: S_1, S_2 of the FLI reducer not bitwise equal to the sali sums (sum and comp separately): $bad_sali; not bitwise equal to the L1 sums of k = 1 runs from column 1 and 2: $bad_l1", bad_sali + bad_l1, "0", pf(bad_sali + bad_l1 == 0))
        ratio_line("AC8", "(b) K = $K: max |FLI_1(n) - running max of n X_L1| / (4 eps max(1, FLI_1)), $npts points", w_b)
        ratio_line("AC8", "(c) K = $K: max (max over handoffs of (-L/2) - FLI_B)_+ / (B_L/2)", w_c1)
        ratio_line("AC8", "(c) K = $K: max |max(S_1, S_2) - (-L/2 + |S_1 - S_2|/2)| / B_L", w_c2)
        for n in (100, 1_000, 10_000)
            info("AC8", "R1 K = $K, n = $n: drawdown D = FLI_1 - n X_L1 (min; median; max) and FLI_B - FLI_1 (median; max)",
                 string(fmt(minimum(drawdown[n])), "; ", fmt(fl_median(drawdown[n])), "; ", fmt(maximum(drawdown[n])), " and ", fmt(fl_median(basis_gap[n])), "; ", fmt(maximum(basis_gap[n]))))
        end
        if label !== :half
            # R1 label demo at n = 1e3: orbits whose label differs between "FLI_B above ln n" and "n X_L1 above ln n" (Maffione2011a thresholds, Q5)
            n = 1_000
            fb = [r[1].snaps[findfirst(s -> s.n == n, r[1].snaps)].fli for r in run]
            xl = [l1[1][j][1].snaps[n].lambda[1] * n for j in 1:m]
            info("AC8", "R1 K = $K, n = 1e3: orbits (of $m) with the label FLI_B > ln n / n X_L1 > ln n: $(count(>(log(n)), fb)) / $(count(>(log(n)), xl)); labels differ", count(j -> (fb[j] > log(n)) != (xl[j] > log(n)), 1:m))
        end
    end
    # the label demonstration for the regular class too
    c = vcase(:half, T, 20)
    run = fl_run((FLI(),), c.sys, c.x0, c.W0, 1_000, 1; steps = [1_000])
    l1 = fl_run((Benettin(),), c.sys, c.x0, [SMatrix{2, 1, T}(W[:, 1]) for W in c.W0], 1_000, 1; steps = [1_000])
    fb = [r[1].snaps[1].fli for r in run]
    xl = [r[1].snaps[1].lambda[1] * 1_000 for r in l1]
    info("AC8", "R1 K = 0.5 regular class, n = 1e3: orbits (of 20) with FLI_B > ln n / n X_L1 > ln n: $(count(>(log(1000)), fb)) / $(count(>(log(1000)), xl)); labels differ", count(j -> (fb[j] > log(1000)) != (xl[j] > log(1000)), 1:20))
end

# ---------------------------------------------------------------------------------------------------
# AC9: growth laws (ensemble level)
section("AC9") do
    m = 20
    c = vcase(:half, Float64, m)
    ns = unique(round.(Int, 10 .^ range(4, 6; length = 21)))
    run = fl_run((FLI(),), c.sys, c.x0, c.W0, 1_000_000, 1; steps = ns)
    first_decade = findall(<=(100_000), ns)
    transient = passed = below = 0
    slopes = Float64[]
    offsets = Float64[]
    for r in run
        sn = r[1].snaps
        slope = fl_slope(log.(Float64.(ns)), Float64[s.fli for s in sn])
        early = fl_slope(log.(Float64.(ns[first_decade])), Float64[sn[i].fli for i in first_decade])
        push!(slopes, slope)
        push!(offsets, fl_mean([sn[i].fli - log(ns[i]) for i in eachindex(ns)]))
        below += sn[1].fli < log(1e4)
        if !(0.7 <= early <= 1.3)
            transient += 1
        elseif 0.7 <= slope <= 1.3
            passed += 1
        end
    end
    emit("AC9", "(a) regular class K = 0.5, 20 starts, n in [1e4,1e6]: starts with the LS slope of FLI_B against ln n in [0.7,1.3] (slopes: min $(fmt(minimum(slopes))), median $(fmt(fl_median(slopes))), max $(fmt(maximum(slopes))))", passed, ">= $(ceil(Int, 0.8 * m))", pf(passed >= ceil(Int, 0.8 * m)))
    emit("AC9", "(a) starts in the transient (slope over [1e4,1e5] outside [0.7,1.3], D72 item 3)", transient, "<= $(floor(Int, 0.2 * m))", pf(transient <= 0.2 * m))
    info("AC9", "(a) reported: starts with FLI_B below ln n at n = 1e4 (the detector of the spec before D72), and FLI_B - ln n averaged over [1e4,1e6] per start (min; median; max)", string(below, " and ", fmt(minimum(offsets)), "; ", fmt(fl_median(offsets)), "; ", fmt(maximum(offsets))))
    DATA[:regular_slopes] = slopes
    mm = 1000
    sanity = Dict(:five => log(5 / 2), :ten => log(10 / 2), :two_pi => 1.172)
    for (ci, label) in enumerate((:two_pi, :five, :ten))
        seed = vseed(9, ci)
        sys = StandardMap(fl_K(Float64, label))
        x0 = [fl_cast(Float64, s) for s in fl_starts(StableRNG(seed), mm)]
        W0 = [fl_frame(seed, j, Float64) for j in 1:mm]
        w = collect(20:200)
        run = fl_run((FLI(),), sys, x0, W0, 1_000, 1; steps = vcat(w, [1_000]))
        l1 = fl_run((Benettin(),), sys, x0, [SMatrix{2, 1, Float64}(W[:, 1]) for W in W0], 10_000, 1; steps = [10_000])
        slopes = [fl_slope(Float64.(w), Float64[s.fli for s in r[1].snaps[1:181]]) for r in run]
        X = [(r[1].snaps[1].sum[1] + r[1].snaps[1].comp[1]) / 10_000 for r in l1]
        se = sqrt(var(slopes) / mm + var(X) / mm)
        d = abs(mean(slopes) - mean(X))
        K = vlabel(label)
        ratio_line("AC9", "(b) K = $K, 1000 starts: |mean(slope of FLI_B over n in [20,200]) - mean(X(1e4))| / (4 combined SE + 4/180); mean slope $(fmt(mean(slopes))), mean X $(fmt(mean(X)))", d / (4 * se + 4 / 180))
        info("AC9", "(b) K = $K: realized |difference| of the means; combined SE; z = difference / SE (allowance of the spec 4/180 = $(fmt(4 / 180)); the sali form 2/180 = $(fmt(2 / 180)))", string(fmt(d), "; ", fmt(se), "; ", fmt(d / se)))
        info("AC9", "(b) K = $K: R3 sanity rate (Chirikov ln(K/2); K = 2 pi: the benchmarks row 1.172, loose, K18) against mean slope", string(fmt(sanity[label]), " vs ", fmt(mean(slopes))))
        if label !== :two_pi
            a1 = [r[1].snaps[182].fli1 / 1_000 for r in run]
            aB = [r[1].snaps[182].fli / 1_000 for r in run]
            for (nm, a) in (("FLI_1", a1), ("FLI_B", aB))
                se2 = sqrt(var(a) / mm + var(X) / mm)
                ratio_line("AC9", "(c) K = $K: |mean($nm(1e3)/1e3) - mean(X(1e4))| / (4 combined SE + 4/n); means $(fmt(mean(a))) and $(fmt(mean(X)))", abs(mean(a) - mean(X)) / (4 * se2 + 4 / 1_000))
            end
        end
        hit = [(h = findfirst(s -> s.fli >= log(1e20), r[1].snaps); h === nothing ? nothing : r[1].snaps[h].n) for r in run]
        hv = [n for n in hit if n !== nothing]
        info("AC9", "R4 K = $K: first handoff at which FLI_B reaches ln 1e20 = 46.05 (min; median; max), and starts that never reach it within n = 1e3", string(minimum(hv), "; ", fl_median(hv), "; ", maximum(hv), " and ", count(isnothing, hit)))
    end
end

# ---------------------------------------------------------------------------------------------------
# AC10: tau dependence
section("AC10") do
    for T in FL_TYPES, (label, N, m) in VPANEL
        label === :ten && continue
        epsT = fl_eps(T)
        lns = log(fl_sigma(Float64(fl_K(T, label))))
        O = vorbit(label, T, m, N)
        base = Dict((r.j, r.n) => r for r in vrows(label, T, 1, m, N))
        for tau in (2, 5, 10)
            rows = [r for r in vrows(label, T, tau, m, N) if r.n in (100, 1_000)]
            wa = wb1 = wb2 = 0.0
            gap = Float64[]
            gapr = Float64[]
            for r in rows
                r1 = base[(r.j, r.n)]
                o = O.per[r.j]
                wa = max(wa, abs(Float64(r.S1 - r1.S1)) / (16 * epsT * fl_Psi(o, 1, r.n)), abs(Float64(r.S2 - r1.S2)) / (16 * epsT * fl_Psi(o, 2, r.n)))
                for (a, a1, B, B1) in ((r.fli, r1.fli, r.BB, r1.BB), (r.fli1, r1.fli1, r.B1, r1.B1))
                    wb1 = max(wb1, rat(Float64(a - a1), B + B1))
                    wb2 = max(wb2, rat(Float64(a1 - a) - (tau ÷ 2) * lns, B + B1))
                end
                push!(gap, Float64(r1.fli1 - r.fli1))
                push!(gapr, Float64(r1.fli1 - r.fli1) / ((tau ÷ 2) * lns))
            end
            K = vlabel(label)
            emit("AC10", "K = $K, $T, tau = $tau: rows (n in {100,1e3}, $m starts)", length(rows), "$(2m)", pf(length(rows) == 2m))
            ratio_line("AC10", "(a) K = $K, $T, tau = $tau: max |S_i^tau - S_i^1| / (16 eps Psi_n^(i))", wa)
            ratio_line("AC10", "(b) K = $K, $T, tau = $tau: max (FLI_tau - FLI_1)_+ / (B^tau + B^1), both forms", wb1)
            ratio_line("AC10", "(b) K = $K, $T, tau = $tau: max (FLI_1 - FLI_tau - floor(tau/2) ln sigma(K))_+ / (B^tau + B^1), both forms", wb2)
            info("AC10", "K = $K, $T, tau = $tau: FLI_1 - FLI_tau (median; max) and the same over the bound floor(tau/2) ln sigma = $(fmt((tau ÷ 2) * lns)) (max)", string(fmt(fl_median(gap)), "; ", fmt(maximum(gap)), " and ", fmt(maximum(gapr))))
        end
    end
end

# ---------------------------------------------------------------------------------------------------
# AC11: precision ladder and representation
section("AC11") do
    for label in (:two_pi, :five)
        a, b = v_ens(Float64, label, 1000, 1_000), v_ens(Double64, label, 1000, 1_000)
        ratio_line("AC11", "(a) K = $(vlabel(label)), 1000 starts, n = 1e3: |mean(FLI_B/n) Float64 - Double64| / (4 combined SE); means $(fmt(mean(a))) and $(fmt(mean(b))), z = $(fmt(zdiff(a, b)))", zdiff(a, b) / 4)
    end
    for label in (:five, :ten)
        a, b = v_ens(Float64, label, 200, 10_000), v_ens(Double64, label, 200, 10_000)
        ratio_line("AC11", "(a) K = $(vlabel(label)), 200 starts, n = 1e4: same; means $(fmt(mean(a))) and $(fmt(mean(b))), z = $(fmt(zdiff(a, b)))", zdiff(a, b) / 4)
    end
    setprecision(BigFloat, FL_BITS) do
        rows = vrows(:two_pi, BigFloat, 1, 50, 60; bits = 512)
        emit("AC11", "(b) BigFloat(256), K = 2 pi, 50 starts, n in {10,..,60}: rows", length(rows), "$(50 * length(fl_checkpoints(60)))", pf(length(rows) == 50 * length(fl_checkpoints(60))))
        ratio_line("AC11", "(b) max |FLI_1 - FLI_ref(512 bits)| / B_FLI (own eps $(fmt(fl_eps(BigFloat))))", maximum(r.ratio1 for r in rows))
        ratio_line("AC11", "(b) max |FLI_B - FLI_ref(512 bits)| / B_FLI", maximum(r.ratioB for r in rows))
        for label in (:two_pi, :five)
            a, b = v_ens(BigFloat, label, 100, 300), v_ens(Double64, label, 100, 300)
            ratio_line("AC11", "(c) K = $(vlabel(label)), 100 starts, n = 300: BigFloat(256) vs Double64 |difference of means| / (4 combined SE); means $(fmt(mean(a))) and $(fmt(mean(b))), z = $(fmt(zdiff(a, b)))", zdiff(a, b) / 4)
        end
    end
    O = vorbit(:ten, Float64, 5, 10_000)
    c = O.case
    runs = fl_run((FLI(),), c.sys, c.x0, c.W0, 10_000, 1; steps = [10_000])
    w = 0.0
    okd = true
    mn, mx = Inf, 0.0
    for j in 1:5
        s = runs[j][1].snaps[1]
        r = only(fl_ref_at(O.per[j], 1, [10_000]))
        okd &= isfinite(s.fli) && s.fli > 709.78 && fli_unlogged(s.fli) === nothing
        mn, mx = min(mn, s.fli), max(mx, s.fli)
        w = max(w, fl_err(s.fli, r.FB) / fl_BFLI(EPS64, O.per[j], s.colB, s.stepB, r.cB, r.jB))
    end
    emit("AC11", "(d) K = 10, n = 1e4, Float64, 5 starts: FLI_B finite and above 709.78 (range $(fmt(mn)) to $(fmt(mx))), not clipped, fli_unlogged == nothing (display reports the overflow)", okd, "true", pf(okd))
    ratio_line("AC11", "(d) max |FLI_B - FLI_ref| / B_FLI at n = 1e4, K = 10", w)
    okf = fli_unlogged(1.0) ≈ exp(1.0) && fli_unlogged(700.0) ≈ exp(700.0) && fli_unlogged(0.0) == 1.0 && fli_unlogged(709.0) isa Float64 && fli_unlogged(710.0) === nothing
    emit("AC11", "(d) fli_unlogged on small values returns exp(F); above 709.78 returns nothing, never Inf", okf, "true", pf(okf))
end

# ---------------------------------------------------------------------------------------------------
# AC12: health and frames
section("AC12") do
    seed = vseed(12, 1)
    x0s(flagged) = [SVector{3, Float64}(0, 0, j in flagged ? 1 : 0) for j in 1:6]
    zero_map, col2_map = fl_Fault(:zero, 0), fl_Fault(:col2, 0)
    oka = true
    for chain in ((FLI(),), (FLI(replace = false), SALI()))
        ok = propagate(zero_map, x0s(()), 12; k = 2, tau = 3, seed = seed, reducers = chain, chunk = 2)
        res = nothing
        err = fl_health_error(() -> (res = propagate(zero_map, x0s((4,)), 12; k = 2, tau = 3, seed = seed, reducers = chain, chunk = 2)))
        res2 = nothing
        err2 = fl_health_error(() -> (res2 = propagate(zero_map, x0s((5, 2)), 12; k = 2, tau = 3, seed = seed, reducers = chain, chunk = 2)))
        eye = [SMatrix{2, 2, Float64}(1, 0, 0, 1) for _ in 1:6]
        err3 = fl_health_error(() -> propagate(col2_map, x0s((5,)), 12; k = 2, tau = 3, seed = seed, reducers = chain, chunk = 2, W0 = eye))
        oka &= length(ok.finals) == 6 && err isa NumericalHealthError && res === nothing && err.what == "deviation matrix" && err.orbit == 4 && err.step == 3 && err.backend === :cpu &&
               err.column == 1 && err.value == 0.0 && err2 isa NumericalHealthError && err2.orbit == 2 && res2 === nothing && err3 isa NumericalHealthError && err3.orbit == 5 &&
               err3.column == 2 && err3.step == 3
    end
    x = SVector{3, Double64}(0, 0, 1)
    err4 = fl_health_error(() -> propagate_orbit(zero_map, x, fl_frame(seed, 1, Double64), 12, 3, (FLI(),); orbit = 4))
    err5 = fl_health_error(() -> propagate_orbit(col2_map, x, SMatrix{2, 2, Double64}(1, 0, 0, 1), 12, 3, (FLI(),); orbit = 4))
    e1 = try
        propagate(fl_Boom(), [SVector(0.0, 0.0)], 4; k = 1, tau = 2, seed = seed, reducers = (FLI(),)); nothing
    catch ex
        ex
    end
    oka &= err4 isa NumericalHealthError && err4.orbit == 4 && err4.column == 1 && err5 isa NumericalHealthError && err5.column == 2 && e1 isa ArgumentError && occursin("FLI needs k = 2, got k = 1", e1.msg)
    emit("AC12", "(a) zero tangent on a flagged orbit (propagate chunk 2, chains {FLI} and {FLI(replace = false), SALI}; Double64 per-orbit): kind deviation matrix, orbit, step 3, backend cpu, column 1 (2 for diag(1,0)), value 0, nothing stored, lowest id first; k = 1 ArgumentError", oka, "true", pf(oka))
    okb = true
    for chain in ((FLI(),), (FLI(replace = false), SALI()))
        resn = nothing
        errn = fl_health_error(() -> (resn = propagate(fl_Fault(:nan, 3), x0s((3,)), 20; k = 2, tau = 5, seed = vseed(12, 2), reducers = chain, chunk = 2)))
        okb &= errn isa NumericalHealthError && resn === nothing && errn.what == "deviation matrix" && errn.orbit == 3 && errn.step == 5 && errn.backend === :cpu && errn.column == 0 && errn.threshold == Inf && !isfinite(errn.value)
    end
    errn2 = fl_health_error(() -> propagate_orbit(fl_Fault(:nan, 3), x, fl_frame(vseed(12, 2), 1, Double64), 20, 5, (FLI(),); orbit = 3))
    okb &= errn2 isa NumericalHealthError && errn2.column == 0 && errn2.step == 5
    emit("AC12", "(b) NaN in the tangent: kind deviation matrix, column 0, threshold Inf, orbit 3, step 5 (not masked), both chains and Double64 per-orbit", okb, "true", pf(okb))
    seed = vseed(12, 3)
    r = FLI()
    for T in FL_TYPES
        epsT = fl_eps(T)
        big1 = T(1) + 1000 * T(epsT)
        rejected_ok = true
        for (W, col) in ((SMatrix{2, 2, T}(2, 0, 0, 1), 1), (SMatrix{2, 2, T}(1, 0, 0, 3), 2), (SMatrix{2, 2, T}(1, 0, 0, 0), 2), (SMatrix{2, 2, T}(0, 0, 0, 1), 1),
                         (SMatrix{2, 2, T}(big1, 0, 0, 1), 1), (SMatrix{2, 2, T}(1, 0, 0, -big1), 2))
            e = fl_health_error(() -> check_frame(r, W; orbit = 7))
            rejected_ok &= e isa NumericalHealthError && e.what == "initial frame" && e.orbit == 7 && e.step == 0 && e.backend === :cpu && e.column == col && 0 < e.threshold < 1e-6 && e.value > e.threshold
        end
        maxdef = 0.0
        accepted_ok = true
        for W in (SMatrix{2, 2, T}(1, 0, 0, 1), SMatrix{2, 2, T}(1, 0, 1, 0), SMatrix{2, 2, T}(0, 1, 0, -1), SMatrix{2, 2, T}(1 + T(epsT), 0, 0, 1))
            delta = nothing
            accepted_ok &= fl_health_error(() -> (delta = check_frame(r, W; orbit = 3))) === nothing && delta !== nothing && 0 <= delta <= 2 * epsT
            delta === nothing || (maxdef = max(maxdef, Float64(delta)))
        end
        argok = try; check_frame(r, SMatrix{2, 3, T}(1, 0, 0, 1, 0, 0)); false; catch ex; ex isa ArgumentError; end
        nid = T === Float64 ? 10_000 : 1_000
        defects = [Float64(check_frame(r, fl_frame(seed, j, T); orbit = j)) for j in 1:nid]
        emit("AC12", "(c) $T: columns of norm 2, 3, 0, 1 + 1000 eps rejected (kind initial frame, step 0, orbit, column); unit, dependent and 1 + eps columns accepted (largest defect $(fmt(maxdef / epsT)) eps; threshold 8 eps); 2-column ArgumentError", rejected_ok && accepted_ok && argok, "true", pf(rejected_ok && accepted_ok && argok))
        emit("AC12", "(c) $T: seeded frame accepted for all ids 1..$nid; largest norm defect in units of eps (threshold 8)", maximum(defects) / epsT, "<= 8, 0 rejected", pf(maximum(defects) <= 8 * epsT))
    end
    W0 = [fl_frame(vseed(12, 4), 1, Float64), SMatrix{2, 2, Float64}(2, 0, 0, 1)]
    erb = fl_health_error(() -> begin
        foreach(j -> check_frame(FLI(), W0[j]; orbit = j), eachindex(W0))
        propagate_orbit(fl_Boom(), SVector(0.0, 0.0), W0[2], 4, 2, (FLI(),))
    end)
    boom = try
        propagate_orbit(fl_Boom(), SVector(0.0, 0.0), W0[1], 4, 2, (FLI(),)); false
    catch ex
        ex isa ErrorException
    end
    okr = erb isa NumericalHealthError && erb.what == "initial frame" && erb.orbit == 2 && erb.step == 0 && boom
    emit("AC12", "(c) caller-side frame check rejects the non-unit frame before any map call (fl_Boom throws if called); an accepted frame reaches the map (D67 form)", okr, "true", pf(okr))
    mu = exp(3.0)
    sys = LinearSymplecticMap(SMatrix{2, 2, Float64}(mu, 0, 0, 1 / mu))
    Wd = fl_frame(vseed(12, 5), 1, Float64)
    okd = true
    fb200 = Float64[]
    for chain in ((FLI(replace = false), SALI()), (FLI(),))
        e200 = fl_health_error(() -> propagate_orbit(sys, SVector(0.0, 0.0), Wd, 200, 200, chain))
        _, _, st = propagate_orbit(sys, SVector(0.0, 0.0), Wd, 200, 200, chain)
        f = reducer_final(FLI(), st[1])
        e300 = fl_health_error(() -> propagate_orbit(sys, SVector(0.0, 0.0), Wd, 300, 300, chain))
        push!(fb200, f.fli)
        okd &= e200 === nothing && f.n == 200 && 599 < f.fli < 601 && e300 isa NumericalHealthError && e300.what == "deviation matrix" && e300.step == 300 && e300.backend === :cpu
    end
    emit("AC12", "(d) A = diag(e^3, e^-3), Float64, both chains: tau = 200 no error at n = 200 (FLI_B = $(fmt(fb200[1])), $(fmt(fb200[2]))); tau = 300 error at handoff 300", okd, "true", pf(okd))
    res = propagate(StandardMap(2.0), [SVector(1.0, 2.0)], 8; k = 2, tau = 2, seed = vseed(12, 6), reducers = (FLI(),))
    f = res.finals[1][1]
    oke = reducer_provenance(FLI()).indicator == "fli" && reducer_provenance(FLI(replace = false)).indicator == "fli" && res.provenance.indicator == "fli" && res.provenance.k == 2 && res.provenance.tau == 2 &&
          res.provenance.precision == "Float64" && f.n == 8 && f.sum isa SVector{2, Float64} && f.comp isa SVector{2, Float64} && f.fli isa Float64 && f.fli1 isa Float64 && f.stepB isa Integer && f.colB isa Integer && f.step1 isa Integer
    emit("AC12", "(e) run record on {FLI} alone: indicator fli, k, tau, precision; final field types (the chain position and frame source are recorded by Experiments, D72 item 1)", oke, "true", pf(oke))
end

# ---------------------------------------------------------------------------------------------------
# AC13: chain, orbit and execution invariance
section("AC13") do
    seed = vseed(13)
    T = Float64
    bad = ncmp = 0
    bad_b = 0
    for tau in (1, 4)
        sys = StandardMap(fl_K(T, :two_pi))
        for (j, s) in enumerate(fl_starts(StableRNG(seed + 1), 5))
            x0 = fl_cast(T, s)
            W2 = fl_frame(seed, j, T)
            W1 = SMatrix{2, 1, T}(W2[:, 1])
            xa, _, sta = propagate_orbit(sys, x0, W1, 400, tau, (fl_Trace(),); orbit = j)
            xb, _, stb = propagate_orbit(sys, x0, W1, 400, tau, (Benettin(), fl_Trace()); orbit = j)
            xc, Wc, stc = propagate_orbit(sys, x0, W2, 400, tau, (FLI(), fl_Trace()); orbit = j)
            xd, Wd, std = propagate_orbit(sys, x0, W2, 400, tau, (FLI(replace = false), SALI(), fl_Trace()); orbit = j)
            xe, _, ste = propagate_orbit(sys, x0, W2, 400, tau, (fl_Trace(),); orbit = j)
            xf, Wf, stf = propagate_orbit(sys, x0, W2, 400, tau, (SALI(), fl_Trace()); orbit = j)
            tra, trb = reducer_final(fl_Trace(), sta[1]).xs, reducer_final(fl_Trace(), stb[2]).xs
            trc, trd = reducer_final(fl_Trace(), stc[2]).xs, reducer_final(fl_Trace(), std[3]).xs
            tre, trf = reducer_final(fl_Trace(), ste[1]).xs, reducer_final(fl_Trace(), stf[2]).xs
            bad += !(isequal(tra, trb) && isequal(tra, trc) && isequal(tra, trd) && isequal(tra, tre) && isequal(tra, trf) && xa === xb && xa === xc && xa === xd && xa === xe && xa === xf)
            fa, fd1 = reducer_final(FLI(), stc[1]), reducer_final(FLI(replace = false), std[1])
            bad_b += !(isequal(fa, fd1) && isequal(reducer_final(SALI(), std[2]), reducer_final(SALI(), stf[1])) && isequal(Wc, Wf) && isequal(Wd, Wf))
            ncmp += 1
        end
    end
    emit("AC13", "(a) orbit points at every handoff differ across the chains {}, {L1} (k=1), {FLI}, {FLI(replace=false), SALI}, {SALI}, {} (k=2); cases = (5 starts x tau in {1,4})", bad, "0 of $ncmp", pf(bad == 0))
    emit("AC13", "(b) in {FLI(replace=false), SALI}: FLI outputs differ from the stand-alone {FLI}, or SALI outputs from the stand-alone {SALI}, or final W differs among the chains with a replacer", bad_b, "0 of $ncmp", pf(bad_b == 0))
    # (b) propagate equals the serial per-orbit body
    for (reducers, taus) in (((FLI(replace = false), SALI()), (1, 4)), ((FLI(),), (4,))), tau in taus
        sys = StandardMap(fl_K(T, :two_pi))
        m, N = 4096, 1000
        x0 = [fl_cast(T, s) for s in fl_starts(StableRNG(seed + 2), m)]
        W0 = [fl_frame(seed, j, T) for j in 1:m]
        serial = [propagate_orbit(sys, x0[j], W0[j], N, tau, reducers; orbit = j) for j in 1:m]
        xs, Ws = [s[1] for s in serial], [s[2] for s in serial]
        fs = [map(reducer_final, reducers, s[3]) for s in serial]
        for chunk in (64, 1000)
            res = propagate(sys, x0, N; k = 2, tau = tau, seed = seed, reducers = reducers, chunk = chunk, W0 = W0)
            okb = isequal(res.x, xs) && isequal(res.W, Ws) && isequal(res.finals, fs)
            emit("AC13", "(b) propagate equals the serial map of propagate_orbit bitwise (x, W, finals), chain {$(join(map(r -> r isa FLI ? (r.replace ? "FLI" : "FLI(replace=false)") : "SALI", reducers), ", "))}, 4096 orbits, K = 2 pi, n = 1e3, tau = $tau, chunk = $chunk", okb, "true", pf(okb))
        end
    end
    # (c) 1 vs N threads in child processes
    hashes = Dict{Int, Dict{String, String}}(Threads.nthreads() => thread_hashes())
    for nt in (1, 2, 4)
        haskey(hashes, nt) && continue
        cmd = addenv(`$(Base.julia_cmd()) -t $nt --project=$(dirname(Base.active_project())) $(@__FILE__)`, "FL_CHILD" => "1")
        txt = read(pipeline(cmd; stderr = devnull), String)
        d = Dict{String, String}(m[1] => m[2] for m in eachmatch(r"HASH (.+) = (\d+)\n", txt))
        mt = match(r"THREADS (\d+)", txt)
        hashes[nt] = (mt !== nothing && parse(Int, mt.captures[1]) == nt) ? d : Dict("child failed" => first(txt, 100))
    end
    keys_ = sort(collect(keys(hashes[Threads.nthreads()])))
    ok_all = true
    for k in keys_
        hs = [get(hashes[nt], k, "missing") for nt in sort(collect(keys(hashes)))]
        ok = length(unique(hs)) == 1
        ok_all &= ok
        ok || emit("AC13", "(c) hashes differ for $k across threads $(sort(collect(keys(hashes)))): $hs", 0, "equal", "FAIL")
    end
    # all chunk sizes also agree with each other
    for cname in ("FLI", "FLI,SALI"), tau in (1, 4)
        hs = unique([hashes[Threads.nthreads()]["$cname tau=$tau chunk=$chunk"] for chunk in (64, 1000)])
        ok_all &= length(hs) == 1
    end
    emit("AC13", "(c) x, W and every final (FLI, FLI_1, sums, compensations, maximisers, SALI) bitwise identical for thread counts $(sort(collect(keys(hashes)))) (hashes in separate processes) and chunk sizes 64 and 1000: $(length(keys_)) cases (chains {FLI}, {FLI(false), SALI} x tau in {1,4} x chunk), 4096 orbits, K = 2 pi, n = 1e3; example hash $(first(values(hashes[1])))", ok_all, "true", pf(ok_all))
    for T in FL_TYPES, tau in (1, 4), reducers in ((FLI(),), (FLI(replace = false), SALI()))
        a = fl_alloc_bytes(StandardMap(fl_K(T, :two_pi)), fl_cast(T, SVector(1.0, 2.0)), fl_frame(seed, 1, T), 10_000, tau, reducers)
        emit("AC13", "(d) bytes allocated by propagate_orbit after warm-up, n = 1e4, tau = $tau, chain length $(length(reducers)), $T", a, "0", pf(a == 0))
    end
    # (e) column exchange and sign change
    badx = badn = 0
    ncase = 0
    for T in FL_TYPES, label in (:two_pi, :half), tau in (1, 4)
        c = vcase(label, T, 5)
        W0 = c.W0
        swap(W) = SMatrix{2, 2, T}(W[1, 2], W[2, 2], W[1, 1], W[2, 1])
        neg1(W) = SMatrix{2, 2, T}(-W[1, 1], -W[2, 1], W[1, 2], W[2, 2])
        neg2(W) = SMatrix{2, 2, T}(W[1, 1], W[2, 1], -W[1, 2], -W[2, 2])
        base = fl_run((FLI(),), c.sys, c.x0, W0, 400, tau; all = true)
        sw = fl_run((FLI(),), c.sys, c.x0, swap.(W0), 400, tau; all = true)
        n1 = fl_run((FLI(),), c.sys, c.x0, neg1.(W0), 400, tau; all = true)
        n2 = fl_run((FLI(),), c.sys, c.x0, neg2.(W0), 400, tau; all = true)
        for j in 1:5
            a, b = base[j][1].snaps, sw[j][1].snaps
            ncase += 1
            badx += !(isequal([s.fli for s in a], [s.fli for s in b]) && isequal([s.sum[1] + s.comp[1] for s in a], [s.sum[2] + s.comp[2] for s in b]) &&
                      isequal([s.fli1 for s in b], accumulate(max, [s.sum[2] + s.comp[2] for s in a]; init = zero(T))))
            for other in (n1, n2)
                o = other[j][1].snaps
                badn += !(isequal([s.fli for s in a], [s.fli for s in o]) && isequal([s.fli1 for s in a], [s.fli1 for s in o]) && isequal([s.sum for s in a], [s.sum for s in o]) && isequal([s.comp for s in a], [s.comp for s in o]))
            end
        end
    end
    emit("AC13", "(e) exchange of the two columns: FLI_B and the exchanged sums not bitwise equal, or FLI_1 != running max of the other sum: $badx of $ncase; sign change of either column changing anything: $badn of $(2ncase) (K = 2 pi, 0.5; Float64, Double64; tau in {1,4}; 5 starts; 400 steps, every handoff)", badx + badn, "0", pf(badx + badn == 0))
end

# ---------------------------------------------------------------------------------------------------
# AC14: wrapped versus lifted
function v_rho_M(Kf::Float64, x0, w0, n::Int)
    return setprecision(BigFloat, FL_BITS) do
        K = BigFloat(Kf)
        x, p = fl_big(x0[1]), fl_big(x0[2])
        a, b, c_, d = one(K), zero(K), zero(K), one(K)
        M = max(1.0, abs(Kf))
        for _ in 1:n
            M = max(M, Float64(abs(x)), Float64(abs(p)))
            c = K * cos(x)
            a, b, c_, d = (1 + c) * a + c_, (1 + c) * b + d, c * a + c_, c * b + d
            p = p + K * sin(x)
            x = x + p
        end
        w1, w2 = fl_big(w0[1]), fl_big(w0[2])
        y1, y2 = a * w1 + b * w2, c_ * w1 + d * w2
        F = a^2 + b^2 + c_^2 + d^2
        dt = a * d - b * c_
        smax = sqrt((F + sqrt(F^2 - 4 * dt^2)) / 2)
        M = max(M, Float64(abs(x)), Float64(abs(p)))
        (Float64(smax / sqrt(y1^2 + y2^2)), M)
    end
end

section("AC14") do
    seed = vseed(14)
    rng = StableRNG(seed)
    dy() = rand(rng, 0:2048) / 1024
    T = Float64
    nbad = 0
    for tau in (1, 3)
        zs = [SVector{2, T}(T(dy()), T(dy())) for _ in 1:50]
        W0 = [fl_frame(seed, j, T) for j in 1:50]
        rl = fl_run((FLI(),), StandardMap(exact(T, 0)), zs, W0, 300, tau; all = true)
        rw = fl_run((FLI(),), WrappedStandardMap(exact(T, 0)), torus_state.(zs), W0, 300, tau; all = true)
        nbad += count(j -> !isequal(rl[j][1].snaps, rw[j][1].snaps), 1:50)
    end
    emit("AC14", "(a) K = 0, dyadic starts: starts (of 2 x 50, tau in {1,3}, all handoffs to 300) whose FLI, sums, compensations, maximisers differ between wrapped and lifted", nbad, "0", pf(nbad == 0))
    m = 50
    c = vcase(:half, T, m)
    for n in (10, 100)
        rl = fl_run((FLI(),), c.sys, c.x0, c.W0, n, 1; steps = [n])
        rw = fl_run((FLI(),), WrappedStandardMap(0.5), torus_state.(c.x0), c.W0, n, 1; steps = [n])
        excluded = 0
        w = 0.0
        diffmax = 0.0
        for j in 1:m
            a, b = rl[j][1].snaps[1], rw[j][1].snaps[1]
            rho(i, step) = v_rho_M(0.5, c.x0[j], c.W0[j][:, i], step)
            bound(step) = step == 0 ? (0.0, 0.0) : (r1 = rho(1, step); r2 = rho(2, step); (64 * 0.5 * r1[2] * step^2 * eps(Float64) * (r1[1] + r2[1]), max(r1[1], r2[1])))
            pairs = ((a.stepB, b.stepB), (a.step1, b.step1))
            forms = ((a.fli, b.fli), (a.fli1, b.fli1))
            for (k, (sa, sb)) in enumerate(pairs)
                ba, bb = bound(sa), bound(sb)
                if max(ba[2], bb[2]) > 50
                    excluded += 1
                    continue
                end
                d = abs(forms[k][1] - forms[k][2])
                diffmax = max(diffmax, d)
                w = max(w, rat(d, max(ba[1], bb[1])))
            end
        end
        ratio_line("AC14", "(b) K = 0.5 regular class, n = $n, 50 starts, both forms: max |FLI_w - FLI_l| / (64 K M_state j^2 eps (rho1 + rho2)) at the larger maximiser bound; max |difference| $(fmt(diffmax)); excluded (max rho > 50) $excluded of 100", w)
        emit("AC14", "(b) n = $n: excluded comparisons (of 100)", excluded, "<= 10", pf(excluded <= 10))
    end
end

# ---------------------------------------------------------------------------------------------------
# AC15: ChaosTools GALI_2 oracle (one-sided line)
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

# B_oracle of the sali criterion 16 (A15 form with eps = eps(Float64)) from the sali helpers, at the steps `cps`, for start j of a case.
function v_boracle(c, j, N, cps)
    pts = fl_pts(c.sys, c.x0[j], N)
    ref = sl_ref(sl_mats_std(c.K, pts), c.W0[j], N; keep = cps)
    Ms = lb_tangents(c.sys, c.x0[j], N)
    ps = sl_psis(Ms, c.W0[j], cps)
    return [(saliref = Float64(sl_at(ref, n).S), B = sl_BSALI(EPS64, sl_Lambda(ref.lnsali, n), Float64(sl_at(ref, n).S), ps[i][1], ps[i][2])) for (i, n) in enumerate(cps)]
end

section("AC15") do
    meas = 4 * EPS64
    DATA[:ct] = NamedTuple[]
    for (label, m, n) in ((:half, 20, 10_000), (:two_pi, 100, 12), (:five, 100, 12))
        starts = label === :half ? vcase(:half, Float64, 20).x0 : vcase(label, Float64, 1000).x0[1:m]
        K = Float64(fl_K(Float64, label))
        tot = mism = 0
        maxd = 0.0
        for z0 in starts
            mm, dd = ct_orbit_check(K, z0, n)
            mism += mm; maxd = max(maxd, dd); tot += n
        end
        emit("AC15", "(c) K = $K: oracle orbit vs forward(StandardMap) over $(length(starts)) starts x $n steps: steps with a bitwise difference; largest coordinate difference $(fmt(maxd))", mism, "0 of $tot (else add the orbit divergence term)", pf(mism == 0))
    end
    # (a) K = 0.5 regular class, 20 starts, n in {10, 100, 1e3, 1e4}
    c = vcase(:half, Float64, 20)
    cps = [10, 100, 1000, 10_000]
    rows = vrows(:half, Float64, 1, 20, 10_000)
    gs = [ct_gali(0.5, c.x0[j], c.W0[j], 10_000) for j in 1:20]
    bor = [v_boracle(c, j, 10_000, cps) for j in 1:20]
    w = -Inf
    minmargin = Inf
    npairs = nexc = 0
    for r in rows
        r.n in cps || continue
        go = gs[r.j][r.n + 1]
        bo = bor[r.j][findfirst(==(r.n), cps)].B
        lower = -log(go) / 2
        w = max(w, (lower - r.fli) / (r.BB + bo / (2 * go)))
        minmargin = min(minmargin, Float64(r.fli) - lower)
        nexc += r.fli < lower
        npairs += 1
        push!(DATA[:ct], (K = 0.5, n = r.n, fli = Float64(r.fli), lower = lower, bound = r.BB + bo / (2 * go), resolved = true))
    end
    emit("AC15", "(a) K = 0.5 regular, 20 starts, n in {10,100,1e3,1e4}, $npairs pairs: max (-ln g_oracle/2 - FLI_B) / (B_FLI + B_oracle/(2 g_oracle)) (one-sided; negative means FLI_B above the oracle value); min (FLI_B - (-ln g_oracle/2)) = $(fmt(minmargin)); pairs with FLI_B below the oracle value: $nexc", w, "<= 1", pf(w <= 1))
    for label in (:two_pi, :five)
        K = Float64(fl_K(Float64, label))
        c = vcase(label, Float64, 1000)
        cps = [2, 5, 8, 12]
        rows = vrows(label, Float64, 1, 1000, 12; cps = cps)
        tot = excl = nexc = 0
        w = -Inf
        minmargin = Inf
        gexc = NaN
        for j in 1:1000
            bo = v_boracle(c, j, 12, cps)
            gj = ct_gali(K, c.x0[j], c.W0[j], 12)
            for r in (r for r in rows if r.j == j)
                i = findfirst(==(r.n), cps)
                go = gj[r.n + 1]
                tot += 1
                if !(bo[i].saliref >= 10 * meas && go >= 10 * meas)
                    excl += 1
                    continue
                end
                lower = -log(go) / 2
                w = max(w, (lower - r.fli) / (r.BB + bo[i].B / (2 * go)))
                Float64(r.fli) - lower < minmargin && (gexc = go)
                minmargin = min(minmargin, Float64(r.fli) - lower)
                nexc += r.fli < lower
                push!(DATA[:ct], (K = K, n = r.n, fli = Float64(r.fli), lower = lower, bound = r.BB + bo[i].B / (2 * go), resolved = true))
            end
        end
        emit("AC15", "(b) K = $(fmt(K)), 1000 starts, resolved handoffs n in {2,5,8,12}, $(tot - excl) pairs: max (-ln g_oracle/2 - FLI_B) / (B_FLI + B_oracle/(2 g_oracle)); min (FLI_B - (-ln g_oracle/2)) = $(fmt(minmargin)) (at g_oracle = $(fmt(gexc))); pairs with FLI_B below the oracle value: $nexc", w, "<= 1", pf(w <= 1))
        emit("AC15", "(b) K = $(fmt(K)): pairs excluded as unresolved ($excl of $tot)", excl / tot, "<= 0.10", pf(excl / tot <= 0.10))
    end
end

# ---------------------------------------------------------------------------------------------------
# Figures
section("figure") do
    seed = vseed(20)
    fig = CM.Figure(size = (1700, 520))
    # (1) regular orbits: FLI_B against ln n (linear law, slope 1)
    ax1 = CM.Axis(fig[1, 1]; xlabel = "ln n", ylabel = "FLI_B (log form)", title = "regular orbits, K = 0.5: FLI_B = ln n + O(1), slope 1 (Q8)")
    c = vcase(:half, Float64, 20)
    nn = unique(round.(Int, 10 .^ range(1, 6; length = 60)))
    run = fl_run((FLI(),), c.sys, c.x0[1:8], c.W0[1:8], 1_000_000, 1; steps = nn)
    for r in run
        CM.lines!(ax1, log.(Float64.(nn)), [Float64(s.fli) for s in r[1].snaps]; linewidth = 1)
    end
    CM.lines!(ax1, [log(10.0), log(1e6)], [log(10.0), log(1e6)]; color = :black, linestyle = :dash, label = "ln n (slope 1)")
    CM.axislegend(ax1; position = :lt, labelsize = 10)
    # (2) chaotic orbits: FLI_B against n with the law lambda_1 n (lambda_1 = ensemble mean of X(1e4) of the L1 runs on the identical orbits)
    ax2 = CM.Axis(fig[1, 2]; xlabel = "n", ylabel = "ensemble mean of FLI_B(n) (200 starts)", title = "chaotic orbits: FLI_B = lambda_1 n + O(1) (Q8)")
    for (K, col) in ((:five, :seagreen), (:two_pi, :darkorange), (:ten, :purple))
        sd = vseed(21, findfirst(==(K), (:five, :two_pi, :ten)))
        sys = StandardMap(fl_K(Float64, K))
        x0 = [fl_cast(Float64, s) for s in fl_starts(StableRNG(sd), 200)]
        W0 = [fl_frame(sd, j, Float64) for j in 1:200]
        nl = collect(1:5:300)
        r = fl_run((FLI(),), sys, x0, W0, 300, 1; steps = nl)
        l1 = fl_run((Benettin(),), sys, x0, [SMatrix{2, 1, Float64}(W[:, 1]) for W in W0], 10_000, 1; steps = [10_000])
        lam = mean(q[1].snaps[1].lambda[1] for q in l1)
        CM.lines!(ax2, Float64.(nl), [mean(Float64(q[1].snaps[i].fli) for q in r) for i in eachindex(nl)]; color = col, label = "K = $(Float64(fl_K(Float64, K))): mean FLI_B")
        CM.lines!(ax2, Float64.(nl), lam .* nl; color = col, linestyle = :dash, label = "lambda_1 n, lambda_1 = $(fmt(lam))")
    end
    CM.axislegend(ax2; position = :lt, labelsize = 9)
    # (3) bitwise agreement with the L1 sum: FLI_1 against the running maximum of n X_L1(n), regular and chaotic orbit
    ax3 = CM.Axis(fig[1, 3]; xlabel = "running maximum of n X_L1(n) (Benettin, k = 1)", ylabel = "FLI_1 (FLI reducer)", title = "FLI_1 = running maximum of the L1 sum (criterion 8)")
    worst = 0.0
    for (K, col) in ((:half, :seagreen), (:two_pi, :darkorange))
        cc = vcase(K, Float64, 20)
        run = fl_run((FLI(),), cc.sys, cc.x0[1:3], cc.W0[1:3], 1000, 1; all = true)
        l1 = fl_run((Benettin(),), cc.sys, cc.x0[1:3], [SMatrix{2, 1, Float64}(W[:, 1]) for W in cc.W0[1:3]], 1000, 1; all = true)
        for j in 1:3
            rm = accumulate(max, [i * s.lambda[1] for (i, s) in enumerate(l1[j][1].snaps)]; init = 0.0)
            f1 = [Float64(s.fli1) for s in run[j][1].snaps]
            worst = max(worst, maximum(abs.(f1 .- rm) ./ (4 * EPS64 .* max.(1, f1))))
            CM.scatter!(ax3, rm, f1; color = (col, 0.5), markersize = 4, label = j == 1 ? "K = $(Float64(fl_K(Float64, K)))" : nothing)
        end
    end
    CM.lines!(ax3, [0, 1000], [0, 1000]; color = :black, linestyle = :dash)
    CM.axislegend(ax3; position = :lt, labelsize = 10)
    ax3.title = "FLI_1 against max of n X_L1; worst |diff|/(4 eps max(1,FLI_1)) = $(fmt(worst))"
    CM.save(joinpath(OUT, "fli-growth-laws.png"), fig)
    say("# figure saved: fli-growth-laws.png (panel 3 worst ratio $(fmt(worst)))")
    haskey(DATA, :ct) || return
    fig2 = CM.Figure(size = (1200, 480))
    pts = DATA[:ct]
    ax4 = CM.Axis(fig2[1, 1]; xlabel = "-ln(g_oracle)/2 (ChaosTools gali, k = 2)", ylabel = "FLI_B (Chaoslab)", title = "AC15: FLI_B above the oracle's lower value")
    for (K, col) in ((0.5, :seagreen), (Float64(fl_K(Float64, :two_pi)), :darkorange), (5.0, :purple))
        sel = [p for p in pts if p.K == K]
        CM.scatter!(ax4, [p.lower for p in sel], [p.fli for p in sel]; color = (col, 0.5), markersize = 4, label = "K = $(fmt(K))")
    end
    mx = maximum(p.fli for p in pts)
    CM.lines!(ax4, [0, mx], [0, mx]; color = :black)
    CM.axislegend(ax4; position = :lt, labelsize = 10)
    ax5 = CM.Axis(fig2[1, 2]; xlabel = "-ln(g_oracle)/2", ylabel = "(lower - FLI_B) and the bound B_FLI + B_oracle/(2 g)", yscale = log10, title = "AC15: margin FLI_B - lower (points) against the bound (crosses)")
    for (K, col) in ((0.5, :seagreen), (Float64(fl_K(Float64, :two_pi)), :darkorange), (5.0, :purple))
        sel = [p for p in pts if p.K == K]
        CM.scatter!(ax5, [p.lower for p in sel], [max(p.fli - p.lower, 1e-18) for p in sel]; color = (col, 0.5), markersize = 4, label = "K = $(fmt(K)): FLI_B - lower (margin)")
        CM.scatter!(ax5, [p.lower for p in sel], [p.bound for p in sel]; color = col, marker = :xcross, markersize = 5, label = "K = $(fmt(K)): bound")
    end
    CM.axislegend(ax5; position = :lt, labelsize = 8)
    CM.save(joinpath(OUT, "fli-chaostools.png"), fig2)
    say("# figure saved: fli-chaostools.png")
end

say(@sprintf("# verdicts: %d PASS, %d FAIL, %d other", count(==("PASS"), VERDICTS), count(==("FAIL"), VERDICTS), count(v -> v != "PASS" && v != "FAIL", VERDICTS)))
say(@sprintf("# total wall time %.1f s", time() - T0))
close(LOG)
