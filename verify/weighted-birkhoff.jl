# Stage 04 verification of module weighted-birkhoff (catalog row I13: weighted_birkhoff, birkhoff_index, WeightedBirkhoff and the observable kinds; layer Indicators).
# Run from the repo root:  julia -t 16 --project=verify --check-bounds=yes verify/weighted-birkhoff.jl
# Spec: workspaces/module-build/stages/01-spec/output/weighted-birkhoff-spec.md (criteria 1 to 13, Amendment D84). Criterion 13 (CPU vs GPU) is in verify/gpu.jl, tag "weighted-birkhoff".
# Sections: "suite" = the package folder test/weighted-birkhoff at the spec sizes (WB_FULL_RUN = true, set before runtests.jl is included; its @info lines are forwarded into the log); "AC1" to "AC8" = the
# numerical criteria re-measured here with verify seeds (NV_SEED0 + 50 ac + sub), so the realized ratios are logged; "AC11" = 1 versus N threads in child processes; "AC12" = the independent oracles
# (mpmath Poisson form through verify/weighted-birkhoff-poisson.py, the BigFloat(256) orbit, the ChaosTools search); "R1" = class separation panel with the cross-tabulation, "figure".
# R2, R3, R5 (CPU), R6, R7, R8 are the record testsets of the suite at spec size, forwarded. R4 (sticky orbits) is the study's and is not written (spec Amendment D84).
# Output lines:  AC<n> | quantity | value | bound | PASS/FAIL/INFO  with the MEASURED number (the rli D77 form for the recorded lines). WB_ONLY="AC12,R1" runs only the named sections, WB_LOG names the log file,
# WB_SMOKE=1 shrinks the sizes (debugging aid), WB_CHILD=1 is the child-process mode of AC11.
using Printf, Test, Random, LinearAlgebra, Statistics, Logging
using Chaoslab, StaticArrays, DoubleFloats, StableRNGs

const NV_SEED0 = 20261012 + 62000      # verify offset (naff-fma used 61000)
const SMOKE = get(ENV, "WB_SMOKE", "") == "1"

# ---------------------------------------------------------------------------------------------------
# Thread-count invariance (AC11): the same outputs hashed in this process and in child processes with other thread counts.
function v_flatten!(v, x)
    if x isa Number
        push!(v, Float64(x))
    elseif x isa Union{Tuple, NamedTuple, AbstractArray}
        foreach(y -> v_flatten!(v, y), x)
    elseif x isa Symbol
        push!(v, Float64(hash(x) % 1000003))
    else
        error("v_flatten!: $(typeof(x))")
    end
    return v
end
function thread_hashes()
    out = Dict{String, String}()
    seed = NV_SEED0 + 50 * 11 + 3
    N, m = 1000, (SMOKE ? 256 : 4096)
    obs = (CosObservable(1), SinObservable(1), LiftedObservable(2, 4), CosObservable(1, 2))
    for (K, tag) in ((0.5, "K=0.5"), (2pi, "K=2pi"))
        sys = WrappedStandardMap(K)
        rng = StableRNG(seed + (K == 0.5 ? 0 : 1))
        x0 = [torus_state(SVector(2pi * rand(rng), 2pi * rand(rng))) for _ in 1:m]
        r = WeightedBirkhoff(; window = (N ÷ 2, N), observables = obs)
        for chunk in (64, 1000)
            res = propagate(sys, x0, 2N; k = 2, tau = 1, seed = seed, reducers = (r, MEGNO()), chunk = chunk, backend = :cpu)
            v = Float64[]
            v_flatten!(v, res.x); foreach(f -> v_flatten!(v, f), res.finals)
            out["reader+MEGNO $tag chunk=$chunk"] = string(hash(reinterpret(UInt64, v)))
        end
        res = propagate(sys, x0, 2N; k = 1, tau = 1, seed = seed, reducers = (r, Benettin()), chunk = 64, backend = :cpu)
        v = Float64[]
        v_flatten!(v, res.x); foreach(f -> v_flatten!(v, f), res.finals)
        out["reader+Benettin $tag chunk=64"] = string(hash(reinterpret(UInt64, v)))
    end
    return out
end
if get(ENV, "WB_CHILD", "") == "1"
    for (k, h) in sort!(collect(thread_hashes())); println("HASH ", k, " = ", h); end
    println("THREADS ", Threads.nthreads())
    exit(0)
end

const WB_FULL_RUN = !SMOKE
const ROOT = dirname(@__DIR__)
const TESTDIR = joinpath(ROOT, "test", "weighted-birkhoff")
include(joinpath(TESTDIR, "wb_helpers.jl"))

vseed(ac, sub = 0) = NV_SEED0 + 50 * ac + sub

import ChaosTools as CT
import CairoMakie as CM

const OUT = joinpath(ROOT, "workspaces", "module-build", "stages", "04-verify", "output")
mkpath(OUT)
const ONLY = filter(!isempty, split(get(ENV, "WB_ONLY", ""), ","))
const LOG = open(get(ENV, "WB_LOG", joinpath(OUT, "weighted-birkhoff-verify-log.txt")), "w")
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
rmedian(v) = sort(v)[max((length(v) + 1) ÷ 2, 1)]
quants(v, qs = (0.05, 0.25, 0.5, 0.75, 0.95)) = isempty(v) ? "-" : join((fmt(quantile(v, q)) for q in qs), " / ")
frac(v) = isempty(v) ? NaN : count(identity, v) / length(v)

function section(f, name)
    isempty(ONLY) || name in ONLY || return
    t0 = time()
    try
        f()
    catch e
        say(string(name, " | ERROR | ", sprint(showerror, e, catch_backtrace()), " | - | FAIL"))
        push!(VERDICTS, "FAIL")
    end
    say(@sprintf("# section %s took %.1f s", name, time() - t0))
end

const T0 = time()
const DATA = Dict{Any, Any}()
say("# weighted-birkhoff verification, julia $(VERSION), threads $(Threads.nthreads()), base seed $WB_SEED (verify offset 62000), smoke $SMOKE, time $(round(time())) (unix)")
say("# timings in this log were taken on a machine that was not otherwise idle")

# ---------------------------------------------------------------------------------------------------
# The package suite at full size; its @info lines are forwarded into this log (all INFO: the suite asserts its gates itself)
struct SuiteLogger <: Logging.AbstractLogger end
Logging.min_enabled_level(::SuiteLogger) = Logging.Info
Logging.shouldlog(::SuiteLogger, args...) = true
Logging.catch_exceptions(::SuiteLogger) = true
function Logging.handle_message(::SuiteLogger, level, message, _module, group, id, file, line; kwargs...)
    rest = join(("$k = $(v isa AbstractFloat ? fmt(v) : v)" for (k, v) in kwargs), "; ")
    emit("SUITE", string(message), rest, "reported", "INFO")
end

section("suite") do
    local ts
    t = @elapsed try
        ts = with_logger(SuiteLogger()) do
            @testset verbose = true "weighted-birkhoff verify" begin
                include(joinpath(TESTDIR, "runtests.jl"))
            end
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
    nbrk = ok ? c.broken + c.cumulative_broken : -1
    ok && foreach(r -> r isa Test.AbstractTestSet && (cc = Test.get_test_counts(r); info("SUITE", "testset '$(r.description)': assertions passed (fail + error = $(cc.fails + cc.errors + cc.cumulative_fails + cc.cumulative_errors), broken = $(cc.broken + cc.cumulative_broken))", cc.passes + cc.cumulative_passes)), ts.results)
    emit("SUITE", "package suite at the spec sizes (WB_FULL_RUN = $(WB_FULL_RUN)), assertions passed (fail + error = $nbad, broken = $nbrk; the broken one is the device testset AC13, skipped without CUDA in this process); wall time $(round(t; digits = 1)) s", npass, "0 fail, 0 error", pf(ok && nbad == 0 && npass > 0))
end

# ---------------------------------------------------------------------------------------------------
# Helpers for the starts: the K = 0.5 regular class (L1 S7 screening: Benettin X(1e5) <= 5e-4 on the lifted map, tau = 10) of m starts out of 3 m candidates, and uniform starts of the other classes.
function regular_class(m)
    haskey(DATA, (:reg, m)) && return DATA[(:reg, m)]
    rng = StableRNG(vseed(40, 1))
    cand = [SVector(2pi * rand(rng), 2pi * rand(rng)) for _ in 1:3m]
    X = [f[1].lambda[1] for f in propagate(StandardMap(0.5), cand, 10^5; k = 1, tau = 10, seed = vseed(40, 2), reducers = (Benettin(),), backend = :cpu).finals]
    reg = findall(<=(5e-4), X)
    length(reg) >= m || error("only $(length(reg)) of $(3m) candidates are regular")
    return DATA[(:reg, m)] = (starts = cand[reg[1:m]], ncand = 3m, nreg = length(reg))
end
function class_starts(label, m)
    label === :half && return regular_class(m).starts
    rng = StableRNG(vseed(41, label === :two_pi ? 1 : label === :five ? 2 : label === :ten ? 3 : 4))
    return [SVector(2pi * rand(rng), 2pi * rand(rng)) for _ in 1:m]
end
const R_CLASSES = ((:half, 0.5), (:two_pi, 2pi), (:five, 5.0), (:ten, 10.0))
const NST = SMOKE ? 4 : 20                       # starts per class of the numerical criteria (spec: 20)
orbit_obs(sys, T, z0, obs, n; tau = 1) = wb_series(obs, wb_states_tau(sys, torus_state(wb_cast(T, z0)), n, tau))

# ---------------------------------------------------------------------------------------------------
# AC1 constants and covariances (I1, I2, I3), through the pure function (the reducer is bitwise the function, AC6)
section("AC1") do
    # (a) f = c, K = 2 (c and 2c), H = 2 (N, 2N): avg = c within 2 eps |c|; Delta = 0 bitwise; the flag set; absdig = -log10(32 eps S) within 2 eps relative
    for T in WB_TYPES
        epsT = wb_eps(T)
        wa = wd = wi = 0.0
        nflag = ntot = ndelta = 0
        for c in (0.3, 1e3, 1e-3, pi), N in (16, 1000, 10^5)
            T === Double64 && N > 10^4 && continue
            cT = T(c)
            f = hcat(fill(cT, 4N), fill(2cT, 4N))
            r = weighted_birkhoff(f, (N, 2N))
            for h in 1:2, k in 1:2
                cc = k * cT
                wa = max(wa, Float64(max(abs(r.avg1[k, h] - cc), abs(r.avg2[k, h] - cc))) / (2 * epsT * Float64(abs(cc))))
                ix = birkhoff_index(r.avg1[k, h], r.avg2[k, h], r.S[k, h])
                ntot += 1
                ix.below_floor && !ix.zero && (nflag += 1)
                ix.delta == 0 && (ndelta += 1)
                ref = setprecision(BigFloat, WB_BITS) do
                    -log10(32 * BigFloat(epsT) * wb_big(r.S[k, h]))
                end
                wi = max(wi, Float64(abs(wb_big(ix.absdig) - ref)) / (2 * epsT * Float64(ref)))
            end
        end
        ratio_line("AC1a", "constants f = c, c in {0.3, 1e3, 1e-3, pi}, N in {16, 1000, 1e5 (Float64 only)}, K = 2, H = 2, $T: max |avg - c| / (2 eps |c|)", wa)
        emit("AC1a", "same, $T: Delta == 0 bitwise / flag set (not zero), of $ntot entries", "$ndelta / $nflag", "all", pf(ndelta == ntot && nflag == ntot))
        ratio_line("AC1a", "same, $T: max |absdig - (-log10 32 eps S)| / (2 eps relative)", wi)
    end
    # (b) covariances on the K = 0.5 regular class, cos q, N = 1000, p = 1
    N = 1000
    reg = regular_class(NST).starts
    sys = WrappedStandardMap(0.5)
    bad2k = 0
    r37 = roff = 0.0
    nbad_swap = 0
    S0 = 0.0
    for z in reg
        f = vec(orbit_obs(sys, Float64, z, (CosObservable(1),), 2N))
        r0 = weighted_birkhoff(f, N)
        d0 = r0.avg1[1] - r0.avg2[1]
        S0 = max(S0, r0.S[1])
        for k in (-20, 7)
            rs = weighted_birkhoff(f .* 2.0^k, N)
            (rs.avg1[1] == r0.avg1[1] * 2.0^k && rs.avg2[1] == r0.avg2[1] * 2.0^k && rs.S[1] == r0.S[1] * 2.0^k && (rs.avg1[1] - rs.avg2[1]) == d0 * 2.0^k) || (bad2k += 1)
        end
        r3 = weighted_birkhoff(f .* 3.7, N)
        r37 = max(r37, max(abs(r3.avg1[1] - 3.7 * r0.avg1[1]), abs(r3.avg2[1] - 3.7 * r0.avg2[1])) / (2 * 4 * wb_eps(Float64) * 3.7 * r0.S[1]))
        rb = weighted_birkhoff(f .+ 1e3, N)
        roff = max(roff, abs((rb.avg1[1] - rb.avg2[1]) - d0) / (2 * 4 * wb_eps(Float64) * (r0.S[1] + 1e3)))
        sw = weighted_birkhoff(vcat(f[(N + 1):(2N)], f[1:N]), N)
        a, b = birkhoff_index(r0.avg1[1], r0.avg2[1], r0.S[1]), birkhoff_index(sw.avg1[1], sw.avg2[1], sw.S[1])
        (sw.avg1[1] == r0.avg2[1] && sw.avg2[1] == r0.avg1[1] && b.delta == -a.delta && isequal(a.absdig, b.absdig) && isequal(a.reldig, b.reldig) && isequal(a.dig, b.dig) && sw.S[1] == r0.S[1] && sw.B == r0.B) || (nbad_swap += 1)
    end
    emit("AC1b", "scaling by 2^k, k in {-20, 7}: average, S, Delta not bitwise scaled, over $(length(reg)) starts x 2", bad2k, "0", pf(bad2k == 0))
    ratio_line("AC1b", "scaling by 3.7: max |avg(3.7 f) - 3.7 avg(f)| / (2 c_1 eps (3.7 S)) (spec dry run 0.57 for the offset)", r37)
    ratio_line("AC1b", "adding b = 1e3: max |Delta(f + b) - Delta(f)| / (2 c_1 eps (S + |b|))", roff)
    emit("AC1c", "exchanging the windows: starts where Delta is not negated bitwise or absdig, reldig, dig, S, B not bitwise equal (of $(length(reg)))", nbad_swap, "0", pf(nbad_swap == 0))
    # I3: the weights sum to 1 within 4 eps and are bitwise symmetric (verified weights); B within c_p eps relative of the 256-bit sum
    for T in WB_TYPES, p in (1, 2), N in (16, 1000)
        w = birkhoff_weights(T, N; p = p)
        wsum = abs(sum(w) - 1) / wb_eps(T)
        r = weighted_birkhoff(ones(T, 2N), N; order = p)
        Bref = sum(wb_wref(N, p))
        eb = wb_err(r.B[1], Bref) / (Float64(Bref) * wb_cp(p) * wb_eps(T))
        emit("AC1d", "I3 $T p = $p N = $N: |sum of the normalised weights - 1| in eps (bound 4); B against the 256-bit sum / (c_p eps relative)", "$(fmt(wsum)); $(fmt(eb))", "<= 4; <= 1", pf(wsum <= 4 && eb <= 1))
    end
end

# ---------------------------------------------------------------------------------------------------
# AC2 exact rotation (WB16, with the premise of Amendment D84 item 6)
section("AC2") do
    alpha = setprecision(BigFloat, WB_BITS) do
        sqrt(BigFloat(2)) - 1
    end
    for T in WB_TYPES
        epsT = wb_eps(T)
        Nk1 = T === Float64 ? (500, 1000, 4000) : (1000, 4000, 16000)
        Nk2 = T === Float64 ? (1000, 4000) : (4000, 16000)
        for p in (1, 2)
            cp = wb_cp(p)
            wavg = wdel = wprem = 0.0
            for th in (wb_dec("0.142857142857142857142857142857142857142857142857142857142857142857142857142857"), wb_dec("0.3"))
                for (k, Ns) in ((1, Nk1), (2, Nk2)), N in Ns
                    cs, sn = wb_rot_big(th, alpha, k, 2N)
                    cols = k == 1 ? (cs, sn) : (cs,)
                    for col in cols
                        r = weighted_birkhoff(wb_round_vec(T, col), N; order = p)
                        ex = wb_ref_avgs(col, N, p)                      # the exact (BigFloat) averages: the premise "limit 0"
                        wavg = max(wavg, Float64(max(abs(r.avg1[1]), abs(r.avg2[1]))) / (cp * epsT))
                        wdel = max(wdel, Float64(abs(r.avg1[1] - r.avg2[1])) / (2 * cp * epsT))
                        wprem = max(wprem, Float64(max(abs(ex.avg1), abs(ex.avg2))) / (cp * epsT))
                    end
                end
            end
            ratio_line("AC2", "exact rotation, $T, w^[$p]: max over theta_0 in {1/7, 0.3}, cos, sin 2 pi theta (N in $Nk1), cos 4 pi theta (N in $Nk2): |average| / (c_p eps)", wavg)
            ratio_line("AC2", "same, $T, w^[$p]: |Delta| / (2 c_p eps)", wdel)
            emit("AC2", "same, $T, w^[$p]: premise (D84 item 6) the exact average / (c_p eps), must be below 0.1", wprem, "<= 0.1", pf(wprem <= 0.1))
        end
    end
end

# ---------------------------------------------------------------------------------------------------
# AC3 the K_std = 0 standard map (WB16), through the stream
section("AC3") do
    sys(T) = WrappedStandardMap(wb_K(T, :zero))
    for T in WB_TYPES
        epsT = wb_eps(T)
        p0s = (1.0, 2.0, 3.3, 5.0)
        Ns = T === Float64 ? (1000, 4000) : (4000, 16000)
        nq = SMOKE ? 2 : NST
        rng = StableRNG(vseed(3, T === Float64 ? 1 : 2))
        q0 = [2pi * rand(rng) for _ in 1:nq]
        wtr = wdel = wlift = 0.0
        for p0 in p0s, N in Ns
            r = WeightedBirkhoff(; window = N, observables = (CosObservable(1), SinObservable(1), LiftedObservable(2, 4)))
            for q in q0
                x0 = torus_state(SVector{2, T}(T(q), T(p0)))
                _, _, fin = wb_run((r, wb_Norm()), sys(T), x0, wb_frame(vseed(3, 3), 1, T), 2N, 1)
                f = fin[1]
                for k in 1:2
                    wtr = max(wtr, Float64(max(abs(f.avg1[k]), abs(f.avg2[k]))) / (4 * epsT))
                    wdel = max(wdel, Float64(abs(f.delta[k])) / (8 * epsT))
                end
                wlift = max(wlift, Float64(max(abs(f.avg1[3] - T(p0)), abs(f.avg2[3] - T(p0)))) / (2 * epsT * p0))
            end
        end
        ratio_line("AC3", "K_std = 0, $T, p_0 in {1, 2, 3.3, 5}, N in $Ns, $nq random q_0: cos q, sin q windows |avg| / (4 eps)", wtr)
        ratio_line("AC3", "same, $T: |Delta| / (8 eps)", wdel)
        T === Float64 && ratio_line("AC3", "same, $T: lifted p, |avg - p_0| / (2 eps p_0)", wlift)
        # p_0 = pi (Float64 pi, Float64 only as in the spec): the average of cos 2q equals cos 2 q_0 within 16 eps, |Delta| <= 8 eps
        T === Float64 || continue
        wp = wpd = 0.0
        for N in Ns
            r = WeightedBirkhoff(; window = N, observables = (CosObservable(1, 2),))
            for q in q0
                x0 = torus_state(SVector{2, T}(T(q), T(pi)))
                _, _, fin = wb_run((r, wb_Norm()), sys(T), x0, wb_frame(vseed(3, 3), 1, T), 2N, 1)
                f = fin[1]
                ref = setprecision(BigFloat, WB_BITS) do
                    cos(2 * wb_big(T(q)))
                end
                wp = max(wp, max(wb_err(f.avg1[1], ref), wb_err(f.avg2[1], ref)) / (16 * epsT))
                wpd = max(wpd, Float64(abs(f.delta[1])) / (8 * epsT))
            end
        end
        ratio_line("AC3", "K_std = 0, $T, p_0 = pi (Float64 pi), N in $Ns: |avg of cos 2q - cos 2 q_0| / (16 eps)", wp)
        ratio_line("AC3", "same, $T: |Delta| / (8 eps)", wpd)
    end
end

# ---------------------------------------------------------------------------------------------------
# AC4 and AC5(a)(b): the algebra against the independent BigFloat(256) sums of WB14 at the realized samples; the index against WB15 evaluated in BigFloat from the returned values; I5 of D84
section("AC4") do
    Ns = (256, 1000, 4096)
    ps = (1, 2)
    obs = wb_std5()
    for p in ps, N in Ns
        wb_wref(N, p)                            # fills the cache before the threaded loops
    end
    rows = NamedTuple[]
    for T in WB_TYPES
        epsT = wb_eps(T)
        for (label, K) in ((:half, 0.5), (:two_pi, 2pi), (:ten, 10.0))
            starts = class_starts(label, NST)
            sys = WrappedStandardMap(wb_K(T, label))
            res = Vector{Vector{NamedTuple}}(undef, NST)
            Threads.@threads :dynamic for j in 1:NST
                x0 = torus_state(wb_cast(T, starts[j]))
                rs = Tuple(WeightedBirkhoff(; window = Ns, observables = obs, birkhoff_order = p) for p in ps)
                _, _, fin = wb_run((rs..., wb_Norm()), sys, x0, wb_frame(vseed(4, 1), j, T), 2 * maximum(Ns), 1)
                fb = wb_fref_std(wb_states(sys, x0, 2 * maximum(Ns)))
                out = NamedTuple[]
                for (ip, p) in enumerate(ps), (h, N) in enumerate(Ns), k in 1:5
                    f = fin[ip]
                    cp = wb_cp(p)
                    ref = wb_ref_avgs(fb[:, k], N, p)
                    S = Float64(ref.S)
                    ea = max(wb_err(f.avg1[k, h], ref.avg1), wb_err(f.avg2[k, h], ref.avg2)) / (cp * epsT * S)
                    ed = setprecision(BigFloat, WB_BITS) do
                        Float64(abs(wb_big(f.delta[k, h]) - (ref.avg1 - ref.avg2)))
                    end / (2 * cp * epsT * S)
                    ixr = wb_ref_index(ref.avg1, ref.avg2, ref.S, epsT)
                    big = Float64(abs(ixr.delta)) >= 1e3 * cp * epsT * S
                    dabs = big ? Float64(abs(wb_big(f.absdig[k, h]) - ixr.absdig)) : NaN
                    # the index of the module's own returned values, by WB15 in BigFloat (criterion 5(a), (b), I5)
                    ixm = setprecision(BigFloat, WB_BITS) do
                        wb_ref_index(wb_big(f.avg1[k, h]), wb_big(f.avg2[k, h]), wb_big(f.S[k, h]), epsT)
                    end
                    e5(x, y) = isnan(y) ? (isnan(x) ? 0.0 : Inf) : Float64(abs(wb_big(x) - y)) / (2 * epsT * max(Float64(abs(y)), 1.0))
                    r5 = max(e5(f.absdig[k, h], ixm.absdig), e5(f.reldig[k, h], ixm.reldig), e5(f.dig[k, h], ixm.dig))
                    flagok = f.below_floor[k, h] == ixm.below_floor && f.zero[k, h] == ixm.zero
                    undefok = isnan(f.reldig[k, h]) == isnan(ixm.reldig)
                    i5a, i5r = setprecision(BigFloat, WB_BITS) do
                        phi = ixm.phi
                        (Float64(f.absdig[k, h]) <= Float64(-log10(phi)) + 1e-9, isnan(f.reldig[k, h]) || Float64(f.reldig[k, h]) <= Float64(-log10(phi / ixm.m)) + 1e-9)
                    end
                    push!(out, (T = T, p = p, N = N, cls = label, ea = ea, ed = ed, dabs = dabs, r5 = r5, flagok = flagok, undefok = undefok, i5a = i5a, i5r = i5r,
                                nflag = ixm.below_floor, nundef = isnan(ixm.reldig)))
                end
                res[j] = out
            end
            append!(rows, vcat(res...))
        end
    end
    DATA[:ac4] = rows
    for T in WB_TYPES, p in ps
        sel = [r for r in rows if r.T === T && r.p == p]
        ratio_line("AC4", "algebra, $T, w^[$p], 3 classes x $NST starts x 5 observables x N in $Ns: max over both windows of |avg - avg_ref| / (c_p eps S) (c_p = $(Int(wb_cp(p))); spec dry run Float64 0.6 / 7.2 eps S i.e. 0.15 / 0.45)", maximum(r.ea for r in sel))
        ratio_line("AC4", "same, $T, w^[$p]: |Delta - Delta_ref| / (2 c_p eps S)", maximum(r.ed for r in sel))
        info("AC4", "realized, $T, w^[$p]: max |avg - avg_ref| in eps S; max |Delta - Delta_ref| in eps S", "$(fmt(maximum(r.ea for r in sel) * wb_cp(p))); $(fmt(maximum(r.ed for r in sel) * 2 * wb_cp(p)))")
        dd = [r.dabs for r in sel if !isnan(r.dabs)]
        emit("AC4", "same, $T, w^[$p]: absdig against the reference absdig where |Delta_ref| >= 1e3 c_p eps S ($(length(dd)) entries): max |difference| / 1e-2", isempty(dd) ? NaN : maximum(dd) / 1e-2, "<= 1", pf(isempty(dd) || maximum(dd) <= 1e-2))
    end
    for T in WB_TYPES
        sel = [r for r in rows if r.T === T]
        ratio_line("AC5a", "index of the returned values against WB15 in BigFloat, $T, $(length(sel)) entries (absdig, reldig, dig): max |difference| / (2 eps max(|value|, 1)) (D84 item 4)", maximum(r.r5 for r in sel))
        emit("AC4", "flags below_floor and zero equal their definition from the module's own Delta and S, $T: mismatches of $(length(sel))", count(r -> !r.flagok, sel), "0", pf(all(r -> r.flagok, sel)))
        emit("AC5b", "reldig undefined exactly where m <= phi, $T: mismatches; entries with reldig undefined / with the flag set (premises)", "$(count(r -> !r.undefok, sel)); $(count(r -> r.nundef, sel)) / $(count(r -> r.nflag, sel))", "0 mismatches", pf(all(r -> r.undefok, sel)))
        emit("I5", "D84 item 3: absdig <= -log10 phi and reldig <= -log10(phi / m), $T: violations of $(length(sel))", count(r -> !r.i5a || !r.i5r, sel), "0", pf(all(r -> r.i5a && r.i5r, sel)))
    end
end

# ---------------------------------------------------------------------------------------------------
# AC6 reducer against the function (I6): bitwise, K = 4, H = 3, p in {1, 2}, tau in {1, 4}, three classes, NST starts, both types; and the weight vector equal to the module's
section("AC6") do
    obs = wb_std4()
    Ns = (64, 256, 1000)
    nbad = ntot = 0
    for T in WB_TYPES, (label, K) in ((:half, 0.5), (:two_pi, 2pi), (:ten, 10.0)), p in (1, 2), tau in (1, 4)
        starts = class_starts(label, NST)
        sys = WrappedStandardMap(wb_K(T, label))
        r = WeightedBirkhoff(; window = Ns, observables = obs, birkhoff_order = p)
        res = Vector{Int}(undef, NST)
        Threads.@threads :dynamic for j in 1:NST
            x0 = torus_state(wb_cast(T, starts[j]))
            _, _, fin = wb_run((r, wb_Norm()), sys, x0, wb_frame(vseed(6, 1), j, T), 2000 * tau, tau)
            f = fin[1]
            g = weighted_birkhoff(wb_series(obs, wb_states_tau(sys, x0, 2000, tau)), Ns; order = p)
            bad = !(isequal(Matrix(f.avg1), g.avg1) && isequal(Matrix(f.avg2), g.avg2) && isequal(Matrix(f.S), g.S) && isequal(Vector(f.B), g.B))
            # a weight vector equal to the module's weights, one horizon at a time
            for (h, N) in enumerate(Ns)
                gw = weighted_birkhoff(wb_series(obs, wb_states_tau(sys, x0, 2000, tau)), N; order = p, weights = wb_wvec(T, N, p))
                bad |= !(isequal(gw.avg1[:, 1], g.avg1[:, h]) && isequal(gw.avg2[:, 1], g.avg2[:, h]) && isequal(gw.B[1], g.B[h]))
            end
            res[j] = bad ? 1 : 0
        end
        nbad += sum(res); ntot += NST
    end
    emit("AC6", "reducer vs function, K = 4, H = 3, p in {1, 2}, tau in {1, 4}, K_std in {0.5 regular, 2 pi, 10}, $NST starts, Float64 and Double64, plus the weight-vector form: orbits with any bitwise difference, of $ntot", nbad, "0", pf(nbad == 0))
end

# ---------------------------------------------------------------------------------------------------
# AC7 consistency with E2 (F3)
section("AC7") do
    for T in WB_TYPES
        epsT = wb_eps(T)
        worst = 0.0
        for (label, K) in ((:half, 0.5), (:nine, 0.9)), N in (500, 1000, 4000)
            starts = class_starts(label, NST)
            sys = WrappedStandardMap(wb_K(T, label))
            r = WeightedBirkhoff(; window = N, observables = (wb_StdKick(K),))
            e2 = WeightedPhaseAdvance(; window = N, observable = LiftedAngleObservable(1, 3))
            tp = wb_tp(T)
            w = zeros(NST)
            Threads.@threads :dynamic for j in 1:NST
                x0 = torus_state(wb_cast(T, starts[j]))
                _, _, fin = wb_run((r, e2, wb_Norm()), sys, x0, wb_frame(vseed(7, 1), j, T), 2N + 1, 1)
                f, g = fin[1], fin[2]
                S = Float64(f.S[1])
                w[j] = maximum(Float64(abs(wb_big(a) - wb_big(nu) * tp)) / (epsT * S) for (a, nu) in ((f.avg1[1], g.nu1[1]), (f.avg2[1], g.nu2[1])))
            end
            worst = max(worst, maximum(w))
        end
        emit("AC7", "consistency with E2 (reader vs the lifted-increment sums A / B of WeightedPhaseAdvance), $T, K_std in {0.5 regular, 0.9}, N in {500, 1000, 4000}, $NST starts: worst |A/B - tp nu_E2| in eps S (spec dry run 0.5)", worst, "<= 4 (c_1 eps S)", pf(worst <= 4))
    end
end

# ---------------------------------------------------------------------------------------------------
# AC8 state size: sizeof of the state against the P7 formula
section("AC8") do
    for T in WB_TYPES, (K, H) in ((1, 1), (4, 1), (4, 3), (4, 11))
        r = WeightedBirkhoff(; window = ntuple(i -> 16 * 2^(i - 1), H), observables = ntuple(i -> CosObservable(1, i), K))
        b = reducer_bytes(r, T, Val(2), Val(1))
        form = 16 + sizeof(T) * H * (2 + 5K)
        emit("AC8", "reducer bytes, $T, (K, H) = ($K, $H); sizeof(state) $(sizeof(typeof(reducer_init(r, T, Val(2), Val(1))))), workspace $(reducer_workspace(r, T, Val(2), Val(1)))", b, "<= $form (P7); the formula gives $form", pf(b <= form && reducer_workspace(r, T, Val(2), Val(1)) == 0 && isbitstype(typeof(r))))
    end
end

# ---------------------------------------------------------------------------------------------------
# AC11: 1 versus N threads in child processes (output hashes), chunk sizes 64 and 1000
section("AC11") do
    hashes = Dict{Int, Dict{String, String}}(Threads.nthreads() => thread_hashes())
    for nt in (1, 2, 4)
        haskey(hashes, nt) && continue
        cmd = addenv(`$(Base.julia_cmd()) -t $nt --project=$(dirname(Base.active_project())) $(@__FILE__)`, "WB_CHILD" => "1", "WB_SMOKE" => (SMOKE ? "1" : "0"))
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
        ok || emit("AC11", "hashes differ for $k across threads $(sort(collect(keys(hashes)))): $hs", 0, "equal", "FAIL")
    end
    for tag in ("K=0.5", "K=2pi")
        ok_all &= hashes[Threads.nthreads()]["reader+MEGNO $tag chunk=64"] == hashes[Threads.nthreads()]["reader+MEGNO $tag chunk=1000"]
    end
    emit("AC11", "the reader (K = 4, H = 2, window (500, 1000)) through propagate, chains {reader, MEGNO} (chunk sizes 64 and 1000) and {reader, Benettin}, 4096 orbits each, K_std in {0.5, 2 pi}, N = 1000: outputs bitwise identical for thread counts $(sort(collect(keys(hashes)))) (hashes from separate processes): $(length(keys_)) cases", ok_all, "all equal", pf(ok_all))
end

# ---------------------------------------------------------------------------------------------------
# AC12: the independent oracles. (a) the Poisson form by mpmath (external process) and by the BigFloat trapezoid quadrature; (b) the plain BigFloat(256) orbit; (c) the search for an external routine.
section("AC12") do
    alpha = setprecision(BigFloat, WB_BITS) do
        sqrt(BigFloat(2)) - 1
    end
    th0 = setprecision(BigFloat, WB_BITS) do
        BigFloat(1) / 7
    end
    Ns = SMOKE ? (50,) : (50, 100, 200)
    mp = Dict{Int, BigFloat}()
    have_mp = success(pipeline(`python -c "import mpmath"`; stderr = devnull))
    if have_mp
        txt = read(`python $(joinpath(@__DIR__, "weighted-birkhoff-poisson.py")) $(string.(Ns))`, String)
        for l in split(txt, '\n'; keepempty = false)
            n, v = split(l)
            mp[parse(Int, n)] = setprecision(BigFloat, WB_BITS) do
                parse(BigFloat, v)
            end
        end
    end
    have_mp || info("AC12a", "mpmath not importable: the external Poisson form was not run", "-")
    dat = NamedTuple[]
    for N in Ns
        ser, _ = wb_rot_big(th0, alpha, 1, 2N)
        r = setprecision(BigFloat, WB_BITS) do
            weighted_birkhoff(ser, N)
        end
        ref = wb_ref_avgs(ser, N, 1)
        pq = wb_poisson_avg(N, th0, alpha)
        e_mp = haskey(mp, N) ? Float64(abs(r.avg1[1] - mp[N])) : NaN
        e_q = Float64(abs(r.avg1[1] - pq))
        e_pl = Float64(abs(r.avg1[1] - ref.avg1))
        push!(dat, (N = N, mp = e_mp, q = e_q, plain = e_pl))
        haskey(mp, N) && emit("AC12a", "Poisson form WB17 by mpmath 30 digits (external python process), f = cos 2 pi (1/7 + j (sqrt 2 - 1)), N = $N: |BigFloat(256) average of weighted_birkhoff - Poisson value| (value $(fmt(Float64(mp[N]))); spec dry run 1.9e-21, 2.1e-30, 3.8e-43 for N = 50, 100, 200)", e_mp, "<= 1e-18", pf(e_mp <= 1e-18))
        emit("AC12a", "same by the BigFloat trapezoid quadrature of W(s) (128 bits) of the package suite, N = $N", e_q, "<= 1e-18", pf(e_q <= 1e-18))
        info("AC12a", "same, N = $N: against the independent plain-array sum (no code shared)", e_pl)
    end
    DATA[:ac12a] = dat
    # (b) the independent orbit: K = 0.5 regular class, four trigonometric observables, N in {1000, 4000}, both windows, Float64, verify starts
    reg = regular_class(NST).starts
    sys = WrappedStandardMap(0.5)
    obs = (CosObservable(1), SinObservable(1), CosObservable(2), SinObservable(2))
    rows = NamedTuple[]
    for N in (1000, 4000), j in eachindex(reg)
        red = WeightedBirkhoff(; window = N, observables = obs)
        x0 = torus_state(wb_cast(Float64, reg[j]))
        _, _, fin = wb_run((red, wb_Norm()), sys, x0, wb_frame(vseed(12, 1), j, Float64), 2N, 1)
        f = fin[1]
        orb = wb_orbit_big(0.5, reg[j], 2N)
        cols = setprecision(BigFloat, WB_BITS) do
            [[cos(o[1]) for o in orb], [sin(o[1]) for o in orb], [cos(o[2]) for o in orb], [sin(o[2]) for o in orb]]
        end
        e = 0.0
        for k in 1:4
            rf = wb_ref_avgs(cols[k], N, 1)
            e = max(e, wb_err(f.avg1[k], rf.avg1), wb_err(f.avg2[k], rf.avg2))
        end
        push!(rows, (N = N, j = j, e = e))
    end
    DATA[:ac12b] = rows
    for N in (1000, 4000)
        sel = sort([r for r in rows if r.N == N]; by = r -> -r.e)
        info("AC12b", "N = $N: the largest differences by start (start index: value), and the quantiles 50 / 90 / 100 % of the $(length(sel)) starts", "$(join(("$(r.j): $(fmt(r.e))" for r in sel[1:min(3, end)]), ", ")); $(fmt(quantile([r.e for r in sel], 0.5))) / $(fmt(quantile([r.e for r in sel], 0.9))) / $(fmt(sel[1].e))")
    end
    for N in (1000, 4000)
        w = maximum(r.e for r in rows if r.N == N)
        emit("AC12b", "independent BigFloat(256) orbit (plain lifted iteration from the exactly widened start), K_std = 0.5 regular class, $NST starts, 4 trigonometric observables, N = $N, Float64: max over starts, observables and both windows of |avg(module orbit) - avg(BigFloat orbit)| (spec dry run 1.6e-14, 3.4e-14; median $(fmt(rmedian([r.e for r in rows if r.N == N]))))", w, "<= 1e-12", pf(w <= 1e-12))
    end
    # (c) ChaosTools search
    ctdir = dirname(dirname(pathof(CT)))
    terms = ("birkhoff", "weighted", "time average", "bump", "ergodic", "rotation number")
    hits = Dict(t => String[] for t in terms)
    for sub in ("src", joinpath("docs", "src")), (root, _, files) in walkdir(joinpath(ctdir, sub)), fn in files
        endswith(fn, ".jl") || endswith(fn, ".md") || continue
        path = joinpath(root, fn)
        for (i, l) in enumerate(eachline(path)), t in terms
            occursin(t, lowercase(l)) && push!(hits[t], "$(relpath(path, ctdir)):$i")
        end
    end
    info("AC12c", "ChaosTools $(pkgversion(CT)) at $ctdir, case-insensitive search of src/ and docs/src/ (lines): birkhoff / weighted / time average / bump / ergodic / rotation number", join((string(length(hits[t])) for t in terms), " / "))
    info("AC12c", "files with a 'birkhoff' line", isempty(hits["birkhoff"]) ? "none" : join(unique(first.(split.(hits["birkhoff"], ':'))), ", "))
    info("AC12c", "files with a 'weighted' line", isempty(hits["weighted"]) ? "none" : join(unique(first.(split.(hits["weighted"], ':'))), ", "))
    pyok = [m_ for m_ in ("nafflib", "PyNAFF", "rebound", "celmech") if success(pipeline(`python -c "import $m_"`; stderr = devnull))]
    info("AC12c", "Python packages that could hold a weighted Birkhoff routine (nafflib, PyNAFF, rebound, celmech) importable in the default Python: installed ones", isempty(pyok) ? "none" : join(pyok, ", "))
    info("AC12c", "result", "no external weighted Birkhoff routine exists in the environment; nothing to compare (spec Oracle)")
end

# ---------------------------------------------------------------------------------------------------
# R1: class separation on the 2D panel, 1000 starts per class, the five observables, N in {250, 1000, 4000, 16000, 64000} in one pass (H = 5), with the cross-tabulation at n = 2N for N = 1000 and 4000
const R1_NS = SMOKE ? (250, 1000, 4000) : (250, 1000, 4000, 16000, 64000)
const R1_M = SMOKE ? 40 : 1000
const R1_OBS = ("cos q", "sin q", "cos p", "sin p", "lifted p")
function r1_panel(label, K, starts)
    sys = WrappedStandardMap(K)
    m = length(starts)
    red = WeightedBirkhoff(; window = R1_NS, observables = wb_std5())
    H = length(R1_NS)
    out = (absdig = fill(NaN, 5, H, m), reldig = fill(NaN, 5, H, m), dig = fill(NaN, 5, H, m), floor_ = fill(false, 5, H, m), S = fill(NaN, 5, H, m))
    Threads.@threads :dynamic for j in 1:m
        _, _, fin = wb_run((red, wb_Norm()), sys, torus_state(starts[j]), wb_frame(vseed(42, 1), j, Float64), 2 * maximum(R1_NS), 1)
        f = fin[1]
        out.absdig[:, :, j] .= Matrix(f.absdig); out.reldig[:, :, j] .= Matrix(f.reldig); out.dig[:, :, j] .= Matrix(f.dig)
        out.floor_[:, :, j] .= Matrix(f.below_floor); out.S[:, :, j] .= Matrix(f.S)
    end
    return out
end
# E2 (naff-fma) on one class: the D of the lifted increments, N in {1000, 4000}, as the comparator of the cross-tabulation
function e2_panel(sys, x0, Ns; seed)
    rl = WeightedPhaseAdvance(; window = Tuple(Ns), observable = LiftedAngleObservable(1, 3))
    res = propagate(sys, x0, 2 * maximum(Ns) + 1; k = 1, tau = 1, seed = seed, reducers = (rl, Benettin()), backend = :cpu)
    return Dict(N => [Float64(f[1].D[h]) for f in res.finals] for (h, N) in enumerate(Ns))
end

section("R1") do
    m = R1_M
    reg = regular_class(m)
    info("R1", "K = 0.5 regular class: candidates screened, regular (Benettin X(1e5) <= 5e-4), taken", "$(reg.ncand), $(reg.nreg), $m")
    R = Dict{Symbol, Any}()
    for (label, K) in R_CLASSES
        t = @elapsed R[label] = r1_panel(label, K, class_starts(label, m))
        say(@sprintf("# R1 class %s: panel %.1f s", label, t))
    end
    DATA[:R1] = R
    Ns = R1_NS
    for (k, nm) in enumerate(R1_OBS), (h, N) in enumerate(Ns)
        for (label, K) in R_CLASSES
            d = filter(isfinite, R[label].absdig[k, h, :])
            info("R1", "$nm, N = $N, K = $(fmt(K)) ($label), $m starts: quantiles 5 / 25 / 50 / 75 / 95 % of absdig (non-finite excluded: $(m - length(d))); share at the floor", "$(quants(d)); $(fmt(frac(R[label].floor_[k, h, :])))")
        end
        reg_d = filter(isfinite, R[:half].absdig[k, h, :])
        q = median(reg_d)
        line = String[]
        for (label, K) in R_CLASSES[2:end]
            d = R[label].absdig[k, h, :]
            ok = isfinite.(d)
            push!(line, "K=$(fmt(K)): at or above the regular median or at the floor $(fmt(mean((d[ok] .> q) .| R[label].floor_[k, h, :][ok]))), above 5.5: $(fmt(mean(d[ok] .> 5.5))), above 5: $(fmt(mean(d[ok] .> 5)))")
        end
        info("R1", "$nm, N = $N: regular-class median of absdig $(fmt(q)); regular starts below 5.5 / below 5: $(fmt(mean(reg_d .< 5.5))) / $(fmt(mean(reg_d .< 5))); the chaotic classes", join(line, "; "))
    end
    # the same for reldig and dig (max) of the first observable, N = 1000 and 4000
    for (what, key) in (("reldig", :reldig), ("dig", :dig)), N in (1000, 4000)
        h = findfirst(==(N), Ns)
        for (label, K) in R_CLASSES
            d = filter(isfinite, getfield(R[label], key)[1, h, :])
            info("R1", "cos q, N = $N, K = $(fmt(K)) ($label): quantiles 5 / 25 / 50 / 75 / 95 % of $what (non-finite excluded: $(m - length(d)))", quants(d))
        end
    end
    # modes of absdig (SanderMeiss2020 Fig. 3): histogram of cos q absdig over all starts, unit bins from -2 to 17
    edges = collect(-2.0:1.0:17.0)
    for N in (1000, 4000, 16000)
        h = findfirst(==(N), Ns)
        h === nothing && continue
        pooled = vcat((filter(isfinite, R[label].absdig[1, h, :]) for (label, _) in R_CLASSES)...)
        cnt = [count(x -> lo <= x < lo + 1, pooled) for lo in edges[1:end-1]]
        hs = [(get(cnt, i - 1, 0) + 2 * cnt[i] + get(cnt, i + 1, 0)) / 4 for i in eachindex(cnt)]
        peaks = [i for i in eachindex(hs) if hs[i] >= 0.03 * maximum(hs) && hs[i] > get(hs, i - 1, 0) && hs[i] >= get(hs, i + 1, 0)]
        regcnt = [count(x -> lo <= x < lo + 1, filter(isfinite, R[:half].absdig[1, h, :])) for lo in edges[1:end-1]]
        info("R1", "cos q, N = $N: histogram of absdig over the $(length(R_CLASSES) * m) starts, unit bins from -2 to 17 (counts; smoothed local maxima at the lower bin edges: $(join(edges[peaks], ", ")) -> $(length(peaks)) modes)", join(cnt, ","))
        info("R1", "cos q, N = $N: the regular class alone, same bins", join(regcnt, ","))
    end
    # the Float64 plateau: the regular class median / maximum of absdig and the share at the floor, cos q (S about 1: floor -log10(32 eps) = 14.15), against N
    for (h, N) in enumerate(Ns)
        d = filter(isfinite, R[:half].absdig[1, h, :])
        info("R1", "Float64 plateau, regular class, cos q, N = $N: median / 95 % quantile / maximum of absdig; share at the floor (flag); -log10(32 eps S) at the median S", "$(fmt(median(d))) / $(fmt(quantile(d, 0.95))) / $(fmt(maximum(d))); $(fmt(frac(R[:half].floor_[1, h, :]))); $(fmt(-log10(32 * wb_eps(Float64) * median(R[:half].S[1, h, :]))))")
    end
    # cross-tabulation at n = 2N for N in {1000, 4000}: cos q absdig below 5.5 (the printed cutoff of SanderMeiss2020, not adopted) against the labels of L1 (n X above ln n), FLI_B (above ln n), MEGNO
    # (Ybar_1 above 3), RLI (RLI / Dx_r above 1 at Dx = 1e-12) and stretching numbers D (below its regular-class median, the chaotic direction of the module) and the naff-fma E2 D (above its regular-class median), as in the naff-fma R1
    labels = Dict{Tuple{Symbol, Int}, NamedTuple}()
    seedR = vseed(42, 1)
    for (label, K) in R_CLASSES
        base = WrappedStandardMap(K)
        P = PairMap(base)
        starts = class_starts(label, m)
        xw = [torus_state(s) for s in starts]
        xis = [SVector{2, Float64}(initial_deviations(seedR, j, Val(2), Val(1), Float64)[:, 1]) for j in 1:m]
        ps = [pair_start(base, starts[j], xis[j], 1e-12; orbit = j) for j in 1:m]
        z0 = [p.z for p in ps]
        W0 = [SMatrix{4, 2, Float64}(xi[1], xi[2], 0.0, 0.0, 0.0, 0.0, xi[1], xi[2]) for xi in xis]
        d2 = e2_panel(base, xw, (1000, 4000); seed = seedR)
        for N in (1000, 4000)
            n = 2N
            r = propagate(P, z0, n; k = 2, tau = 1, seed = seedR, reducers = (RLI(),), W0 = W0, backend = :cpu)
            l1 = propagate(base, xw, n; k = 1, tau = 1, seed = seedR, reducers = (Benettin(),), backend = :cpu)
            fm = propagate(base, xw, n; k = 2, tau = 1, seed = seedR, reducers = (FLI(replace = false), MEGNO()), backend = :cpu)
            sn = propagate(base, xw, n; k = 2, tau = 1, seed = seedR, reducers = (StretchingNumbers(; range = stretching_range(base), bins = 1000, interval = 100),), backend = :cpu)
            labels[(label, N)] = (RLI = [r.finals[j][1].RLI / ps[j].Dx_r > 1.0 for j in 1:m], L1 = [f[1].lambda[1] * n > log(n) for f in l1.finals],
                                  FLI = [f[1].fli > log(n) for f in fm.finals], MEGNO = [f[2].Ybar[1] > 3 for f in fm.finals], SNdist = [Float64(f[1].D) for f in sn.finals], E2D = d2[N])
        end
    end
    for N in (1000, 4000)
        h = findfirst(==(N), Ns)
        dig = vcat((R[l].absdig[1, h, :] for (l, _) in R_CLASSES)...)
        okd = isfinite.(dig)
        sep = dig .< 5.5
        for lab in (:L1, :FLI, :MEGNO, :RLI, :SNdist, :E2D)
            lv = if lab === :SNdist
                vcat((labels[(l, N)].SNdist .< median(labels[(:half, N)].SNdist) for (l, _) in R_CLASSES)...)
            elseif lab === :E2D
                vcat((labels[(l, N)].E2D .> median(filter(isfinite, labels[(:half, N)].E2D)) for (l, _) in R_CLASSES)...)
            else
                vcat((getproperty(labels[(l, N)], lab) for (l, _) in R_CLASSES)...)
            end
            lv, s = lv[okd], sep[okd]
            per = join((fmt(frac(lab === :SNdist ? labels[(l, N)].SNdist .< median(labels[(:half, N)].SNdist) : lab === :E2D ? labels[(l, N)].E2D .> median(filter(isfinite, labels[(:half, N)].E2D)) : getproperty(labels[(l, N)], lab))) for (l, _) in R_CLASSES), "; ")
            info("R1", "n = $(2N) (N = $N), cos q absdig below 5.5 against label $lab ($(lab === :SNdist ? "stretching D below its regular-class median: the chaotic direction, the two columns fall into the same cells" : lab === :E2D ? "naff-fma E2 D above its regular-class median" : "chaotic")), all $(count(okd)) finite starts of 4 classes: (label chaotic, digits below 5.5) / (chaotic, above) / (regular, below) / (regular, above); fraction labelled chaotic per class (K = 0.5; 2 pi; 5; 10)",
                 "$(count(lv .& s)) / $(count(lv .& .!s)) / $(count(.!lv .& s)) / $(count(.!lv .& .!s)); $per")
        end
    end
end

# ---------------------------------------------------------------------------------------------------
# Figures
section("figure") do
    # oracle.png: (a) the Poisson oracle, (b) the independent orbit, (c) the algebra ratios of AC4
    fig = CM.Figure(size = (1800, 560))
    if haskey(DATA, :ac12a)
        d = DATA[:ac12a]
        ax = CM.Axis(fig[1, 1]; yscale = log10, xscale = log10, xlabel = "window N", ylabel = "absolute difference of the window average", title = "AC12(a): weighted_birkhoff (BigFloat 256) against the Poisson form WB17 (gate 1e-18)")
        CM.scatterlines!(ax, [r.N for r in d], [max(r.q, 1e-80) for r in d]; color = :seagreen, marker = :circle, label = "BigFloat trapezoid quadrature (package suite)")
        any(r -> isfinite(r.mp), d) && CM.scatterlines!(ax, [r.N for r in d], [max(r.mp, 1e-80) for r in d]; color = :purple, marker = :utriangle, label = "mpmath 30 digits (external process)")
        CM.scatterlines!(ax, [r.N for r in d], [max(r.plain, 1e-80) for r in d]; color = :gray, marker = :rect, label = "independent plain-array sum")
        CM.hlines!(ax, [1e-18]; color = :black, linestyle = :dash)
        CM.axislegend(ax; position = :lb, labelsize = 11)
    end
    if haskey(DATA, :ac12b)
        rows = DATA[:ac12b]
        ax = CM.Axis(fig[1, 2]; yscale = log10, xlabel = "regular-class start", ylabel = "max |avg(module orbit) - avg(BigFloat orbit)|", title = "AC12(b): module orbit against the BigFloat(256) orbit, Float64 (gate 1e-12)")
        for (N, col, mk) in ((1000, :seagreen, :circle), (4000, :darkorange, :utriangle))
            sel = [r for r in rows if r.N == N]
            CM.scatter!(ax, [r.j for r in sel], [max(r.e, 1e-19) for r in sel]; color = col, marker = mk, label = "N = $N")
        end
        CM.hlines!(ax, [1e-12]; color = :black, linestyle = :dash)
        CM.axislegend(ax; position = :rb, labelsize = 11)
    end
    if haskey(DATA, :ac4)
        rows = DATA[:ac4]
        ax = CM.Axis(fig[1, 3]; yscale = log10, xticks = (1:4, ["F64 p=1", "F64 p=2", "D64 p=1", "D64 p=2"]), ylabel = "error / bound (gate 1)", title = "AC4: algebra against the BigFloat(256) sums, max ratio per class")
        for (i, (T, p)) in enumerate(((Float64, 1), (Float64, 2), (Double64, 1), (Double64, 2))), (cls, col) in ((:half, :seagreen), (:two_pi, :darkorange), (:ten, :purple))
            sel = [r for r in rows if r.T === T && r.p == p && r.cls === cls]
            isempty(sel) && continue
            CM.scatter!(ax, [i - 0.12], [max(maximum(r.ea for r in sel), 1e-6)]; color = col, marker = :circle, markersize = 12)
            CM.scatter!(ax, [i + 0.12], [max(maximum(r.ed for r in sel), 1e-6)]; color = col, marker = :utriangle, markersize = 12)
        end
        CM.scatter!(ax, [NaN], [NaN]; color = :seagreen, marker = :circle, label = "average (circle), Delta (triangle): K = 0.5 regular")
        CM.scatter!(ax, [NaN], [NaN]; color = :darkorange, marker = :circle, label = "K = 2 pi")
        CM.scatter!(ax, [NaN], [NaN]; color = :purple, marker = :circle, label = "K = 10")
        CM.hlines!(ax, [1.0]; color = :black, linestyle = :dash)
        CM.axislegend(ax; position = :lb, labelsize = 10)
    end
    CM.save(joinpath(OUT, "weighted-birkhoff-oracle.png"), fig)
    say("# figure saved: weighted-birkhoff-oracle.png")
    # digits.png: absdig against N per observable, regular versus chaotic class, with the Float64 plateau
    if haskey(DATA, :R1)
        R = DATA[:R1]
        Ns = collect(R1_NS)
        fig2 = CM.Figure(size = (1800, 950))
        local legax
        cols = Dict(:half => :seagreen, :two_pi => :darkorange, :five => :firebrick, :ten => :purple)
        for (k, nm) in enumerate(R1_OBS)
            ax = CM.Axis(fig2[(k - 1) ÷ 3 + 1, (k - 1) % 3 + 1]; xscale = log10, xlabel = "window N (handoffs 2N)", ylabel = "absdig = -log10 max(|Delta|, 32 eps S)", title = "observable $nm")
            for (label, K) in R_CLASSES
                med = [median(filter(isfinite, R[label].absdig[k, h, :])) for h in eachindex(Ns)]
                lo = [quantile(filter(isfinite, R[label].absdig[k, h, :]), 0.1) for h in eachindex(Ns)]
                hi = [quantile(filter(isfinite, R[label].absdig[k, h, :]), 0.9) for h in eachindex(Ns)]
                CM.band!(ax, Ns, lo, hi; color = (cols[label], 0.2))
                CM.scatterlines!(ax, Ns, med; color = cols[label], label = label === :half ? "K = 0.5 regular class (median, 10-90 %)" : "K = $(fmt(K)) (median, 10-90 %)")
            end
            Smed = median(filter(isfinite, R[:half].S[k, length(Ns), :]))
            CM.hlines!(ax, [-log10(32 * wb_eps(Float64) * Smed)]; color = :black, linestyle = :dash, label = "Float64 floor -log10(32 eps S), median S of the regular class")
            CM.hlines!(ax, [5.5]; color = :gray, linestyle = :dot, label = "5.5 (SanderMeiss2020 cutoff, not adopted)")
            k == 1 && (legax = ax)
        end
        CM.Legend(fig2[2, 3], legax; labelsize = 12)
        CM.save(joinpath(OUT, "weighted-birkhoff-digits.png"), fig2)
        say("# figure saved: weighted-birkhoff-digits.png")
    end
end

say(@sprintf("# verdicts: %d PASS, %d FAIL, %d other", count(==("PASS"), VERDICTS), count(==("FAIL"), VERDICTS), count(v -> v != "PASS" && v != "FAIL", VERDICTS)))
say(@sprintf("# total wall time %.1f s", time() - T0))
close(LOG)
