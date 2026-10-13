# Stage 04 verification of the precision spec (workspaces/module-build/stages/01-spec/output/precision-spec.md).
# Run from the repo root: julia --project=verify -t 4 verify/precision.jl
# Reuses the independent reference helpers of the package suite (BigFloat orbit, M_n, Gamma_n, n_H).
# Line format: AC<n> | <quantity> | <value> | <bound> | PASS/FAIL/INFO/STAGE04-GPU/DROPPED
using Chaoslab, StaticArrays, DoubleFloats, StableRNGs, Statistics, Printf
using CairoMakie

include(joinpath(@__DIR__, "..", "test", "precision", "helpers.jl"))   # pr_* helpers, PR_SEED = 20261012
const OUT = joinpath(@__DIR__, "..", "workspaces", "module-build", "stages", "04-verify", "output")
mkpath(OUT)

row(n, q, v, b, verdict) = println("AC$n | $q | $v | $b | $verdict")
fmt(x) = @sprintf("%.4g", x)
pf(c) = c ? "PASS" : "FAIL"
println("# precision verify; threads = $(Threads.nthreads()), julia $(VERSION), seed base $PR_SEED")

# ---- AC1 parameter conversion (Float32 dropped, D24) --------------------------------------------------
let
    refK(n) = n == "pi/2" ? BigFloat(pi) / 2 : n == "2pi" ? BigFloat(pi) * 2 : n == "0.1" ? parse(BigFloat, "0.1") : parse(BigFloat, "0.001")
    build(T, n) = n == "pi/2" ? exact(T, pi) / 2 : n == "2pi" ? exact(T, pi) * 2 : n == "0.1" ? exact(T, "0.1") : exact(T, "0.001")
    worst = Dict{Any,Float64}()
    for n in ("pi/2", "2pi", "0.1", "1e-3"), T in (Float64, Double64, BigFloat)
        Kref = setprecision(() -> refK(n), BigFloat, 512)
        K = setprecision(() -> build(T, n), BigFloat, 256)
        epsT = setprecision(() -> Float64(eps(T)), BigFloat, 256)
        setprecision(BigFloat, 512) do
            err = abs(pr_big(K) - Kref)
            bound = T === Float64 ? BigFloat(eps(K)) / 2 : 2 * BigFloat(epsT) * Kref
            worst[T] = max(get(worst, T, 0.0), Float64(err / bound))
        end
    end
    for T in (Float64, Double64, BigFloat)
        row(1, "max |K_T-K_ref|/bound, T=$T", fmt(worst[T]), "1", pf(worst[T] <= 1))
    end
end

# ---- AC2 exact widening --------------------------------------------------------------------------------
let
    rng = StableRNG(PR_SEED + 2); mant = 0x000fffffffffffff; vals = Float64[]
    for i in 1:10^4
        r = i % 4
        v = r == 0 ? (rand(rng) * 2 - 1) * 1e3 :
            r == 1 ? reinterpret(Float64, (rand(rng, UInt64) & mant) | 0x1) * (rand(rng, Bool) ? 1 : -1) :
            r == 2 ? 2pi + (rand(rng) - 0.5) * 1e-3 : reinterpret(Float64, rand(rng, UInt64))
        push!(vals, isfinite(v) ? v : 1.0)
    end
    bd = 0; bb = 0
    setprecision(BigFloat, 256) do
        for v in vals
            d = exact(Double64, v); b = exact(BigFloat, v)
            bd += !(Float64(d) === v && pr_big(d) == BigFloat(v)); bb += !(Float64(b) === v && b == BigFloat(v))
        end
    end
    row(2, "round-trip failures of 10^4 values, Double64 / BigFloat(256)", "$bd / $bb", "0 / 0", pf(bd == 0 && bb == 0))
end

# ---- AC3 transcendentals in Double64, effective epsilon --------------------------------------------------
let
    rng = StableRNG(PR_SEED + 3); epsD = Float64(eps(Double64)); N = 10^4
    sys = StandardMap{Double64}(one(Double64))
    worst = Dict("sin" => 0.0, "cos" => 0.0, "log" => 0.0)       # max err / max(1,|f|), absolute units
    setprecision(BigFloat, 256) do
        for (R, cnt) in ((2pi, N ÷ 2), (1e3, N ÷ 2)), _ in 1:cnt
            x = pr_rand_d64(rng, -R, R); z = SVector{2,Double64}(x, zero(Double64)); xb = pr_big(x)
            for (nm, got, want) in (("sin", forward(sys, z)[2], sin(xb)), ("cos", tangent(sys, z)[2, 1], cos(xb)))
                worst[nm] = max(worst[nm], Float64(abs(pr_big(got) - want)) / max(1.0, Float64(abs(want))))
            end
        end
        for _ in 1:N
            x = pr_d64_from(rng, exp10(-3 + 6 * rand(rng))); want = log(pr_big(x))
            worst["log"] = max(worst["log"], Float64(abs(pr_big(log(x)) - want)) / max(1.0, Float64(abs(want))))
        end
    end
    for nm in ("sin", "cos", "log")
        w = worst[nm]
        row(3, "$nm max err/max(1,|f|) in units of eps_D=2^-104 (n=$N points)", fmt(w / epsD), "16", pf(w / epsD <= 16))
        println("AC3 | $nm effective eps = max err/max(1,|f|) | 2^$(@sprintf("%.2f", log2(w))) = $(fmt(w)) | 2^-100 (amendment) | ", pf(w <= 2.0^-100))
    end
end

# ---- AC4 pathwise orbit agreement (+ data for the figure) ---------------------------------------------------
fig = Dict{Symbol,Any}()
let
    ics = pr_ics(PR_SEED + 4, 200)
    refs = setprecision(BigFloat, 256) do
        [pr_reference(BigFloat(5), x0, p0, PR_NMAX) for (x0, p0) in ics]
    end
    MG = [pr_bounds(r) for r in refs]
    fig[:ics] = ics; fig[:MG] = MG
    for (T, fl) in ((Float64, 15), (Double64, 40))
        epsT = Float64(eps(T)); nHs = Int[]; worst = 0.0; low = 0
        curves = Vector{Vector{Float64}}(); errs = Vector{Vector{Float64}}()
        setprecision(BigFloat, 256) do
            for (i, (x0, p0)) in enumerate(ics)
                ref = refs[i]; M, G = MG[i]
                nH = pr_nH(M, G, epsT); push!(nHs, nH)
                # full-length run for the figure (to n = 100), ratios only counted up to nH
                zs = pr_run(T, exact(T, 5), x0, p0, PR_NMAX)
                e = [Float64(pr_err(zs[n+1], ref.xs[n+1], ref.ps[n+1])) for n in 1:PR_NMAX]
                push!(errs, e); push!(curves, [e[n] / (epsT * M[n+1] * G[n+1]) for n in 1:PR_NMAX])
                nH < fl && (low += 1; continue)
                worst = max(worst, maximum(curves[end][1:nH]))
            end
        end
        fig[Symbol(T, :_errs)] = errs; fig[Symbol(T, :_ratios)] = curves; fig[Symbol(T, :_nH)] = nHs
        capped = count(==(PR_NMAX), nHs)
        row(4, "$T realised n_H min/median/max (capped at $PR_NMAX: $capped orbits)", "$(minimum(nHs))/$(median(nHs))/$(maximum(nHs))", ">= $fl", minimum(nHs) >= fl ? "PASS" : "see exclusions")
        row(4, "$T orbits with n_H < $fl (excluded)", low, "<= 10 (5 percent)", pf(low <= 10))
        row(4, "$T max |e_n|/(eps M_n Gamma_n), n <= n_H", fmt(worst), "8", pf(worst <= 8))
    end
end

# ---- AC5 BigFloat 256 vs 512 ---------------------------------------------------------------------------------
let
    ics = pr_ics(PR_SEED + 5, 200); eps256 = Float64(setprecision(() -> eps(BigFloat), BigFloat, 256))
    refs = setprecision(BigFloat, 512) do
        [pr_reference(BigFloat(5), x0, p0, PR_NMAX) for (x0, p0) in ics]
    end
    runs = setprecision(BigFloat, 256) do
        [pr_run(BigFloat, exact(BigFloat, 5), x0, p0, PR_NMAX) for (x0, p0) in ics]
    end
    worst = 0.0; e10 = 0.0
    setprecision(BigFloat, 512) do
        for (ref, zs) in zip(refs, runs)
            M, G = pr_bounds(ref); ow = 0.0
            for n in 1:PR_NMAX
                e = Float64(pr_err(zs[n+1], ref.xs[n+1], ref.ps[n+1]))
                ow = max(ow, e / (eps256 * M[n+1] * G[n+1])); n == 10 && (e10 = max(e10, e))
            end
            worst = max(worst, ow)
        end
    end
    row(5, "max |e_n|/(eps_256 M_n Gamma_n), n <= 100", fmt(worst), "8", pf(worst <= 8))
    row(5, "max |e_10|", fmt(e10), "1e-60", pf(e10 <= 1e-60))
end

# ---- AC6 tangent growth (amendment: at most 10 of 200 orbits may exceed) -------------------------------------
let
    ics = pr_ics(PR_SEED + 6, 200); rng = StableRNG(PR_SEED + 60); epsD = Float64(eps(Double64))
    ratios = Float64[]; nmaxs = Int[]
    setprecision(BigFloat, 256) do
        for (x0, p0) in ics
            a, b = rand(rng) - 0.5, rand(rng) - 0.5; u0 = (a, b) ./ hypot(a, b)
            ref = pr_reference(BigFloat(5), x0, p0, PR_NMAX); M, G = pr_bounds(ref)
            nmax = min(pr_nH(M, G, epsD), 40); push!(nmaxs, nmax)
            Sd = pr_growth_lib(Double64, exact(Double64, 5), x0, p0, u0, nmax); Sr = pr_growth_ref(BigFloat(5), x0, p0, u0, nmax)
            push!(ratios, maximum(Float64(abs(pr_big(Sd[n+1]) - Sr[n+1])) / (64 * 6 * epsD * M[n+1] * G[n+1]) for n in 1:nmax))
        end
    end
    s = sort(ratios)
    row(6, "realised horizon min(n_H(D64),40): min/median/max", "$(minimum(nmaxs))/$(median(nmaxs))/$(maximum(nmaxs))", "-", "INFO")
    row(6, "max |S_n(D64)-S_n(ref)|/(64(1+K) eps M_n Gamma_n)", fmt(maximum(ratios)), "1 (per orbit)", maximum(ratios) <= 1 ? "PASS" : "see count")
    row(6, "99th percentile / median ratio", "$(fmt(s[198])) / $(fmt(median(ratios)))", "-", "INFO")
    row(6, "orbits with ratio > 1 (of 200)", count(>(1), ratios), "<= 10 (amendment)", pf(count(>(1), ratios) <= 10))
end

# ---- AC7 fallback matrix, AC8 recording ----------------------------------------------------------------------
let
    cr = "kernel does not compile for the element type"; bad = 0; cases = 0
    for T in (Float64, Double64, BigFloat), req in (:auto, :gpu, :cpu)
        probe = (T === BigFloat && req !== :cpu) ? pr_probe_boom : pr_probe_noload
        r = resolve_backend(req, T; probe = probe)
        want = req === :cpu ? "" : T === BigFloat ? cr : "CUDA.jl not loadable"
        cases += 1; bad += !(r.backend === :cpu && r.reason == want && parse_precision(precision_label(T))[1] === T)
    end
    for T in (Float64, Double64), req in (:auto, :gpu)
        cases += 3
        r = resolve_backend(req, T; probe = pr_probe_nofunc);    bad += !((r.backend, r.reason) == (:cpu, "no functional device"))
        r = resolve_backend(req, T; probe = pr_probe_nocompile); bad += !((r.backend, r.reason) == (:cpu, cr))
        r = resolve_backend(req, T; probe = pr_probe_good);      bad += !((r.backend, r.reason) == (:gpu, ""))
    end
    row(7, "fallback-matrix mismatches over $cases cases (9 type x request cases + probe cases; Float32 rows dropped D24)", bad, "0", pf(bad == 0))
    bad8 = 0
    bad8 += !(precision_label(Float64) == "Float64" && precision_label(Double64) == "Double64")
    for b in (64, 128, 256, 512)
        l = setprecision(() -> precision_label(BigFloat), BigFloat, b)
        bad8 += !(l == "BigFloat($b)" && occursin(r"^(Float64|Double64|BigFloat\(\d+\))$", l) && parse_precision(l) == (BigFloat, b))
    end
    row(8, "precision-string failures (Float64, Double64, BigFloat(64/128/256/512))", bad8, "0", pf(bad8 == 0))
end

# ---- AC9 / AC10 GPU, AC11 CPU part ---------------------------------------------------------------------------
row(9, "GPU vs CPU orbit agreement (Float64, Double64)", "see verify/gpu.jl", "2 c eps M_n Gamma_n", "STAGE04-GPU")
row(10, "Float32 GPU gate (>= 4x Float64 CPU)", "Float32 dropped from M2", "-", "DROPPED D24")

function ens!(out, ::Type{T}, K, x0, p0, nsteps, threaded) where {T}
    sys = StandardMap{T}(K)
    body(i) = (z = SVector{2,T}(T(x0[i]), T(p0[i])); for _ in 1:nsteps; z = forward(sys, z); end; out[i] = z[1])
    if threaded
        Threads.@threads :static for i in eachindex(x0); body(i); end
    else
        for i in eachindex(x0); body(i); end
    end
    out
end
let
    NORB, NST = 4000, 500; rng = StableRNG(PR_SEED + 11)
    x0 = 2pi .* rand(rng, NORB); p0 = 2pi .* rand(rng, NORB)
    times = Dict{Tuple{Any,Bool},Float64}()
    for T in (Float64, Double64, BigFloat), thr in (false, true)
        out = Vector{T}(undef, NORB); K = exact(T, 5)
        ens!(out, T, K, x0, p0, 5, thr)                                   # warm-up
        times[(T, thr)] = minimum(@elapsed(ens!(out, T, K, x0, p0, NST, thr)) for _ in 1:3)
    end
    for thr in (false, true)
        nt = thr ? Threads.nthreads() : 1
        tf = times[(Float64, thr)]
        row(11, "time ratio Double64/Float64, $nt thread(s), $(NORB)x$(NST) steps (Float64 $(fmt(tf*1e3)) ms)", fmt(times[(Double64, thr)] / tf), "reported, no threshold", "PASS")
        row(11, "time ratio BigFloat(256)/Float64, $nt thread(s)", fmt(times[(BigFloat, thr)] / tf), "reported, no threshold", "PASS")
    end
    row(11, "thread scaling 1 -> $(Threads.nthreads()) threads, speedup Float64 / Double64 / BigFloat", join((fmt(times[(T, false)] / times[(T, true)]) for T in (Float64, Double64, BigFloat)), " / "), "-", "INFO")
    row(11, "GPU Double64/Float64 ratio", "see verify/gpu.jl", "-", "STAGE04-GPU")
end

# ---- AC12 closure and allocation -----------------------------------------------------------------------------
let
    bad = 0; allocmax = 0
    for T in (Float64, Double64, BigFloat)
        setprecision(BigFloat, 256) do
            sys = StandardMap{T}(exact(T, 5)); z = SVector{2,T}(exact(T, 1.25), exact(T, -0.75))
            zf = forward(sys, z); M = tangent(sys, z); acc = add(CompensatedSum{T}(), log(pr_norm(zf)))
            bad += !(zf isa SVector{2,T} && M isa SMatrix{2,2,T} && acc isa CompensatedSum{T} && total(acc) isa T)
        end
    end
    for T in (Float64, Double64)
        sys = StandardMap{T}(exact(T, 5)); z = SVector{2,T}(exact(T, 1.25), exact(T, -0.75)); acc = CompensatedSum{T}(); x = exact(T, 0.3)
        pr_loop_forward(sys, z, 2); pr_loop_tangent(sys, z, 2); pr_loop_add(acc, x, 2)
        allocmax = max(allocmax, pr_allocs(pr_loop_forward, sys, z, 100), pr_allocs(pr_loop_tangent, sys, z, 100), pr_allocs(pr_loop_add, acc, x, 100))
    end
    row(12, "closure violations (Float64, Double64, BigFloat)", bad, "0", pf(bad == 0))
    row(12, "max allocated bytes in 100-step hot loops (Float64, Double64)", allocmax, "0", pf(allocmax == 0))
end

# ---- AC13 conversion consequence -----------------------------------------------------------------------------
let
    ics = pr_ics(PR_SEED + 13, 100); s1 = StandardMap{Double64}(exact(Double64, pi) * 2); s2 = StandardMap{Double64}(Double64(2pi))
    ns = Int[]
    for (x0, p0) in ics
        z1 = SVector{2,Double64}(Double64(x0), Double64(p0)); z2 = z1; n = 0
        while n < 200
            z1 = forward(s1, z1); z2 = forward(s2, z2); n += 1
            maximum(abs.(z1 - z2)) > 1e-8 && break
        end
        push!(ns, n)
    end
    row(13, "median n* (separation > 1e-8), min/max", "$(median(ns)) ($(minimum(ns))/$(maximum(ns)))", "reported", "INFO")
    row(13, "fraction of ICs with n* <= 60", fmt(count(<=(60), ns) / length(ns)), ">= 0.8", pf(count(<=(60), ns) >= 80))
end

# ---- Figure --------------------------------------------------------------------------------------------------
let
    MG = fig[:MG]; ns = 1:PR_NMAX
    nHf = fig[:Float64_nH]; i0 = argmin(abs.(nHf .- median(nHf)))    # orbit with median Float64 horizon
    M, G = MG[i0]
    f = Figure(size = (1100, 480))
    ax = Axis(f[1, 1]; yscale = log10, xlabel = "step n", ylabel = "max-norm orbit error vs BigFloat(256)",
              title = "Orbit $i0 (K = 5): error and bound c eps M_n Gamma_n, c = 8")
    for (T, col) in ((Float64, :royalblue), (Double64, :crimson))
        e = max.(fig[Symbol(T, :_errs)][i0], 1e-300); epsT = Float64(eps(T))
        lines!(ax, ns, e; color = col, label = "$T error")
        lines!(ax, ns, [8 * epsT * M[n+1] * G[n+1] for n in ns]; color = col, linestyle = :dash, label = "$T bound")
        vlines!(ax, [fig[Symbol(T, :_nH)][i0]]; color = col, linestyle = :dot)
    end
    hlines!(ax, [1e-4]; color = :gray, linestyle = :dot, label = "1e-4 horizon level")
    ylims!(ax, 1e-34, 1e3); axislegend(ax; position = :rb)
    ax2 = Axis(f[1, 2]; yscale = log10, xlabel = "step n", ylabel = "|e_n| / (eps M_n Gamma_n)",
               title = "Ratio to the bound, 200 orbits, n <= n_H (bound = 8)")
    for (T, col) in ((Float64, (:royalblue, 0.25)), (Double64, (:crimson, 0.25)))
        rs = fig[Symbol(T, :_ratios)]; nH = fig[Symbol(T, :_nH)]
        for i in eachindex(rs)
            nH[i] >= 1 && lines!(ax2, 1:nH[i], max.(rs[i][1:nH[i]], 1e-6); color = col)
        end
        lines!(ax2, [NaN], [NaN]; color = col[1], label = string(T))
    end
    hlines!(ax2, [8.0]; color = :black, linestyle = :dash, label = "c = 8")
    ylims!(ax2, 1e-4, 1e2); axislegend(ax2; position = :rb)
    save(joinpath(OUT, "precision-orbit-error.png"), f)
    println("# figure: ", normpath(joinpath(OUT, "precision-orbit-error.png")))
end
