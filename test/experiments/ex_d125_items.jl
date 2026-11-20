# Amendment D125, stage 02: D125c (systems), D125j (items, hashes, frames, plan), D125k (full spectrum), D125l (diagnostics), D125m (indicators).
# Seeds: d125_seed('c' | 'j' | 'k' | 'l' | 'm', sub). Fixtures are tiny (orbits of the hand panel, horizons 100 and 160). RED reasons: missing keys of the new
# kinds (the reader refuses them), missing EX functions (lazy proxy), missing system kinds in plan items.

@testset "D125c systems: exact parameters, decoupled limits, wrapped lift, flow loop, wrapped flow refused (E44, I41)" begin
    rng = StableRNG(d125_seed('c'))
    # (1) every parameter equals the exact value rounded once, 0 ulp against the 1024-bit evaluation, for the three precisions
    for label in ("Float64", "Double64", "BigFloat(256)"), kind in ("froeschle", "chain")
        ok = ex_with_prec(label) do T
            ps = kind == "froeschle" ? ["0.5", "0.7", "0.3"] : ["3", "1.1", "0.2", "open"]
            s = d125_system(kind, ps, T)
            vals = kind == "froeschle" ? [s.nu, s.kappa, s.mu] : [s.K, s.eps]
            texts = kind == "froeschle" ? ps : ps[2:3]
            all(zip(vals, texts)) do (v, t)
                ref = setprecision(() -> d125_big(t), BigFloat, 1024)
                setprecision(() -> T(ref), BigFloat, label == "BigFloat(256)" ? 256 : 53) == v
            end
        end
        @test ok
    end
    # (2) Y3 at mu = 0 and the chain at eps = 0 equal the standard map, bitwise per coordinate pair, 1000 steps (conjugation C of the coupled-maps spec for Y3)
    for t in 1:3
        z = [2pi * rand(rng) for _ in 1:4]
        F = Chaoslab.FroeschleMap(0.5, 0.7, 0.0); S1 = Chaoslab.StandardMap(-0.5); S2 = Chaoslab.StandardMap(-0.7)
        x = Chaoslab.SVector{4, Float64}(z...); a = Chaoslab.SVector(z[1] + z[3], z[3]); b = Chaoslab.SVector(z[2] + z[4], z[4])
        same = true
        for _ in 1:1000
            x = Chaoslab.forward(F, x); a = Chaoslab.forward(S1, a); b = Chaoslab.forward(S2, b)
            same &= (x[1] + x[3], x[3]) == (a[1], a[2]) && (x[2] + x[4], x[4]) == (b[1], b[2])
        end
        @test same
        C = Chaoslab.StandardMapChain(3, 1.1, 0.0; boundary = :open); Ss = Chaoslab.StandardMap(1.1)
        w = Chaoslab.SVector{6, Float64}((2pi .* rand(rng, 6))...); pairs = [Chaoslab.SVector(w[j], w[3 + j]) for j in 1:3]
        same = true
        for _ in 1:1000
            w = Chaoslab.forward(C, w); pairs = [Chaoslab.forward(Ss, p) for p in pairs]
            same &= all(j -> (w[j], w[3 + j]) == (pairs[j][1], pairs[j][2]), 1:3)
        end
        @test same
    end
    # (3) the wrapped orbit lifted equals the lifted orbit within the wrapped bound of the coupled-maps spec (10 steps, 16 eps per step and wrap count, reading U8)
    for kind in ("froeschle", "chain")
        ps = kind == "froeschle" ? ["0.5", "0.7", "0.3"] : ["2", "1.1", "0.2", "open"]
        sl = d125_system(kind, ps, Float64); sw = d125_system(kind, ps, Float64; wrapped = true); d = kind == "froeschle" ? 4 : 4
        z = Chaoslab.SVector{d, Float64}((2pi .* rand(rng, d))...); u = Chaoslab.torus_state(sw, z)
        worst = 0.0
        for n in 1:10
            z = Chaoslab.forward(sl, z); u = Chaoslab.forward(sw, u)
            worst = max(worst, maximum(abs.(Chaoslab.lift(sw, u) .- z)))
        end
        @test worst <= 16 * EX_EPS * 10 * 8 * 2pi * 10
    end
    # (4) the one-step map of the flow equals an independent loop of the integrator step, bitwise, 1000 steps, each integrator of the fixtures
    for row in d125_integrators()
        x0 = first(d125_flow_starts(EX_SEED, "hh1", "1/12", ("-0.3", "0.4", "-0.25", "0.25"), 1)[1])
        m = d125_flow_map(row, Float64)
        got = d125_states(m, x0, [1000], Float64)[1]
        meth = d125_method(row["method"], row["scheme"], Float64); h = d125_f64(row["step"]); z = Chaoslab.SVector{4, Float64}(x0...)
        for _ in 1:1000
            z = Chaoslab.Integrators.advance(meth, Chaoslab.HenonHeilesFlow(), z, h)
        end
        @test got == Float64.(collect(z))
    end
    # (5) a wrapped request on a flow is refused through the module (the new reader: coordinate_mode wrapped on a cell of kind henon_heiles)
    run = d125_block_spectra(); run["arm"][2]["coordinate_mode"] = ["wrapped"]
    @test d125_block_refused(run, ["wrapped", "henon_heiles", "hh1"]; tag = "wrapped-flow")
end

@testset "D125j items, hashes, frames and plan on the fixtures (R35, E55, E56, I45)" begin
    fx = d125_fix("std"); plan = d125_plan(fx)
    # hand count: maps arm 6 orbits x 5 methods, flow arm 3 orbits x 3 methods x 2 integrators, controls arm 2 x 3 methods
    @test length(plan.items) == 6 * 5 + 3 * 3 * 2 + 2 * 3 && length(plan.chunks) == 4 * 5 + 2 * 3 + 3
    # the text of every item equals the hand-written text byte for byte (reading U5 on the line order); two fixtures: the standard one and B
    for it in plan.items
        it.system in ("froeschle", "chain", "henon_heiles") || continue
        @test EX.item_text(it) == d125_item_text(it; integ = it.system == "henon_heiles" ? d125_integ_text(d125_row(it)) : nothing)
        @test EX.item_hash(it) == ex_sha256(d125_item_text(it; integ = it.system == "henon_heiles" ? d125_integ_text(d125_row(it)) : nothing))
    end
    # a changed frame label, integrator step or cell gives a different hash of the same position
    base = Set(EX.item_hash(it) for it in plan.items)
    p2 = d125_variant("j-frame", r -> r["arm"][1]["frames"] = ["b"])
    p3 = d125_variant("j-step", r -> r["integrator"][1]["step"] = "0.03125")
    p4 = d125_variant("j-cell", r -> r["arm"][1]["cells"] = ["y3a"])
    h2 = [EX.item_hash(it) for it in p2.items if it.system in ("froeschle", "chain")]
    @test !any(in(base), h2) && length(unique(h2)) == length(h2)
    @test count(in(base), [EX.item_hash(it) for it in p3.items if it.system == "henon_heiles"]) == 9               # only the leapfrog items change (9 of 18)
    @test length(p4.items) < length(plan.items) && issubset(Set(EX.item_hash(it) for it in p4.items), base)
    # five frames of five arms give five distinct items; the same label in two arms merges
    five = ["a", "b", "c", "d", "e"]
    arms = Any[d125_arm("f$l", "P", ["L5"], [100]; cells = ["y3a"], vectors = "d", frames = [l]) for l in five]
    r5 = d125_run(arms); p5 = d125_plan(d125_fix("j-five", r5))
    @test length(p5.items) == 5 * 3 && length(unique(EX.item_hash(it) for it in p5.items)) == 15
    push!(arms, d125_arm("dup", "P", ["L5"], [100]; cells = ["y3a"], vectors = "d", frames = ["a"]))
    @test length(d125_plan(d125_fix("j-dup", d125_run(arms))).items) == 15
    # the frame of label "a" is the project's seeded rule, bitwise, k = 2, 4, d; the first k columns of the frame of d; two labels two frames
    for (d, ks) in ((4, (2, 4)), (6, (2, 4, 6)), (10, (2, 4, 10)))
        run_seed = EX_RUN_SEED
        full = EX.seeded_frame(Float64, run_seed, 3, "a", Val(d), Val(d))                                       # RED: missing function
        for k in ks
            @test EX.seeded_frame(Float64, run_seed, 3, "a", Val(d), Val(k)) == Chaoslab.initial_deviations(run_seed, 3, Val(d), Val(k), Float64)
            @test EX.seeded_frame(Float64, run_seed, 3, "a", Val(d), Val(k)) == full[:, 1:k]
        end
        id = parse(UInt64, ex_sha256("chaoslab-frame-1|3|b")[1:16]; base = 16)                                   # first 8 bytes, big-endian
        @test EX.seeded_frame(Float64, run_seed, 3, "b", Val(d), Val(2)) == Chaoslab.initial_deviations(run_seed, id, Val(d), Val(2), Float64)
        @test EX.seeded_frame(Float64, run_seed, 3, "b", Val(d), Val(d)) != full
    end
    # the plan is the same for 1, 4 and all threads and for blocks that list keys in another order
    ref = EX.plan_text(plan)
    @test all(n -> EX.plan_text(d125_plan(fx; ntasks = n)) == ref, (1, 4, Threads.nthreads()))
    shuf = d125_block_spectra(); shuf["arm"][1]["cells"] = reverse(shuf["arm"][1]["cells"]); shuf["arm"][2]["integrators"] = reverse(shuf["arm"][2]["integrators"])
    ps = d125_plan(d125_fix("j-shuf", shuf))
    @test sort([EX.item_hash(it) for it in ps.items]) == sort([EX.item_hash(it) for it in plan.items])
    # old items unchanged when a new arm is added: ordinals, keys, hashes of the chunks of the existing items
    more = d125_block_spectra(); push!(more["arm"], d125_arm("extra", "P", ["L5"], [100]; cells = ["y3a"], vectors = "d", frames = ["c"]))
    pm = d125_plan(d125_fix("j-more", more))
    @test [EX.chunk_text(c) for c in pm.chunks[1:length(plan.chunks)]] == [EX.chunk_text(c) for c in plan.chunks]
    @test [EX.item_hash(it) for it in pm.items[1:length(plan.items)]] == [EX.item_hash(it) for it in plan.items]
    # a plan holding zero items and a spec or estimator used by no arm are errors naming them (R35 item 4)
    @test d125_block_refused(d125_run(Any[d125_arm("none", "P", ["L5"], [100]; cells = ["y3a"], vectors = "d", frames = ["a"], integrators = ["lf"])];
                                      estimators = Any[d125_cqr(id = "unused")]), ["unused"]; tag = "unused-estimator")
end

@testset "D125k full spectrum: prefix, spread, Benettin gap, identical orbit, setup check, hot loop (R36, I46, I47)" begin
    run = d125_campaign("std")
    @test run.res !== nothing
    byorb = Dict{Tuple{Int, String}, Dict{String, Any}}()
    rng = StableRNG(d125_seed('k'))
    nchecked = 0
    d125_each(run) do it, raw, c, file
        it.system in ("froeschle", "chain", "henon_heiles") || return
        o = it.orbit; pre = "orbits/$o/"
        lam = raw[pre * "lyap/" * it.method]; t = raw[pre * "t"]; st = raw[pre * "state"]
        d = size(st, 2); k = it.method == "L1" ? 1 : d
        @test size(lam) == (length(t), k) && size(st) == (length(t), d)                                           # n x d state, n x k exponents
        h = d125_item_h(it)
        # (a) a stored value equals the direct run to that step, 0 difference; the first, the last and three seeded rows
        for r in unique([1, length(t), (1 + floor(Int, rand(rng) * length(t)) for _ in 1:3)...])
            n = Int(t[r])
            direct = d125_lambda(d125_item_system(it), d125_start(it), n, it.method, k, parse(Int, it.run_seed), o; tau = it.tau)
            @test lam[r, :] == direct ./ (it.tau * h)
            nchecked += 1
        end
        get!(byorb, (o, it.integrator), Dict{String, Any}())[it.method] = (lam = lam, state = st, h = h, N = Int(t[end]), tau = it.tau)
    end
    @test nchecked >= 50
    for ((o, integ), ms) in byorb
        # (d) the state datasets of all methods of one orbit are bitwise equal
        @test length(unique(reinterpret(UInt64, vec(m.state)) |> hash for m in values(ms))) == 1
        d = size(first(values(ms)).state, 2)
        d > 0 || continue
        haskey(ms, "L5") || continue
        h = ms["L5"].h; N = ms["L5"].N; tau = ms["L5"].tau
        # (b) spread of lambda_i over L3, L4, L4b, L5 at most B_QR(d, N) = 16 d^2 eps / (tau h) (calibrated, plan table K24)
        qr = [ms[m].lam[end, :] for m in ("L3", "L4", "L4b", "L5") if haskey(ms, m) && size(ms[m].lam, 2) == d]
        length(qr) >= 2 && @test maximum(maximum(v[i] for v in qr) - minimum(v[i] for v in qr) for i in 1:d) <= 16 * d^2 * EX_EPS / (tau * h)
        # (c) |lambda_1^L1 - lambda_1^L5| <= 16 eps (N + G_1) / (N h) with G_1 = N h lambda_1 (reading U9: B_S of the benettin spec with D35)
        if haskey(ms, "L1")
            G1 = abs(N * h * ms["L5"].lam[end, 1])
            @test abs(ms["L1"].lam[end, 1] - ms["L5"].lam[end, 1]) <= 16 * EX_EPS * (N + G1) / (N * h)
        end
    end
    # (e) the extended precision items agree with the Float64 item over the first setup_check_steps steps within setup_check_tol; a small run on one cell
    arms = Any[d125_arm("prec", "P", ["L5"], [100]; cells = ["y3a"], vectors = "d", frames = ["a"], precisions = ["Float64", "Double64", "BigFloat(256)"])]
    rp = d125_campaign("k-prec", d125_run(arms))
    pt = Dict{String, Matrix{Float64}}()
    d125_each(rp) do it, raw, c, file
        it.orbit == 1 && (pt[it.precision] = raw["orbits/1/state"])
    end
    @test length(pt) == 3 && all(p -> maximum(abs.(pt[p][1:3, :] .- pt["Float64"][1:3, :])) <= 1e-10, ("Double64", "BigFloat(256)"))
    # (f) the hot loop of one orbit allocates 0 bytes after warm-up, for each kind at d = 4 and d = 10 (assumed function EX.hot_loop_bytes, plan)
    for (kind, d) in (("froeschle", 4), ("chain", 4), ("chain", 10), ("henon_heiles", 4))
        @test EX.hot_loop_bytes(kind, d) == 0
    end
end

@testset "D125l diagnostics for d: eps_pair, eps_sum, thresholds, flags, drift (E57, E58)" begin
    rng = StableRNG(d125_seed('l'))
    for d in (4, 6, 10), t in 1:333
        lam = sort(randn(rng, d); rev = true) .* 10.0 .^ rand(rng, -3:2)
        pair = maximum(abs(lam[i] + lam[d + 1 - i]) for i in 1:(d ÷ 2)); sm = abs(sum(lam))
        @test EX.eps_pair(lam) == pair && EX.eps_sum(lam) == sm                                                   # 1000 spectra in all, 0 difference
    end
    # thresholds: the diagnostics spec values for d with the factor d for sum (E58, reading U10: G_max = 10, tau = 1, M_max = d in the fixture); pair for d > 2 is the
    # finite-time form of D126, pair_const / (n tau h) + the rounding term of sum (reading U10, D126; n = the first recorded grid step, here 200; h = 1 or 1/20)
    T = Float64
    for d in (4, 6, 10), h in (1.0, 0.05)
        th = EX.diagnostic_thresholds(d, T; G_max = 10.0, tau = 1, M_max = d, n = 200, h = h, pair_const = D126_PAIR_CONST)
        @test th.det == Chaoslab.det_threshold(4.0, T, d) && th.symp == Chaoslab.symp_threshold(4.0, T, d)
        @test th.sum == d * Chaoslab.sum_threshold(16.0, T, 10.0, 1, d)
        @test th.pair == D126_PAIR_CONST / (200 * 1 * h) + d * Chaoslab.sum_threshold(16.0, T, 10.0, 1, d)
    end
    # synthetic series above and below each threshold flag as expected (EX.diagnostics_flag(name, values, threshold) -> first exceedance step or 0)
    for d in (4, 10), name in (:det, :symp, :orth, :sum, :pair)
        th = getproperty(EX.diagnostic_thresholds(d, T; G_max = 10.0, tau = 1, M_max = d, n = 200, h = 1.0, pair_const = D126_PAIR_CONST), name)
        below = fill(0.5 * th, 20); above = vcat(fill(0.5 * th, 7), 2.0 * th, fill(0.5 * th, 12))
        @test EX.diagnostics_flag(name, below, th) == 0 && EX.diagnostics_flag(name, above, th) == 8
    end
    # drift reads the momentum components: Y3 (p1, p2), chain (p_i), none for flows; a position beyond the bound does not raise it
    bound = 1e3
    for (kind, d, mom) in (("froeschle", 4, 3:4), ("chain", 6, 4:6), ("henon_heiles", 4, 1:0))
        st = zeros(5, d); st[:, 1] .= 10 * bound                                                                  # a position far beyond the bound
        @test !EX.drift_flag(kind, st, bound)
        isempty(mom) && continue
        st2 = zeros(5, d); st2[4, mom[end]] = 2 * bound
        @test EX.drift_flag(kind, st2, bound)
    end
    # calibration table (K24): the realized ratios are in the plan section; here the bounds are the pinned values, not changed after the data
end

@testset "D125m indicators on the new systems: columns equal the library finals, chain rules, refusals (R37)" begin
    spec = d86_spec("tangent"; id = "tg", members = ["FLI", "MEGNO", "SALI"])
    gali = d86_spec("tangent"; id = "ga", members = ["GALI"], gali_k = 4)
    ofl = d86_spec("tangent"; id = "of", members = ["FLI", "OFLI", "SALI"])
    arms = Any[d125_arm("ind-maps", "P", nothing, [100]; inds = ["tg", "ga"], cells = ["y3a", "ch2"], frames = ["a"]),
               d125_arm("ind-flow", "P", nothing, [160]; inds = ["tg", "ga", "of"], cells = ["hh1"], frames = ["a"], integrators = ["lf"])]
    run = d125_campaign("m-ind", d125_run(arms; specs = Any[spec, gali, ofl]))
    n = 0
    d125_each(run) do it, raw, c, file
        it.system in ("froeschle", "chain", "henon_heiles") || return
        o = it.orbit; d = it.system == "chain" ? 4 : 4; sys = d125_item_system(it); x = Chaoslab.SVector{d, Float64}(d125_start(it)...)
        N = it.horizon; col(name) = raw["orbits/$o/indicators/$(it.method)/$name"]
        if it.method == "ga"
            W = Chaoslab.initial_deviations(parse(Int, it.run_seed), o, Val(d), Val(4), Float64)
            fin = Chaoslab.Tangent.reducer_final(Chaoslab.GALI(), Chaoslab.Tangent.propagate_orbit(sys, x, W, N, it.tau, (Chaoslab.GALI(),))[3][1])
            for j in 2:4
                @test col("GALI.gali$j")[end] == Float64(fin.gali[j - 1])
            end
            @test col("GALI.L")[end] == Float64(fin.L)                                                             # K = d, the log form is recorded
        elseif it.method == "of"
            W = Chaoslab.initial_deviations(parse(Int, it.run_seed), o, Val(d), Val(2), Float64)
            reds = (Chaoslab.FLI(; replace = false), Chaoslab.OFLI(Chaoslab.HenonHeilesFlow(); replace = false), Chaoslab.SALI())
            fins = Chaoslab.Tangent.propagate_orbit(sys, x, W, N, it.tau, reds)[3]
            fo = Chaoslab.Tangent.reducer_final(reds[2], fins[2])
            @test col("OFLI.ofli")[end] == Float64(fo.ofli) && col("OFLI.ofli1")[end] == Float64(fo.ofli1)
        else
            W = Chaoslab.initial_deviations(parse(Int, it.run_seed), o, Val(d), Val(2), Float64)
            reds = (Chaoslab.FLI(; replace = false), Chaoslab.MEGNO(; replace = false), Chaoslab.SALI())
            fins = Chaoslab.Tangent.propagate_orbit(sys, x, W, N, it.tau, reds)[3]
            @test col("SALI.sali")[end] == Float64(Chaoslab.Tangent.reducer_final(reds[3], fins[3]).sali)
            @test col("FLI.fli")[end] == Float64(Chaoslab.Tangent.reducer_final(reds[1], fins[1]).fli)
        end
        n += 1
    end
    @test n >= 12
    # chain rules of E23 extended and the refusals
    bad = Dict(
        "sali-not-last" => d86_spec("tangent"; id = "x1", members = ["SALI", "FLI"]),
        "ofli-on-map" => d86_spec("tangent"; id = "x2", members = ["OFLI"]),
        "gali-with-member" => d86_spec("tangent"; id = "x3", members = ["GALI", "SALI"], gali_k = 3),
        "gali-k-over-d" => d86_spec("tangent"; id = "x4", members = ["GALI"], gali_k = 5),
        "gali-k-missing" => d86_spec("tangent"; id = "x5", members = ["GALI"]),
        "gali-k-no-gali" => d86_spec("tangent"; id = "x6", members = ["SALI"], gali_k = 3))
    for (nm, sp) in bad
        cells = nm == "ofli-on-map" ? ["y3a"] : nm == "gali-k-over-d" ? ["y3a"] : ["hh1"]
        r = d125_run(Any[d125_arm("b", "P", nothing, [100]; inds = [sp["id"]], cells = cells, frames = ["a"], integrators = ["lf"])]; specs = Any[sp])
        @test d125_block_refused(r, [sp["id"], "members", "gali_k", "OFLI", "GALI"]; tag = "ind-" * nm)
    end
    # a frame with a column parallel to the field is refused by OFLI before the first chunk (frame check, D67): an explicit hand frame is not in the block,
    # so this refusal is tested at the check function EX.check_ofli_frame(flow, start, frame)
    z = first(d125_flow_starts(EX_SEED, "hh1", "1/12", ("-0.3", "0.4", "-0.25", "0.25"), 1)[1])
    f = Chaoslab.vector_field(Chaoslab.HenonHeilesFlow(), Chaoslab.SVector{4, Float64}(z...)); fr = hcat(collect(f) ./ sqrt(sum(abs2, f)), [0.0, 1, 0, 0])
    @test ex_rejects(() -> EX.check_ofli_frame(Chaoslab.HenonHeilesFlow(), z, fr), ["OFLI", "parallel"])
end
