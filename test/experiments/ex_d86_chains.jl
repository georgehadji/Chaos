# Amendment D86, criteria D86e (chains), D86f (checks before the stream), D86g (tangent family), D86h (pair).
# Spec: experiments-spec.md, "Amendment D86": E23, E24, E30, R13, R14, I19 to I22, O16, O18. API assumed: experiments-test-plan.md, "API assumed by the D86 tests".

const D86_TG_CHAINS = [["SALI"], ["FLI", "SALI"], ["MEGNO", "SALI"], ["FLI", "MEGNO", "SALI"], ["FLI"], ["MEGNO"], ["FLI", "MEGNO"]]
const D86_PR_CHAINS = [["RLI"], ["RLI", "MEGNO"], ["MEGNO", "RLI"]]

# Every chain of E23 as its own spec on the 9 orbits of S_D (Float64, lifted), one chain in wrapped mode, Double64 on the 2 extra starts,
# the pairs in both modes and with both frame sources, L1 twins of every orbit setting (I21).
function d86_chain_block()
    specs = Dict{String, Any}[]
    for (i, m) in enumerate(D86_TG_CHAINS)
        push!(specs, d86_spec("tangent"; id = "t$i", members = m))
    end
    for (i, m) in enumerate(D86_PR_CHAINS)
        push!(specs, d86_spec("pair"; id = "p$i", members = m))
    end
    push!(specs, d86_spec("pair"; id = "p4", members = ["RLI"], frame_source = "direction"))
    push!(specs, d86_spec("pair"; id = "p5", members = ["RLI", "MEGNO"], frame_source = "direction"))
    push!(specs, d86_spec("stretching"; id = "s1"))
    tg = ["t$i" for i in 1:7]
    arms = [d86_arm("tg", "S_D", [tg; "s1"], ["Float64"], [100]),
            d86_arm("tgw", "S_D", ["t4"], ["Float64"], [100]; coordinate_mode = ["wrapped"]),
            d86_arm("tg2", "S_D", ["t4"], ["Float64"], [80]; tau = [2, 4]),
            d86_arm("pr", "S_D", ["p1", "p2", "p3", "p4", "p5"], ["Float64"], [100]; coordinate_mode = ["lifted", "wrapped"]),
            d86_arm("l1", "S_D", String[], ["Float64"], [100]; methods = ["L1"], coordinate_mode = ["lifted", "wrapped"]),
            d86_arm("l1t", "S_D", String[], ["Float64"], [80]; methods = ["L1"], tau = [2, 4]),
            d86_arm("d64", "h4", tg, ["Double64"], [100]; coordinate_mode = ["wrapped"]),
            d86_arm("d64l1", "h4", String[], ["Double64"], [100]; methods = ["L1"], coordinate_mode = ["wrapped"])]
    for a in arms
        isempty(a["indicators"]) && delete!(a, "indicators")
    end
    return d86_block(arms, specs)
end
d86_chain_run() = d86_run("chains", d86_chain_block())

# provenance of a member as stored: the HDF5 group `provenance/<member>`, one attribute per published key (R15; correction 4 of the D86 stage 02 checkpoint)
d86_prov(raw, m) = Dict(String(k[(length("provenance/$m@") + 1):end]) => v for (k, v) in raw if startswith(k, "provenance/$m@"))

@testset "D86e chains: forms, record, provenance, mutation" begin
    run = d86_chain_run()
    @test run.res.status == "complete"
    plan = run.plan
    chains = Dict("t$i" => m for (i, m) in enumerate(D86_TG_CHAINS))
    merge!(chains, Dict("p$i" => m for (i, m) in enumerate(D86_PR_CHAINS)), Dict("p4" => ["RLI"], "p5" => ["RLI", "MEGNO"], "s1" => ["stretching"]))
    # the derived forms of E23 for every chain of the block, from the independent rule (last member replaces)
    forms = Dict("t1" => [true], "t2" => [false, true], "t3" => [false, true], "t4" => [false, false, true], "t5" => [true], "t6" => [true], "t7" => [false, true],
                 "p1" => [true], "p2" => [false, true], "p3" => [false, true], "p4" => [true], "p5" => [false, true], "s1" => [true])
    @test all(id -> d86_forms(chains[id]) == forms[id], keys(chains))
    seen = Set{String}()
    d86_each(run) do it, orbit, raw, c, file
        it.method in keys(chains) || return
        push!(seen, it.method)
        m = chains[it.method]
        spec = d86_spec_of(plan, it)
        @test raw["@chain"] == [string(x, "(replace=", f, ")") for (x, f) in zip(m, forms[it.method])]
        @test raw["@members"] == m && raw["@kind"] == it.kind && raw["@method"] == it.method
        @test raw["@spec"] == d86_spec_text(spec) && raw["@spec_sha256"] == d86_spec_hash(spec) == it.spec
        orbit == first(plan.items[c.items]).orbit || return              # the provenance is a root attribute: check it once per file
        for mm in m
            mm == "stretching" && continue
            p = d86_prov(raw, mm)
            want = ex_Tg.reducer_provenance(d86_red(mm; replace = false))
            @test !isempty(p) && all(k -> string(p[String(k)]) == string(getfield(want, k)), keys(want))
        end
    end
    @test seen == Set(keys(chains))
    # the key `indicator` is not merged (correction 4, E23, D86e): each member's value is read from its own group, so the three members of t4 show three values
    f4 = first(f for f in ex_chunk_files(run.store) if plan.items[plan.chunks[ex_ordinal(f)].items[1]].method == "t4")
    raw4 = d86_raw(f4)
    vals = [get(d86_prov(raw4, mm), "indicator", nothing) for mm in ("FLI", "MEGNO", "SALI")]
    @test all(!isnothing, vals) && length(unique(vals)) == 3
    # read back equal: the package reader returns the same attributes and datasets as the raw read
    ch = EX.read_chunk(f4)
    @test ch.attrs["chain"] == raw4["@chain"] && ch.attrs["members"] == raw4["@members"] && ch.attrs["spec"] == raw4["@spec"]
    # mutation check: a build that gives SALI the replacing form while first, or two replacing members, changes the later members; D86g would catch it
    it = first(i for i in plan.items if i.method == "t4" && i.coordinate_mode == "lifted" && i.tau == 1)
    n = 100
    good = d86_stand_alone(it, "FLI", n)
    run_chain(red) = ex_with_prec(it.precision) do T
        sys, x, W = d86_orbit_args(it, T)
        ex_Tg.propagate_orbit(sys, x, W, n, it.tau, red)[3]
    end
    right = run_chain((Chaoslab.FLI(replace = false), Chaoslab.MEGNO(replace = false), Chaoslab.SALI()))
    @test ex_Tg.reducer_final(Chaoslab.FLI(), right[1]).fli == good.fli                              # the correct chain: bitwise the stand-alone member
    sab1 = run_chain((Chaoslab.SALI(), Chaoslab.FLI(), Chaoslab.MEGNO()))                           # SALI replacing and first
    @test ex_Tg.reducer_final(Chaoslab.FLI(), sab1[2]).fli != good.fli
    sab2 = run_chain((Chaoslab.FLI(), Chaoslab.MEGNO(replace = false), Chaoslab.SALI()))           # two replacing members
    @test ex_Tg.reducer_final(Chaoslab.MEGNO(), sab2[2]).Y != d86_stand_alone(it, "MEGNO", n).Y
end

@testset "D86f checks before the stream" begin
    # (a) a dependent frame supplied through the test seam is refused before the first chunk, naming kind initial frame, with 0 files
    fx = d86_fx("f-a", d86_one(d86_spec("tangent", members = ["FLI", "MEGNO", "SALI"])))
    dep(info, W) = (W2 = Matrix(W); W2[:, 2] .= W2[:, 1]; typeof(W)(W2))
    store = ex_store("d86-f-a")
    r = try
        ex_campaign(fx; store = store, hooks = EX.RunHooks(clock = ex_clock(info -> 0.0), frame = (kind, info, W) -> kind === :frame ? dep(info, W) : W))
        nothing
    catch e
        e
    end
    @test r !== nothing && occursin("initial frame", sprint(showerror, r))
    @test d86_nofiles(store)
    # the frames the module builds pass the frame checks of every member that has one, for 10 000 seeded directions
    rng = StableRNG(d86_seed('f', 1))
    bad = 0
    for _ in 1:10_000
        txt = string(2pi * rand(rng))
        W = SMatrix{2, 2, Float64}(EX.deviation_frame(Float64, txt, 2))
        try
            Chaoslab.check_frame(Chaoslab.SALI(), W; orbit = 0); Chaoslab.check_frame(Chaoslab.FLI(), W; orbit = 0)
            Chaoslab.check_frame(Chaoslab.MEGNO(), W; orbit = 0)
            Chaoslab.check_frame(Chaoslab.StretchingNumbers(range = (-1.0, 1.0), bins = 10, interval = 1, min_per_cell = 0), W; orbit = 0)
        catch
            bad += 1
        end
    end
    @test bad == 0
    # (b) a pair block frame with unequal vectors in its blocks is refused for the pair by the RLI frame check: the chunk fails and names the pair id
    fb = d86_fx("f-b", d86_one(d86_spec("pair"; members = ["RLI"]); orbits = "S_L"))
    unequal(kind, info, W) = kind === :pair_frame ? (W2 = Matrix(W); W2[3:4, 2] .*= 0.5; typeof(W)(W2)) : W
    storeb = ex_store("d86-f-b")
    resb = ex_campaign(fb; store = storeb, hooks = EX.RunHooks(clock = ex_clock(info -> 0.0), frame = unequal))
    chb = EX.read_chunk(first(ex_chunk_files(storeb)))
    @test chb.attrs["status"] == "failed" && any(r -> occursin("initial frame", r) && occursin("orbit 1", r), chb.attrs["failure_reasons"])
    # (c) at plan time: named by spec id and orbit id, 0 files
    orbs = ex_panel_orbits(:std)
    only(o for o in orbs if o["id"] == 5)["system"] = "henon"            # a system without a derived range
    pp = ex_write_panel(ex_path("d86-f-c-panel.toml"); orbits = orbs)
    pr = ex_write_prereg_run(ex_path("d86-f-c-prereg.md"), d86_one(d86_spec("stretching")))
    @test d86_rej(() -> EX.plan_run(pr, pp; code_version = EX_CODE), "st", "5", "range")
    orbs2 = ex_panel_orbits(:std)
    only(o for o in orbs2 if o["id"] == 6)["x0_dec"] = "1e999"          # a start whose observables are not finite
    pp2 = ex_write_panel(ex_path("d86-f-c2-panel.toml"); orbits = orbs2)
    pr2 = ex_write_prereg_run(ex_path("d86-f-c2-prereg.md"), d86_one(d86_spec("birkhoff"); mode = "wrapped"))
    @test d86_rej(() -> EX.plan_run(pr2, pp2; code_version = EX_CODE), "bk", "6")
    # a run shorter than 2 W + 1 handoffs for phase_advance (the horizon 64 of windows 16 and 32) is refused before anything runs
    short = d86_one(d86_spec("phase_advance"))
    short["arm"][1]["horizons"] = [64]
    @test d86_rej(() -> EX.plan_run(d86_fx("f-c3", short).prereg, d86_fx("f-c3", short).panel; code_version = EX_CODE), "horizons")
    # an observable index beyond the dimension: the grammar allows coordinates 1 and 2 only, so the run check of the reducer is the gate (Indicators level). The limit is 2 d
    # (the wrapped state of a map of dimension d has 2 d entries, naff_reducer.jl: "observable indices must lie in 1:2d"), so the index 5 is beyond for d = 2 and 3 is not
    @test_throws ArgumentError Chaoslab.check_run(Chaoslab.WeightedPhaseAdvance(window = 16, observable = Chaoslab.PlanarObservable(1, 5)), 1, 33; d = 2)
    @test Chaoslab.check_run(Chaoslab.WeightedPhaseAdvance(window = 16, observable = Chaoslab.PlanarObservable(1, 3)), 1, 33; d = 2) === nothing
    # (d) a pair with Dx_r outside [Dx / 2, 2 Dx] fails its chunk with kind separation; the rerun repeats it. Independent of the direction u (A8): both lifted coordinates of
    # the start are 1e9, where one ulp is 1.2e-7 > 2 Dx = 2e-8, so B = A + Dx u rounds back to A in both coordinates whatever u is (Dx_r = 0)
    orbs3 = ex_panel_orbits(:std)
    o1 = only(o for o in orbs3 if o["id"] == 1)
    o1["x0_dec"] = "1.0e9"; o1["p0_dec"] = "1.0e9"; o1["x0"] = 1.0e9; o1["p0"] = 1.0e9
    pp3 = ex_write_panel(ex_path("d86-f-d-panel.toml"); orbits = orbs3)
    pr3 = ex_write_prereg_run(ex_path("d86-f-d-prereg.md"), d86_one(d86_spec("pair"; members = ["RLI"]); orbits = "S_L"))
    fxd = (; run = nothing, panel = pp3, prereg = pr3, tag = "f-d", orbits = orbs3)
    store = ex_store("d86-f-d")
    resd = ex_campaign(fxd; store = store)
    chd = EX.read_chunk(first(ex_chunk_files(store)))
    @test chd.attrs["status"] == "failed" && chd.attrs["attempt"] == 2
    @test any(r -> occursin("separation", r) && occursin("orbit 1", r), chd.attrs["failure_reasons"]) && any(r -> occursin("separation", r), chd.attrs["first_attempt_reasons"])
    @test length(resd.failed) == 1
end

# ---------------------------------------------------------------- D86g
@testset "D86g tangent family against the stand-alone members, the L1 orbit and the closed forms" begin
    run = d86_chain_run()
    plan = run.plan
    rng = StableRNG(d86_seed('g', 1))
    l1raw = Dict{Tuple, Any}()
    d86_each(run; methods = ["L1"]) do it, orbit, raw, c, file
        l1raw[(orbit, ex_label(it.precision), it.tau, it.coordinate_mode)] = (raw, it)
    end
    checked = Ref(0)
    bad = Dict{String, Int}()
    miss!(tag, ok) = ok || (bad[tag] = get(bad, tag, 0) + 1)
    tgids = Set("t$i" for i in 1:7)
    d86_each(run; methods = collect(tgids)) do it, orbit, raw, c, file
        m = D86_TG_CHAINS[parse(Int, it.method[2:end])]
        ser = d86_ds(raw, orbit, it.method, "series")
        cols = d86_at(raw, orbit, it.method, "series", "columns")
        t = raw["orbits/$orbit/t"]
        grid = ex_ref_grid(it.horizon, it.tau; ppd = 5, first_steps = 10, fracs = [0.1, 0.5])
        miss!("grid", t == grid)
        miss!("columns", cols == d86_col_names(m) && size(ser) == (length(t), length(cols)) && ser isa Matrix{Float64})
        miss!("finite", all(isfinite, ser))
        "FLI" in m && miss!("fli nondecreasing", all(j -> issorted(ser[:, j]), findall(c -> startswith(c, "FLI."), cols)))
        if "SALI" in m
            j = findfirst(==("SALI.sali"), cols)
            miss!("sali range", all(x -> 0 <= x <= sqrt(2) * (1 + 2 * EX_EPS), ser[:, j]))
        end
        # stand-alone members at 5 seeded grid steps and at the last step, bitwise
        idx = sort(unique([rand(rng, 1:length(t), 5); length(t)]))
        for mm in m, i in idx
            fin = d86_stand_alone(it, mm, t[i])
            for (name, v) in d86_member_cols(mm, fin)
                j = findfirst(==(name), cols)
                checked[] += 1
                miss!("stand-alone $mm", j !== nothing && reinterpret(UInt64, ser[i, j]) == reinterpret(UInt64, v))
            end
        end
        # the orbit and the monitor columns equal those of the L1 item (I21), bitwise
        key = (orbit, ex_label(it.precision), it.tau, it.coordinate_mode)
        if haskey(l1raw, key)
            raw1, it1 = l1raw[key]
            for name in ("t", "state", "diagnostics/det", "diagnostics/symp")
                miss!("L1 $name", d86_same(raw["orbits/$orbit/$name"], raw1["orbits/$orbit/$name"]))
            end
            it.coordinate_mode == "wrapped" && miss!("L1 winding", d86_same(raw["orbits/$orbit/winding"], raw1["orbits/$orbit/winding"]))
        else
            miss!("no L1 twin", false)
        end
        miss!("no lyap", !haskey(raw, "orbits/$orbit/lyap/L1") && !any(k -> startswith(k, "orbits/$orbit/lyap"), keys(raw)))
    end
    @test isempty(bad)
    @test checked[] > 1000
    @test Set(ex_label(it.precision) for it in plan.items if it.method in tgids) == Set(["Float64", "Double64"])
end

@testset "D86g closed forms through the recorder: the K = 0 shear and Y5 (O18)" begin
    # K = 0 shear with the identity frame (direction 0): tangent [[1, 1], [0, 1]] at every point, w_1 = (1, 0) stays, w_2(n) = (n, 1) (sali A10, fli Q9a, megno M10a).
    orbs = ex_panel_orbits(:std)
    for o in orbs
        o["control"] && continue
        o["K_text"] = "0"; o["K"] = 0.0
    end
    specs = [d86_spec("tangent"; id = "t4", members = ["FLI", "MEGNO", "SALI"])]
    blk = d86_one(specs...; orbits = "S_L", direction = ["0"], horizon = 100)
    blk["arm"][1]["direction"] = ["0"]
    pp = ex_write_panel(ex_path("d86-g-shear-panel.toml"); orbits = orbs)
    pr = ex_write_prereg_run(ex_path("d86-g-shear-prereg.md"), blk)
    fx = (; run = blk, panel = pp, prereg = pr, tag = "g-shear", orbits = orbs)
    store = ex_store("d86-g-shear")
    res = ex_campaign(fx; store = store)
    @test res.status == "complete"
    raw = d86_raw(first(ex_chunk_files(store)))
    ser = d86_ds(raw, 1, "t4", "series"); cols = d86_at(raw, 1, "t4", "series", "columns"); t = raw["orbits/1/t"]
    col(name) = ser[:, findfirst(==(name), cols)]
    # tolerance 0 where the spec says 0: FLI_1 = 0 (S_1 = ln 1) and MEGNO column 1 (all three series)
    @test all(iszero, col("FLI.fli1")) && all(iszero, col("MEGNO.Y1")) && all(iszero, col("MEGNO.Ybar1")) && all(iszero, col("MEGNO.Ybar20_1"))
    # FLI_B = (1/2) ln(n^2 + 1) and SALI = sqrt(2 - 2 n / sqrt(n^2 + 1)) in BigFloat(256), within B_S(n) = 8 eps (n + G_first), G_first = 3 (squared Frobenius norm of the first window [[1,1],[0,1]] is 3)
    bs(n) = 8 * EX_EPS * (n + 3)
    for (i, n) in enumerate(t)
        fli = setprecision(() -> Float64(log(BigFloat(n)^2 + 1) / 2), BigFloat, 256)
        @test abs(col("FLI.fli")[i] - fli) <= bs(n)
        sal(j) = setprecision(() -> sqrt(2 - 2 * BigFloat(j) / sqrt(BigFloat(j)^2 + 1)), BigFloat, 256)
        # B_SALI(n) = 8 eps [Lambda_n + SALI_ref(n) (Psi1 + Psi2)] + B_meas (sali A15), Lambda_n = sum_j SALI_ref(n) / SALI_ref(j); for the linear case 8 eps Psi = B_S(n);
        # B_meas is taken as the SALI floor 4 eps (the stage 04 value is not known here, see the plan)
        lam = setprecision(() -> Float64(sum(sal(n) / sal(j) for j in 1:n)), BigFloat, 256)
        @test abs(col("SALI.sali")[i] - Float64(sal(n))) <= 8 * EX_EPS * lam + Float64(sal(n)) * 2 * bs(n) + 4 * EX_EPS
    end
    # MEGNO column 2 at tau = 1: megno-spec M10(a): "Y(n) = ln(n^2 + 1) - (1/n) ln prod_{i=0..n-1} (i^2 + 1)", within 6 B_S(n)
    for (i, n) in enumerate(t)
        y = setprecision(() -> Float64(log(BigFloat(n)^2 + 1) - sum(log(BigFloat(j)^2 + 1) for j in 0:(n - 1)) / n), BigFloat, 256)
        @test abs(col("MEGNO.Y2")[i] - y) <= 6 * bs(n)
    end
end

# ---------------------------------------------------------------- D86h
@testset "D86h pair: columns, separation, attributes, closed forms" begin
    run = d86_chain_run()
    plan = run.plan
    rng = StableRNG(d86_seed('h', 1))
    bad = Dict{String, Int}()
    miss!(tag, ok) = ok || (bad[tag] = get(bad, tag, 0) + 1)
    xi = Dict{Tuple, Any}()
    pairs = Set(["p1", "p2", "p3", "p4", "p5"])
    n_items = Ref(0)
    y5rli = Dict{Tuple, Any}()
    l1raw = Dict{Tuple, Any}()
    d86_each(run; methods = ["L1"]) do it, orbit, raw, c, file
        l1raw[(orbit, it.tau, it.coordinate_mode)] = raw
    end
    d86_each(run; methods = collect(pairs)) do it, orbit, raw, c, file
        n_items[] += 1
        spec = d86_spec_of(plan, it)
        members = spec["members"]
        ser = d86_ds(raw, orbit, it.method, "series"); cols = d86_at(raw, orbit, it.method, "series", "columns"); t = raw["orbits/$orbit/t"]
        a(name) = raw["orbits/$orbit@$name"]
        miss!("columns", cols == d86_col_names(members))
        xi0 = a("xi0")
        xi[(it.method, orbit, it.coordinate_mode)] = xi0
        miss!("xi0 unit", length(xi0) == 2 && abs(hypot(xi0...) - 1) <= 4 * EX_EPS)
        miss!("attrs", a("Dx") == spec["separation"] && a("mode") == it.coordinate_mode && haskey(raw, "@frame_source") && raw["@frame_source"] == spec["frame_source"] &&
                       haskey(raw, "@direction_rule") && raw["@direction_rule"] == "u = xi_0")
        if spec["frame_source"] == "direction"
            W = EX.deviation_frame(Float64, it.direction, 2)
            miss!("direction xi0", xi0 == [W[1, 1], W[2, 1]])
        end
        # stand-alone pair run (O16) at 5 seeded grid steps and the last one; the replacing form of each member
        for mm in members, i in sort(unique([rand(rng, 1:length(t), 5); length(t)]))
            r = d86_pair_ref(it, spec, xi0, t[i], mm)
            for (name, v) in d86_member_cols(mm, r.fin)
                j = findfirst(==(name), cols)
                miss!("stand-alone $mm", j !== nothing && reinterpret(UInt64, ser[i, j]) == reinterpret(UInt64, v))
            end
            i == 1 && miss!("Dx_r", a("Dx_r") == Float64(r.start.Dx_r))
        end
        # Dx_r against the BigFloat(256) difference of the two initial points (lifted mode: B = A + Dx u in the element type)
        if it.coordinate_mode == "lifted"
            A = ex_sys_x(it, Float64)[2]
            B = A .+ spec["separation"] .* xi0
            dr = setprecision(() -> sqrt(sum((BigFloat(B[k]) - BigFloat(A[k]))^2 for k in 1:2)), BigFloat, 256)
            miss!("Dx_r bigfloat", abs(BigFloat(a("Dx_r")) - dr) <= BigFloat(EX_EPS) * dr)
        end
        # the orbit point is the L1 item's
        raw1 = get(l1raw, (orbit, it.tau, it.coordinate_mode), nothing)
        raw1 === nothing ? miss!("no L1 twin", false) : miss!("orbit = L1", d86_same(raw["orbits/$orbit/state"], raw1["orbits/$orbit/state"]) && d86_same(t, raw1["orbits/$orbit/t"]))
        # RLI = 0 exactly at every grid step for Y5 (all modes) and, below, for the shear
        if it.orbit == 15 && "RLI" in members
            j = findfirst(==("RLI.RLI"), cols)
            miss!("RLI Y5", all(iszero, ser[:, j]))
            y5rli[(it.method, it.coordinate_mode)] = ser[:, j]
        end
    end
    @test n_items[] == 5 * 9 * 2
    @test isempty(bad)
    # the seeded and the direction frame sources give different recorded xi0 and the same orbit
    diff_xi = count(k -> xi[("p1", k[2], k[3])] != xi[("p4", k[2], k[3])], [k for k in keys(xi) if k[1] == "p1" && k[2] != 15])
    @test diff_xi > 0
    # Y5 (no torus): RLI is 0 in lifted and in wrapped mode, bitwise equal (rli criterion 13(a))
    for mth in ("p1", "p2", "p3", "p4", "p5")
        @test d86_same(y5rli[(mth, "lifted")], y5rli[(mth, "wrapped")])
    end
end

@testset "D86h pair on the K = 0 shear: RLI = 0 exactly, wrapped and lifted bitwise equal" begin
    orbs = ex_panel_orbits(:std)
    for o in orbs
        o["control"] && continue
        o["K_text"] = "0"; o["K"] = 0.0
    end
    s = d86_spec("pair"; members = ["RLI"])
    blk = d86_one(s; orbits = "S_L", coordinate_mode = nothing)
    blk["arm"][1]["coordinate_mode"] = ["lifted", "wrapped"]
    pp = ex_write_panel(ex_path("d86-h-shear-panel.toml"); orbits = orbs)
    pr = ex_write_prereg_run(ex_path("d86-h-shear-prereg.md"), blk)
    fx = (; run = blk, panel = pp, prereg = pr, tag = "h-shear", orbits = orbs)
    store = ex_store("d86-h-shear")
    res = ex_campaign(fx; store = store)
    @test res.status == "complete"
    cols = Dict{String, Any}()
    for f in d86_files(store)
        raw = d86_raw(f)
        mode = raw["orbits/1@mode"]
        ser = d86_ds(raw, 1, "pr", "series"); names = d86_at(raw, 1, "pr", "series", "columns")
        cols[mode] = ser[:, findfirst(==("RLI.RLI"), names)]
        @test all(iszero, cols[mode]) && all(iszero, ser[:, findfirst(==("RLI.DL"), names)])
    end
    @test Set(keys(cols)) == Set(["lifted", "wrapped"]) && d86_same(cols["lifted"], cols["wrapped"])
end
