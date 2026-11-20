# Amendment D125, stage 02: D125p (backend), D125q (store and resume), D125r (memory), D125s (health and controls), D125t (manifest and pilot).
# Seeds: d125_seed('p' | 'q' | 'r' | 's' | 't', sub). Campaigns are the fixture campaigns of ex_d125_items.jl (tag "std" is shared and cached).

@testset "D125p backend: cpu, gpu, auto (with and without a simulated device) store identical datasets; the reason strings (R41)" begin
    R = Chaoslab.FALLBACK_REASONS
    base = d125_campaign("p-cpu", d125_block_spectra())
    ref = d86_store_sha(base.store; skip = D86_BACKEND_ATTRS)
    for (req, dev) in (("cpu", nothing), ("auto", nothing), ("auto", ExDev(true, true, true)), ("gpu", ExDev(true, true, true)), ("gpu", nothing))
        run = d125_block_spectra(); run["backend"]["Float64"] = req
        tag = "p-$req-$(dev === nothing ? "nodev" : "dev")"
        got = ex_with_device(dev === nothing ? ExFakeDevice(false, false) : dev) do
            d125_campaign(tag, run)
        end
        @test got.res.status == "complete" && d86_store_sha(got.store; skip = D86_BACKEND_ATTRS) == ref
        reasons = Set(String(EX.read_chunk(f).attrs["backend_fallback_reason"]) for f in ex_chunk_files(got.store)
                      if (it = got.plan.items[first(got.plan.chunks[ex_ordinal(f)].items)]; it.system in ("froeschle", "chain", "henon_heiles")))
        @test all(f -> EX.read_chunk(f).attrs["backend"] == "cpu", ex_chunk_files(got.store))
        req == "cpu" ? @test(reasons == Set([""])) : @test(R[6] in reasons)
    end
    @test length(R) == 9                                                                                         # no tenth reason string (the note of R41 is not a reason)
end

@testset "D125q store and resume: killed and resumed with another thread count equals the uninterrupted campaign; stale; layout (R40, I48)" begin
    rng = StableRNG(d125_seed('q'))
    ref = d125_campaign("std")
    M = length(ref.plan.chunks)
    refsha = d86_store_sha(ref.store; skip = D86_BACKEND_ATTRS)
    for (j, thr) in enumerate((1, 3))
        k = 1 + floor(Int, rand(rng) * (M - 1))
        store = ex_store("d125-q-resume-$j")
        @test ex_killed(() -> ex_campaign(ref.fx; store = store, threads = 4, hooks = EX.RunHooks(kill = ex_kill_after(k))))
        res = ex_campaign(ref.fx; store = store, threads = thr)
        @test res.status == "complete" && sort(res.ran) == collect((k + 1):M)
        @test d86_store_sha(store; skip = D86_BACKEND_ATTRS) == refsha                                          # SHA-256 of the canonical read-back, all datasets and attributes
    end
    # a changed panel file makes the files stale and stops the run (D53o(b))
    other = d125_panel_orbits(); other[1]["start"][1] += 0.25; other[1]["start_dec"] = string.(other[1]["start"])
    pp = d125_write_panel(ex_path("d125-q-panel2.toml"); orbits = other)
    store = ex_store("d125-q-stale"); cp(ref.store, store)
    err = try
        EX.run_campaign(ref.fx.prereg, pp; store = store, manifest = store * "-manifest.json", threads = 2, thread_override = true, allow_dirty = true, code_version = EX_CODE,
                        capacity = 64, manifest_reset = true); nothing
    catch e
        e
    end
    @test err isa ArgumentError && occursin("stale", err.msg)
    # the files read back by O35 (the raw reader) hold the names, shapes and dtypes of R40
    seen = Dict{String, Int}()
    d125_each(ref) do it, raw, c, file
        it.system in ("froeschle", "chain", "henon_heiles") || return
        o = it.orbit; pre = "orbits/$o/"; n = length(raw[pre * "t"]); d = 4
        k = it.method == "L1" ? 1 : d
        for (nm, ok) in (("state", size(raw[pre * "state"]) == (n, d) && eltype(raw[pre * "state"]) == Float64),
                         ("lyap", size(raw[pre * "lyap/" * it.method]) == (n, k)),
                         ("frame", raw["@frame"] == it.frame && raw["@cell"] == c.cell),
                         ("unit", raw["@exponent_unit"] == (it.system == "henon_heiles" ? "time" : "step")),
                         ("integrator", it.system != "henon_heiles" || (raw["@integrator"] == it.integrator && haskey(raw, "@step_time"))),
                         ("energy", it.system != "henon_heiles" || (size(raw[pre * "diagnostics/energy"]) == (n,) && isnan(raw[pre * "diagnostics/energy@threshold"]) &&
                                                                    raw[pre * "diagnostics/energy@first_exceedance_step"] == 0)),
                         ("system", raw[pre * "@system"] == it.system && raw[pre * "@ic"] == it.start))
            ok || (seen[nm] = get(seen, nm, 0) + 1)
        end
    end
    @test isempty(seen)
end

@testset "D125r memory: bytes(orbit) of five fixtures, capacity, refusal, results independent of capacity (R42 item 5, E62, I24)" begin
    plan = d125_plan(d125_fix("std"))
    # worked examples by hand (E62): M = 8 (s + d k + 2 k) + 16, S = 8 c n_g, c = 1 + d + k + n_diag (lifted: no winding), n_diag = 3 (det, symp, orth) for k < d, 5 for k = d
    items = Dict(m => first(i for i in plan.items if i.system == "froeschle" && i.method == m) for m in ("L1", "L3", "L5"))
    ch = first(i for i in plan.items if i.system == "chain" && i.method == "L5")
    hh = first(i for i in plan.items if i.system == "henon_heiles" && i.method == "L5")
    ng = length(ex_ref_grid(100, 1; ppd = 5, first_steps = 10, fracs = (0.1, 0.5)))
    ngf = length(ex_ref_grid(160, 1; ppd = 5, first_steps = 10, fracs = (0.1, 0.5)))
    function hand(it, ng, s, d, k, ndiag, ndiag_flow)
        Mb = 8 * (s + d * k + 2k) + 16
        c = 1 + d + k + ndiag + ndiag_flow
        return (M = Mb, S = 8 * c * ng)
    end
    for (it, ng_, s, d, k, nd, nf) in ((items["L1"], ng, 4, 4, 1, 2, 0), (items["L3"], ng, 4, 4, 4, 5, 0), (items["L5"], ng, 4, 4, 4, 5, 0), (ch, ng, 4, 4, 4, 5, 0), (hh, ngf, 4, 4, 4, 5, 1))
        b = EX.orbit_bytes(plan, it); h = hand(it, ng_, s, d, k, nd, nf)
        @test b.M == h.M && b.S == h.S && b.total == b.M + b.R + b.S + b.B && b.ng == ng_                              # to the byte; R (the reducer bytes) from the library
    end
    # capacity = floor(f free / bytes) for three free-memory values
    it = items["L5"]; bpo = EX.orbit_bytes(plan, it).total
    for free in (10^6, 3 * 10^7, 10^9)
        @test EX.group_capacity(plan, it; free_bytes = free, f = 0.5) == floor(Int, 0.5 * free / bpo)
    end
    # a start with one orbit that does not fit is refused naming kind, bytes and free memory
    free = bpo - 1
    @test ex_rejects(() -> EX.group_capacity(plan, it; free_bytes = free, f = 1.0), ["bytes", string(bpo), string(free)])
    # results do not depend on the capacity (I24)
    ref = d125_campaign("std")
    small = d125_campaign("r-cap3", d125_block_spectra(); capacity = 3)
    @test d86_store_sha(small.store; skip = D86_BACKEND_ATTRS) == d86_store_sha(ref.store; skip = D86_BACKEND_ATTRS)
end

@testset "D125s health and controls: non-finite values, control bound for three control kinds, zero exponent, rerun rule (R42)" begin
    ref = d125_campaign("std")
    c5 = first(c for c in ref.plan.chunks if ref.plan.items[first(c.items)].system == "chain" && ref.plan.items[first(c.items)].method == "L5")
    # a non-finite value fails its chunk naming the orbit (both attempts); in the first attempt only the rerun passes
    hit = Ref(0)
    inf!(orbits, info) = (hit[] = orbits[1].id; orbits[1].state[2, 1] = Inf)
    store = ex_store("d125-s-inf")
    @test ex_killed(() -> ex_campaign(ref.fx; store = store, hooks = EX.RunHooks(mutate = ex_fault(inf!; ord = c5.ordinal, attempts = (1, 2)), kill = ex_kill_after(c5.ordinal))))
    ch = EX.read_chunk(joinpath(store, "chunk-$(lpad(c5.ordinal, 5, '0')).h5"))
    @test ch.attrs["status"] == "failed" && any(r -> occursin("non-finite", r) && occursin("orbit $(hit[])", r), ch.attrs["failure_reasons"])
    nan!(orbits, info) = (hit[] = orbits[1].id; first(values(orbits[1].lyap))[2, 1] = NaN)
    store = ex_store("d125-s-nan")
    @test ex_killed(() -> ex_campaign(ref.fx; store = store, hooks = EX.RunHooks(mutate = ex_fault(nan!; ord = c5.ordinal), kill = ex_kill_after(c5.ordinal))))
    a = EX.read_chunk(joinpath(store, "chunk-$(lpad(c5.ordinal, 5, '0')).h5")).attrs
    @test a["status"] == "ok" && a["attempt"] == 2 && any(r -> occursin("non-finite", r) && occursin("orbit $(hit[])", r), a["first_attempt_reasons"])
    # control rule: the three control kinds (y5 at d = 4, y5 at d = 6, the chain origin), k = d; the exact list passes, one exponent shifted by twice control_bound stops
    kinds = Any[("y5-4", "linear", "3;5/2", 4, d125_exact_y5(["3", "5/2"])),
                ("y5-6", "linear", "2;3/2;4", 6, d125_exact_y5(["2", "3/2", "4"])),
                ("chain", "chain", "2;1.5;0;open", 4, d125_exact_chain_origin(2, "1.5"))]
    function ctl_orbits(exact_of)
        orbs = Dict{String, Any}[]
        for (nm, sys, par, d, ex) in kinds
            id = length(orbs) + 1
            push!(orbs, Dict{String, Any}("id" => id, "cell" => "controls", "class" => "control", "system" => sys, "param" => par, "start" => zeros(d),
                  "start_dec" => fill("0", d), "control" => true, "control_kind" => startswith(nm, "y5") ? "y5" : "chain_origin", "dimension" => d,
                  "exact_exponents" => exact_of(nm, ex), "in_D" => false, "run_seed" => string(EX_RUN_SEED), "orbit_seed" => string(EX.orbit_seed(EX_RUN_SEED, id))))
        end
        return orbs
    end
    arm = Any[d125_arm("ctl", "controls", ["L1", "L5"], [160]; vectors = "d", frames = ["a"], h1_check = true)]
    bound = 0.05
    ok = d125_campaign("s-ctl-ok", d125_run(arm); orbits = ctl_orbits((nm, ex) -> ex))
    @test ok.res.status == "complete"
    for (j, (nm, _, _, d, _)) in enumerate(kinds), idx in unique([1, d])
        shifted = d125_campaign("s-ctl-$nm-$idx", d125_run(arm); orbits = ctl_orbits((n2, ex) -> n2 == nm ? (v = copy(ex); v[idx] += 2 * bound; v) : ex))
        @test shifted.res.status == "stopped: control bound" && any(r -> occursin("control bound", r) && occursin("orbit $j", r), shifted.res.reasons)
    end
    # a zero exact exponent is checked: the chain origin with K = 0 has exact exponents 0; shifted by twice the bound it stops, at zero it passes
    zk = ctl_orbits((nm, ex) -> ex)[3]; zk["param"] = "2;0;0;open"; zk["exact_exponents"] = zeros(4)
    zr = d125_campaign("s-ctl-zero", d125_run(arm); orbits = [zk])
    @test zr.res.status == "complete"
    zk2 = copy(zk); zk2["exact_exponents"] = [2 * bound, 0.0, 0.0, 0.0]
    @test d125_campaign("s-ctl-zero2", d125_run(arm); orbits = [zk2]).res.status == "stopped: control bound"
end

@testset "D125t manifest and pilot (R43)" begin
    ref = d125_campaign("std")
    man = EX.read_json(ref.res.manifest)
    for k in ("systems", "integrators", "frames", "exponent_unit")
        @test haskey(man, k)
    end
    @test Set(keys(man["systems"])) == Set(["y3a", "ch2", "hh1"]) && man["systems"]["hh1"]["kind"] == "henon_heiles" && man["systems"]["hh1"]["d"] == 4
    @test haskey(man["systems"]["hh1"], "integrator") && !haskey(man["systems"]["y3a"], "integrator")
    @test Set(i["id"] for i in man["integrators"]) == Set(["lf", "kl"]) && all(c -> haskey(c, "cell") && haskey(c, "frame") && haskey(c, "integrator"), man["chunks"])
    # the manifest of a plan without a new item has the keys it had (I40)
    old = ex_small_run("d125-t-old"); mo = EX.read_json(old.res.manifest)
    @test !any(k -> haskey(mo, k), ("systems", "integrators", "estimators", "frames", "exponent_unit"))
    # pilot: starts from the cell sampler on the pilot stream, none within 1e-6 of a panel start, the cost per orbit-step for every kind
    out = ex_path("d125-t-pilot"); rm(out; recursive = true, force = true)
    res = EX.run_pilot(ref.fx.prereg, ref.fx.panel; out = out, threads = 1, thread_override = true, allow_dirty = true, code_version = EX_CODE, capacity = 64,
                       hooks = EX.RunHooks(clock = ex_clock(info -> info.steps * 150e-9; start = 0.0)))
    rep = res.report
    label = ref.fx.run["pilot"]["seed_label"]; per = ref.fx.run["pilot"]["orbits_per_k"]                                          # reading U11: the same key serves the cells
    panel_starts = Dict(c => [o["start"] for o in ref.fx.orbits if o["cell"] == c] for c in ("y3a", "ch2", "hh1"))
    for cell in ("y3a", "ch2", "hh1")
        got = [Float64.(parse.(Float64, s["start"])) for s in rep["starts"] if s["cell"] == cell]
        want = cell == "hh1" ? first(d125_flow_starts(0, cell, "1/12", ("-0.3", "0.4", "-0.25", "0.25"), per; text = d125_pilot_text(label, cell))) :
                               d125_points(0, cell, 4, per; text = d125_pilot_text(label, cell))
        @test got == want
        @test all(g -> all(p -> maximum(abs.(g .- p)) > 1e-6, panel_starts[cell]), got)
    end
    kinds = Set(c["system_kind"] for c in rep["cost"] if haskey(c, "system_kind"))
    @test kinds == Set(["froeschle", "chain", "henon_heiles"]) && all(c -> haskey(c, "thread_seconds_per_orbit_step"), [c for c in rep["cost"] if haskey(c, "system_kind")])
    # CPU reference rerun is bitwise: a second campaign into another store equals the first on every dataset and attribute (all new kinds run on the CPU)
    again = d125_campaign("t-again", d125_block_spectra())
    @test d86_store_sha(again.store; skip = D86_BACKEND_ATTRS) == d86_store_sha(ref.store; skip = D86_BACKEND_ATTRS)
end
