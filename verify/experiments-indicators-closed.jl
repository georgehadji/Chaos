# Included by verify/experiments-indicators.jl (section O18): the exactly known cases the indicator specs gate, now through the recorder (spec O18, D86g, D86h, D86k, D86l).
# A small campaign of its own (scratch folder): the K_std = 0 shear orbits of the real panel (cell K = 0, 12 orbits) and the Y5 control, frame direction 0 (the identity frame), horizon 400;
# frequency and phase advance with windows 64 and 128. Closed forms and bounds are those of the criteria the test comments name (sali A10 and A15, fli-ofli Q9a, megno M10a,
# rli RL9a, naff-fma criteria 7(b) and 8 with the lifted-coordinate term of D87); the recorder adds no tolerance.
let
    dir = joinpath(SCR, "o18"); mkpath(dir)
    t = TOML.parsefile(PANEL)
    k0 = sort([o for o in t["orbit"] if !o["control"] && o["k_index"] == 1]; by = o -> (o["class"], o["draw_index"]))
    pick = vcat([o for o in k0 if o["class"] == "T"][1:6], [o for o in k0 if o["class"] == "R"][1:6])
    ctl = [o for o in t["orbit"] if o["control"] && o["id"] in (871, 872, 883)]
    sel = sort(vcat(pick, ctl); by = o -> o["id"])
    for o in sel
        o["in_SL"] = !o["control"]
    end
    t2 = copy(t); t2["orbit"] = sel
    panel = joinpath(dir, "panel-o18.toml")
    open(io -> TOML.print(io, t2; sorted = true), panel, "w")
    run = block_run()
    sp(id) = first(s for s in run["indicator"] if s["id"] == id)
    fq0 = deepcopy(sp("fq")); fq0["id"] = "fq0"; fq0["windows"] = [64, 128]; fq0["mean_removal"] = false      # the closed form of naff-fma criterion 7(b) is for a window without mean removal (the removal moves the peak of a pure line by about 1e-9)
    pal = deepcopy(sp("pa")); pal["id"] = "pal"; pal["windows"] = [64, 128]; pal["observable"] = "lifted"
    run["indicator"] = [deepcopy(sp("tg")), deepcopy(sp("pr")), fq0, pal]
    A(name, orbits, ind, hz; kw...) = merge(Dict{String, Any}("name" => name, "orbits" => orbits, "indicators" => ind, "precisions" => ["Float64"], "horizons" => [hz]), Dict{String, Any}(string(k) => v for (k, v) in kw))
    run["arm"] = [A("s0w", "P", ["tg", "pr"], 400; direction = ["0"], coordinate_mode = ["wrapped"]), A("s0l", "P", ["tg"], 400; direction = ["0"]),
                  A("y5long", "P", ["pr"], 4000; coordinate_mode = ["wrapped"]),
                  A("f0", "S_L", ["fq0"], 256; coordinate_mode = ["wrapped"]), A("pw", "S_L", ["pal"], 257; coordinate_mode = ["wrapped"]), A("pl", "S_L", ["pal"], 257)]
    run["health"]["stop_failed_fraction"] = 0.5      # 18 chunks only: one failed chunk would be 5.6 percent, over the real stop value 0.05
    pre = make_prereg(joinpath(dir, "prereg-o18.md"), run)
    store, mf = joinpath(dir, "store"), joinpath(dir, "manifest.json")
    rm(store; recursive = true, force = true); rm(mf; force = true)
    plan = EX.plan_run(pre, panel)
    r = EX.run_campaign(pre, panel; store = store, manifest = mf, threads = NT, thread_override = NT != K19, allow_dirty = true, capacity = 1000)
    ylong = [c.ordinal for c in plan.chunks if c.run_type.method == "pr" && c.run_type.horizon == 4000 && c.cell == "controls"]
    emit("O18", "closed-form campaign (12 K_std = 0 shear orbits, controls 871, 872, Y5 883; frame direction 0): status, chunks, chunks failed twice (only the controls chunk of the pair arm y5long at N = 4000, finding F2)", "$(r.status), $(length(plan.chunks)) chunks, failed $(r.failed)", "complete with failed chunks, failed = $ylong", pf(r.status == "complete with failed chunks" && sort(r.failed) == sort(ylong)))
    R = load_store(store)
    let ch = R[only(ylong)]
        info("F2", "pair at N = 4000 on the controls cell (holds the linear control Y5, orbit 883): chunk status $(ch.attrs["status"]), attempt $(ch.attrs["attempt"]); failure_reasons $(ch.attrs["failure_reasons"]) (R18 item 5: a failed pair is deterministic and not retried to success; the study chooses the horizon)")
        emit("F2", "the failing orbit of that chunk is the Y5 control (orbit 883) and the failure is the overflow of the shadow orbit at step 770", occursin("orbit 883", only(ch.attrs["failure_reasons"])) && occursin("step 770", only(ch.attrs["failure_reasons"])), "true (design limit of the pair on Y5, not a defect)", pf(occursin("orbit 883", only(ch.attrs["failure_reasons"]))))
    end
    eps64 = eps(Float64)
    big256(f) = setprecision(f, BigFloat, 256)
    item_of = Dict(it.hash => it for it in plan.items)
    bs(n) = 8 * eps64 * (n + 3)
    zero_bad = fli_bad = sali_bad = meg_bad = 0
    nzero = nfli = nsali = nmeg = 0
    worst = Dict("fli" => 0.0, "sali" => 0.0, "megno" => 0.0)
    for (ord, ch) in R, o in ch.orbits
        c = plan.chunks[ord]
        c.run_type.method == "tg" || continue
        it = item_of[o.attrs["item_hash"]]
        it.param == "0" || continue
        tab = o.indicators["tg"]
        for nm in ("FLI.fli1", "MEGNO.Y1", "MEGNO.Ybar1", "MEGNO.Ybar20_1")
            nzero += length(o.t)
            zero_bad += count(!iszero, colvec(tab, nm))
        end
        for (i, n) in enumerate(o.t)
            fli = big256(() -> Float64(log(BigFloat(n)^2 + 1) / 2))
            nfli += 1
            e = abs(colvec(tab, "FLI.fli")[i] - fli)
            worst["fli"] = max(worst["fli"], e / bs(n)); fli_bad += e <= bs(n) ? 0 : 1
            sal(j) = big256(() -> sqrt(2 - 2 * BigFloat(j) / sqrt(BigFloat(j)^2 + 1)))
            lam = big256(() -> Float64(sum(sal(n) / sal(j) for j in 1:n)))
            bound = 8 * eps64 * lam + Float64(sal(n)) * 2 * bs(n) + 4 * eps64
            nsali += 1
            e = abs(colvec(tab, "SALI.sali")[i] - Float64(sal(n)))
            worst["sali"] = max(worst["sali"], e / bound); sali_bad += e <= bound ? 0 : 1
            y = big256(() -> Float64(log(BigFloat(n)^2 + 1) - sum(log(BigFloat(j)^2 + 1) for j in 0:(n - 1)) / n))
            nmeg += 1
            e = abs(colvec(tab, "MEGNO.Y2")[i] - y)
            worst["megno"] = max(worst["megno"], e / (6 * bs(n))); meg_bad += e <= 6 * bs(n) ? 0 : 1
        end
    end
    emit("O18/D86g", "K_std = 0 shear through the recorder (wrapped and lifted, direction 0): FLI.fli1 and MEGNO column 1 (Y1, Ybar1, Ybar20_1) exactly 0 at every grid step; values / nonzero", "$nzero / $zero_bad", "all exactly 0", pf(zero_bad == 0 && nzero > 1000))
    emit("O18/D86g", "FLI_B against (1/2) ln(n^2 + 1) in BigFloat(256) within B_S(n) = 8 eps (n + 3): grid rows / outside / worst ratio", "$nfli / $fli_bad / $(fmt(worst["fli"]))", "0 outside", pf(fli_bad == 0 && nfli > 1000))
    emit("O18/D86g", "SALI against sqrt(2 - 2 n / sqrt(n^2 + 1)) within 8 eps Lambda_n + SALI_ref 2 B_S + 4 eps (sali A15 with the floor as B_meas): grid rows / outside / worst ratio", "$nsali / $sali_bad / $(fmt(worst["sali"]))", "0 outside", pf(sali_bad == 0 && nsali > 1000))
    emit("O18/D86g", "MEGNO column 2 against ln(n^2 + 1) - (1/n) sum_{j=0}^{n-1} ln(j^2 + 1) (megno M10a) within 6 B_S(n): grid rows / outside / worst ratio", "$nmeg / $meg_bad / $(fmt(worst["megno"]))", "0 outside", pf(meg_bad == 0 && nmeg > 1000))
    # RLI exactly zero on the shear and on Y5
    nz = nbadz = ny5 = 0
    for (ord, ch) in R, o in ch.orbits
        plan.chunks[ord].run_type.method == "pr" && plan.chunks[ord].run_type.horizon == 400 || continue
        it = item_of[o.attrs["item_hash"]]
        (it.param == "0" || o.id == 883) || continue
        tab = o.indicators["pr"]
        nz += length(o.t)
        nbadz += count(!iszero, colvec(tab, "RLI.RLI")) + count(!iszero, colvec(tab, "RLI.DL"))
        o.id == 883 && (ny5 += 1)
    end
    emit("O18/D86h", "pair chain RLI, MEGNO: RLI and DL exactly 0 at every grid step for the shear orbits and for Y5 (orbit 883 present: $(ny5 == 1)); grid rows / nonzero values", "$nz / $nbadz", "all exactly 0", pf(nbadz == 0 && ny5 == 1 && nz > 1000))
    # frequency: nu = frac(p tau / 2 pi), p the (constant) momentum read from the stored orbit
    nf = nfb = 0; wf = 0.0
    pi256 = big256(() -> BigFloat(pi))
    circ(a, b) = (d = abs(a - b); min(d, 1 - d))
    for (ord, ch) in R, o in ch.orbits
        plan.chunks[ord].run_type.method == "fq0" || continue
        it = item_of[o.attrs["item_hash"]]
        pst = o.state[1, 2]
        ref = big256(() -> (v = BigFloat(pst) / (2 * pi256); v - floor(v)))
        tab = o.indicators["fq0"]
        for h in 1:2, c in (1, 2)
            nf += 1
            e = Float64(big256(() -> circ(BigFloat(tab.final[h, c]), ref)))
            wf = max(wf, e / (16 * eps64)); nfb += e <= 16 * eps64 ? 0 : 1
        end
    end
    emit("O18/D86k", "frequency on the K_std = 0 line (wrapped, windows 64 and 128): nu1, nu2 within 16 eps of p/(2 pi) mod 1 (BigFloat(256), p the stored momentum; naff-fma criterion 7(b)): values / outside / worst ratio", "$nf / $nfb / $(fmt(wf))", "0 outside", pf(nfb == 0 && nf >= 40))
    # phase advance with the lifted observable: nu = lifted p / (2 pi) within 16 eps + max|q| eps (naff-fma criterion 8 with the D87 term)
    np = npb = 0; wp = 0.0
    for (ord, ch) in R, o in ch.orbits
        plan.chunks[ord].run_type.method == "pal" || continue
        it = item_of[o.attrs["item_hash"]]
        wrapped = it.coordinate_mode == "wrapped"
        pl = big256(() -> BigFloat(o.state[1, 2]) + (wrapped ? 2 * pi256 * BigFloat(o.winding[1, 2]) : BigFloat(0)))
        ref = big256(() -> pl / (2 * pi256))
        qmax = big256(() -> Float64(abs(BigFloat(parse(Float64, it.x0))) + 257 * abs(pl)))
        tab = o.indicators["pal"]
        for h in 1:2, c in (1, 2)
            np += 1
            e = Float64(big256(() -> abs(BigFloat(tab.final[h, c]) - ref)))
            bound = 16 * eps64 + qmax * eps64
            wp = max(wp, e / bound); npb += e <= bound ? 0 : 1
        end
    end
    emit("O18/D86l", "phase advance, lifted observable, wrapped and lifted mode (windows 64 and 128): nu1, nu2 within 16 eps + max|q| eps of the lifted p/(2 pi) (D87 bound): values / outside / worst ratio", "$np / $npb / $(fmt(wp))", "0 outside", pf(npb == 0 && np >= 80))
    info("O18", "stretching, birkhoff and the Y5 rotation block are covered in the suite (D86i, D86m); the Y5 origin start and the K_std = 0 weighted Birkhoff criterion 3 (window 1000) are not repeated here")
end
