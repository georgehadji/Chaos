# Included by verify/experiments-indicators.jl (section MEM): bytes(orbit) of E28 and the launch capacity of R17, for the reduced plan and for the worked examples of R17 item 4
# (plan only, nothing is run). The reference formula is written here from E28 with the reducer sizes read from the verified reducers (Tangent.reducer_bytes).
const D86_E24 = Dict("SALI" => 2, "FLI" => 2, "MEGNO" => 6, "RLI" => 2)
function ref_bytes(plan, it)
    spec = plan.block.raw["indicator"][findfirst(s -> s["id"] == it.method, plan.block.raw["indicator"])]
    kind = spec["kind"]
    T = it.precision == "Double64" ? Chaoslab.Double64 : Float64
    wrapped = it.coordinate_mode == "wrapped"
    s = (wrapped && it.system != "linear") ? 4 : 2
    pair = kind == "pair"
    d, k, sK = pair ? 4 : 2, it.vectors, pair ? 2s : s
    M = Tg.bytes_per_orbit(Val(d), Val(k), T, sK)
    rb(r) = Tg.reducer_bytes(r, T, Val(d), Val(k))
    mon = 3 * sizeof(T)
    mem = Dict("SALI" => Chaoslab.SALI(), "FLI" => Chaoslab.FLI(replace = false), "MEGNO" => Chaoslab.MEGNO(replace = false), "RLI" => Chaoslab.RLI(replace = false))
    R, cols = if kind in ("tangent", "pair")
        sum(rb(mem[m]) for m in spec["members"]) + mon, sum(D86_E24[m] for m in spec["members"])
    elseif kind == "stretching"
        sys = Chaoslab.StandardMap(1.0)
        rb(Chaoslab.StretchingNumbers(; range = Chaoslab.stretching_range(sys), bins = spec["bins"], interval = spec["interval"], min_per_cell = spec["min_per_cell"])) + mon, 1
    elseif kind == "frequency"
        mon, 0
    elseif kind == "phase_advance"
        rb(Chaoslab.WeightedPhaseAdvance(window = Tuple(spec["windows"]), observable = Chaoslab.PlanarObservable(1, 2), birkhoff_order = spec["birkhoff_order"])) + mon, 0
    else
        obs = map(spec["observables"]) do ob
            c = ob["coordinate"]
            ob["kind"] == "cos" ? Chaoslab.CosObservable(c, ob["frequency"]) : ob["kind"] == "sin" ? Chaoslab.SinObservable(c, ob["frequency"]) : Chaoslab.LiftedObservable(c, c + 2)
        end
        rb(Chaoslab.WeightedBirkhoff(window = Tuple(spec["windows"]), observables = Tuple(obs), birkhoff_order = spec["birkhoff_order"])) + mon, 0
    end
    gr = plan.block.raw["grid"]
    ng = length(EX.record_grid(it.horizon, it.tau; points_per_decade = gr["points_per_decade"], first_steps = gr["first_steps"], extra_fractions = gr["extra_fractions"]))
    c = 1 + 2 + (wrapped ? 2 : 0) + 2 + cols
    S = 8 * c * ng
    B = kind == "frequency" ? (T === Float64 ? 32 : 64) * maximum(spec["windows"]) : 0
    return (; M, R, S, B, total = M + R + S + B, ng, c)
end
ref_capacity(plan, it, free, f) = Int(fld(big(free) * big(f), big(ref_bytes(plan, it).total)))

# ---- the reduced plan: every spec, both precisions, wrapped
let
    nbad = n = 0
    for id in sort(collect(keys(SPECS))), prec in ("Float64", "Double64"), mode in ("wrapped", "lifted")
        its = [it for it in PLAN.items if it.method == id && it.precision == prec && it.coordinate_mode == mode && it.system == "standard"]
        isempty(its) && continue
        it = first(its)
        got, want = EX.orbit_bytes(PLAN, it), ref_bytes(PLAN, it)
        n += 1
        same = (got.M, got.R, got.S, got.B, got.total) == (want.M, want.R, want.S, want.B, want.total)
        nbad += same ? 0 : 1
        info("D86n", "bytes(orbit) of $id $prec $mode: module M $(got.M) R $(got.R) S $(got.S) B $(got.B) total $(got.total); E28 written here: M $(want.M) R $(want.R) S $(want.S) B $(want.B) total $(want.total); $(want.ng) grid steps, $(want.c) columns")
    end
    emit("D86n", "bytes(orbit) of every (spec, precision, mode) of the reduced plan: module against E28 with the reducer bytes read from the reducers; combinations compared / differing", "$n / $nbad", "0 differ", pf(nbad == 0 && n >= 10))
    sizes = Dict("SALI" => 72, "FLI" => 96, "MEGNO" => 200, "RLI" => 88)
    got = Dict(nm => Tg.reducer_bytes(r, Float64, Val(2), Val(2)) for (nm, r) in (("SALI", Chaoslab.SALI()), ("FLI", Chaoslab.FLI(replace = false)), ("MEGNO", Chaoslab.MEGNO(replace = false)), ("RLI", Chaoslab.RLI(replace = false))))
    got["RLI"] = Tg.reducer_bytes(Chaoslab.RLI(replace = false), Float64, Val(4), Val(2))
    emit("D86n", "reducer state bytes at Float64 (SALI, FLI, MEGNO, RLI) read from the reducers", "$(got["SALI"]), $(got["FLI"]), $(got["MEGNO"]), $(got["RLI"])", "72, 96, 200, 88", pf(all(got[k] == v for (k, v) in sizes)))
    st1000 = Chaoslab.StretchingNumbers(; range = (-1.0, 1.0), bins = 1000, interval = 4, min_per_cell = 2)
    ws = Tg.reducer_workspace(st1000, Float64, Val(2), Val(2))
    emit("D86n", "stretching reducer at 1000 bins: reducer_bytes (state 104 plus the workspace of 2 (bins + 2) integers of 4 bytes = $(4 * ws))", "$(Tg.reducer_bytes(st1000, Float64, Val(2), Val(2)))", "104 + 8016 = 8120", pf(Tg.reducer_bytes(st1000, Float64, Val(2), Val(2)) == 104 + 8016 && 4 * ws == 8016))
end

# ---- the worked examples of R17 item 4 (plan only)
let
    run = block_run()
    sp(id) = first(s for s in run["indicator"] if s["id"] == id)
    mk(id, base; kw...) = (s = deepcopy(sp(base)); s["id"] = id; for (k, v) in kw; s[string(k)] = v; end; s)
    specs = [mk("wtg", "tg"), mk("wfq", "fq"; windows = [250, 1000, 4000]), mk("wst", "st"; bins = 1000, counts = "none"), mk("wpa", "pa"; windows = [250, 1000, 4000]),
             mk("wbk", "bk"; windows = [250, 1000, 4000])]
    arm(name, ind, hz, prec = ["Float64"]) = Dict{String, Any}("name" => name, "orbits" => "S_L", "indicators" => [ind], "precisions" => prec, "horizons" => [hz], "coordinate_mode" => ["wrapped"])
    arms = [arm("wtg", "wtg", 4000), arm("wfq", "wfq", 8000), arm("wfqd", "wfq", 8000, ["Double64"]), arm("wst", "wst", 4000), arm("wpa", "wpa", 8001), arm("wbk", "wbk", 8000)]
    r2 = deepcopy(run); r2["indicator"] = specs; r2["arm"] = arms
    pre = make_prereg(joinpath(SCR, "prereg-worked.md"), r2)
    pl = EX.plan_run(pre, PP.panel)
    first_item(id, prec = "Float64") = first(it for it in pl.items if it.method == id && it.precision == prec && it.system == "standard")
    F16, F4 = 16 * 10^9, 4 * 10^9
    cases = [("(a) tangent FLI, MEGNO, SALI, N = 4000", "wtg", "Float64", 18864, 424088, 106022),
             ("(b) frequency, windows 250, 1000, 4000 (run length 8000)", "wfq", "Float64", 136000, 58823, nothing),
             ("(c) stretching, 1000 bins, N = 4000", "wst", "Float64", 16896, 473484, nothing),
             ("(d) phase advance, three windows (run length 8001)", "wpa", "Float64", 8224, 972762, nothing),
             ("(d) weighted Birkhoff, K = 4, H = 3", "wbk", "Float64", 8544, 936329, nothing)]
    for (name, id, prec, tot, cap, cap4) in cases
        it = first_item(id, prec)
        b = EX.orbit_bytes(pl, it)
        c = EX.group_capacity(pl, it; free_bytes = F16, f = 0.5)
        ok = b.total == tot && c == cap && c == ref_capacity(pl, it, F16, 0.5)
        if cap4 !== nothing
            c4 = EX.group_capacity(pl, it; free_bytes = F4, f = 0.5)
            ok &= c4 == cap4
        end
        emit("D86n/R17", "worked example $name: bytes(orbit) M $(b.M) + R $(b.R) + S $(b.S) + B $(b.B) = $(b.total); capacity at F = 16 GB, f = 0.5: $c" * (cap4 === nothing ? "" : "; at F = 4 GB: $(EX.group_capacity(pl, it; free_bytes = F4, f = 0.5))"),
             "$(b.total), $c", "$tot, $cap" * (cap4 === nothing ? "" : ", $cap4 at 4 GB"), pf(ok))
    end
    bd = EX.orbit_bytes(pl, first_item("wfq", "Double64"))
    emit("D86n/R17", "frequency, Double64: the buffer term B for W_max = 4000", bd.B, "256000 (64 W_max)", pf(bd.B == 256000))
    cf = EX.group_capacity(pl, first_item("wfq"); free_bytes = div(F16, 2), f = 0.5)
    emit("D86j(a)", "frequency capacity halves (floor) when the free memory halves: capacity at 16 GB, 8 GB", "$(EX.group_capacity(pl, first_item("wfq"); free_bytes = F16, f = 0.5)), $cf", "58823, 29411", pf(EX.group_capacity(pl, first_item("wfq"); free_bytes = F16, f = 0.5) == 58823 && cf == 29411))
    ch = [c.key for c in pl.chunks]
    emit("D86n", "chunk keys and hashes of the worked-example plan do not depend on the free memory (plan is a function of its inputs only)", "plan built once; keys contain no memory field", "true", pf(!any(k -> occursin("free", k) || occursin("capacity", k), ch)))
end
