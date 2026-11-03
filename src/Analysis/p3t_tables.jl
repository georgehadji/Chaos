# Study p3-abstention (Amendment D95): the tables file `p3-abstention-tables.md` (E71, E73): the header, the decisions, T1 to T11, the data blocks of F2 and F3 with digests, the readings, the
# deviations and the exploratory list. No clock, no absolute path, no thread count.

_c(x::Rational) = denominator(x) == 1 ? string(numerator(x)) : "$(numerator(x))/$(denominator(x))"
_c(x::Nothing) = "n.a."
_c(x::Symbol) = replace(string(x), "_" => " ")
_c(x) = _fmt(x)
mdt(cols, rows) = md_table(cols, [[_c(x) for x in r] for r in rows])

function _p3t_t1(A, rec)
    th = rec["selection"]["theta"]
    L = ["Calibration record, transcribed (not recomputed). Outcome: $(rec["outcome"]["selection"]); theta*: V = $(join(th["V"], ", ")), q_C = $(th["qC"]), q_R = $(th["qR"]), w = $(th["w"]), v = $(th["v"]).", ""]
    rows = [[k, v["pooled"], v["panel"], v["grid"]] for (k, v) in sort!(collect(rec["sets"]); by = first) if v isa AbstractDict && haskey(v, "pooled")]
    append!(L, mdt(["set of Cal", "pooled", "panel", "grid"], rows), ["", "Values of Amendment A1 and pilot costs (the test-side table):", ""])
    append!(L, mdt(["key", "value"], [[string(k), k === :theta_voters ? join(getfield(A.c.t, k), ", ") : getfield(A.c.t, k)] for k in P3T_PARAM_KEYS]))
    return L
end

function _p3t_t2(A)
    O, E = A.O, A.O.rows
    rds = collect(1:A.entries[end].round)      # E80: the rounds up to the closing round of the log; a later round of the panel is not counted anywhere
    po = [o for o in A.c.panel_orbits if Int(o["round"]) in rds]
    rows = [[r, count(o -> Int(o["round"]) == r, po), count(o -> Int(o["round"]) == r && Float64(o["K"]) == 0, po), count(x -> x.round == r, E),
             count(x -> x.round == r && x.W2, E), count(x -> x.round == r && x.Ts, E), count(x -> x.round == r && x.S2, E)] for r in rds]
    rho = [count(x -> x.rho === s, A.st.rows) for s in (:C, :R, :U)]
    L = ["The test draw. Orbits of the rounds up to the closing round: $(length(po)) (the panel file holds $(length(A.c.panel_orbits))); rounds $(length(rds)).", ""]
    append!(L, mdt(["round", "orbits", "K = 0", "E_test", "W2", "T_s", "S2"], rows))
    append!(L, ["", "Labels of the orbits with K > 0 read (rho = C, R, U): $(rho[1]), $(rho[2]), $(rho[3]); excluded $(count(x -> x.excluded, A.st.rows)), drift on L1 or L5 (as the round log counts) $(count(x -> x.drift, A.st.rows)); dropped as drift on any item read (E64 b) $(count(d -> d.reason == "drift", A.O.drops)).", ""])
    cnt = (E_test = length(E), W2 = count(r -> r.W2, E), T_s = count(r -> r.Ts, E), S2 = count(r -> r.S2, E), D = length(O.vd), D_B = length(A.dB), gate = length(O.gate.ids))
    append!(L, mdt(["set", "size", "target"], [["E_test", cnt.E_test, nothing], ["W2", cnt.W2, A.c.t.target_w2], ["T_s", cnt.T_s, A.c.t.target_t_s], ["S2", cnt.S2, nothing], ["D", cnt.D, A.c.t.target_f],
                                              ["D_B", cnt.D_B, nothing], ["gate", cnt.gate, nothing]]))
    return L
end

# The value of the L5 item at the recorded step N_B (the last recorded step not above it), or `nothing`.
function _p3t_lam_nb(it, nB)
    it === nothing && return nothing
    i = max(searchsortedlast(it.t, nB), 1)
    return it.lam[i]
end

# The ladder test (a) in exact rationals: the stored exponents as binary64 values (M5), the bound as its decimal text (M9).
_p3t_ladder(lam, tol) = any(abs(Rational{BigInt}(lam[i]) - Rational{BigInt}(lam[i + 1])) > p3_rat(tol) for i in 1:2)

# E79: the chk items of D_B by content: method L5, the precision label without its bit count, the horizon of the item itself, which must be the key N_B of the P1 table.
function _p3t_chk_items(A, nB)
    out = Dict{Tuple{Int, String}, Any}()
    for (k, it) in A.idx
        (k[2] == "L5" && "chk" in it.arms) || continue
        it.horizon == nB || _fail("the check item of orbit $(k[1]) (L5, $(k[3])) has the horizon $(it.horizon), not the key N_B = $nB of the P1 parameter table")
        out[(k[1], replace(k[3], r"\(\d+\)" => ""))] = it
    end
    return out
end

function _p3t_t3(A, files)
    if files.p1_params === nothing || !isfile(files.p1_params)
        return ["The P1 parameter table was not given: the reference check on D_B is not computed."]
    end
    p1, p2, d1 = read_params(files.p1_params), A.c.p2, A.c.p2.directions[1]
    chk = _p3t_chk_items(A, p1.N_B)
    rows, forced = Vector{Any}[], Dict{Symbol, Int}()
    for id in A.dB
        its = [get(A.idx, p2_key(id, "L5", "Float64", p2.N_1, 1, d1), nothing), get(chk, (id, "Double64"), nothing), get(chk, (id, "BigFloat"), nothing)]
        lam = [_p3t_lam_nb(it, p1.N_B) for it in its]
        ladder = any(x -> x === nothing || !isfinite(x), lam) ? nothing : _p3t_ladder(lam, p1.ladder_tol)
        rho = haskey(A.byid, id) ? A.byid[id].rho : :dropped      # an orbit of D_B may have left E_test (drift on an item read, E64 b)
        ladder === true && (forced[rho] = get(forced, rho, 0) + 1)
        push!(rows, Any[id, p1.N_B, lam[1], lam[2], lam[3], ladder, rho])
    end
    L = ["The reference check on D_B (N_B = $(p1.N_B), ladder bound $(p1.ladder_tol) of the P1 parameter table): the three L5 values at the recorded step N_B and the ladder test (a). No label changes (D1). " *
         "The column labels agree of the preregistration is omitted: the items of the arm chk run at their own horizon N_B and have no steps at the scale of N_1, so the labels cannot be formed in three precisions (study amendment A2 states it).", ""]
    append!(L, mdt(["orbit", "N_B", "L5 Float64", "L5 Double64", "L5 BigFloat", "ladder would force", "rho"], rows))
    append!(L, ["", "Orbits on which the ladder would have made rho unresolved: C $(get(forced, :C, 0)), R $(get(forced, :R, 0))."])
    return L
end

function _p3t_t4(A)
    E = A.O.rows
    rows = Vector{Any}[]
    for (name, f) in (("E_test", _ -> true), ("W2", r -> r.W2), ("T_s", r -> r.Ts), ("S2", r -> r.S2))
        for (k, rd) in vcat([("all", "all")], sort!(unique((r.kindex, r.round) for r in E)))
            os = [r for r in E if f(r) && (k == "all" || (r.kindex == k && r.round == rd))]
            push!(rows, Any[name, k, rd, count(r -> r.out === :chaotic_supported, os), count(r -> r.out === :regular_supported, os), count(r -> r.out === :unresolved, os)])
        end
    end
    return vcat(["Outputs of the rule R_a with theta* at b_P, per set, K index and round (counts).", ""], mdt(["set", "K index", "round", "chaotic supported", "regular supported", "unresolved"], rows))
end

function _p3t_t5(A)
    h = A.ha
    b(x) = x.value === nothing ? "undefined ($(x.dropped) dropped)" : "$(_c(x.value)) ($(x.dropped) dropped)"
    rows = [["1", "e (error among decided)", h.point.e, "upper bound " * b(h.bounds.e_up), "e_max = $(A.qt.e_max)", 1 in h.failed ? "fails" : "ok"],
            ["2", "a(E_test)", h.point.aE, "upper bound " * b(h.bounds.aE_up), "a_max = $(A.qt.a_max)", 2 in h.failed ? "fails" : "ok"],
            ["3", "a(W2) - a(S2)", h.point.dw, "lower bound " * b(h.bounds.dw_lo), "delta_a = $(A.qt.delta_a)", 3 in h.failed ? "fails" : "ok"],
            ["3", "a(T_s) - a(S2)", h.point.dt, "lower bound " * b(h.bounds.dt_lo), "delta_a = $(A.qt.delta_a)", 3 in h.failed ? "fails" : "ok"],
            ["4", "g (gate orbits chaotic-supported)", h.point.g, "", "g_max = $(A.qt.g_max)", 4 in h.failed ? "fails" : "ok"]]
    L = ["HA: $(h.decision). Strata sizes (T_s, W2 minus T_s, S2, other C, R): $(join(h.strata, ", ")); overlap of T_s and S2: $(count(r -> r.Ts && r.S2, A.O.rows)).", ""]
    append!(L, mdt(["item", "statistic", "point", "bound", "threshold", "status"], rows))
    append!(L, ["", "Unresolved share a(A) per set:"], mdt(["set", "a"], [["E_test", h.point.a[:E]], ["W2", h.point.a[:W2]], ["T_s", h.point.a[:Ts]], ["S2", h.point.a[:S2]]]))
    return L
end

# The exact p-value of a comparison as printed in T6: 12 significant digits; the JSON file keeps the rational.
_p3t_pshow(p) = p isa Rational ? Float64(p) : p

function _p3t_t6(A)
    r, hb = A.hbr, A.hbr.hb
    cmp = p3t_cmp_rows(hb, A.M, A.qt, r.bx)
    rows = [[x.set, x.indicator, x.b, x.cost, x.capped, x.cost_unmatched, x.fr_x, x.fr_rule, x.n11, x.n10, x.n01, x.n00, _p3t_pshow(x.p), x.margin, x.pass] for x in cmp]
    L = ["HB: $(hb.decision). Final sizes: T_s $(r.sizes.Ts), W2 $(r.sizes.W2); underpowered: $(hb.underpowered); cost-unmatched LE: true; defects: $(isempty(hb.defects) ? "none" : join(hb.defects, ", ")).", ""]
    append!(L, mdt(["set", "comparator", "b_X", "cost", "capped", "cost-unmatched", "FR_X", "FR_rule", "n11", "n10", "n01", "n00", "p", "margin", "pass"], rows))
    append!(L, ["", "Coverage of the rule (price 1 - cov): T_s $(_c(r.cov[:Ts])) (price $(_c(1 - r.cov[:Ts]))), W2 $(_c(r.cov[:W2])) (price $(_c(1 - r.cov[:W2]))); condition met: $(hb.coverage_ok).",
                "Without LE (12 comparisons, reported, no effect): $(r.without_le.decision)."])
    return L
end

function _p3t_t7(A)
    s = A.sens
    rows = [[k, getfield(s, k).decidable, getfield(s, k).verdict, join(_p3t_fail_names(getfield(s, k).failing), ", "), getfield(s, k).reason] for k in (:round_down, :pilot, :measured, :without_le)]
    L = ["Sensitivities of HB (reported, no effect on any decision):", ""]
    append!(L, mdt(["sensitivity", "decidable", "verdict", "failing", "reason"], rows))
    pw = A.power
    append!(L, ["", "Power of the full HB decision at the achieved sizes: " * (pw.hb === nothing ? "not recomputed (both targets met)" :
               "$(pw.hb.successes) of $(pw.hb.runs), power $(_fmt(pw.hb.power)), standard error $(_fmt(pw.hb.se)), planned $(_fmt(pw.hb.planned))")])
    return L
end

function _p3t_t8(A)
    h, pw = A.hc, A.power
    L = ["HC: $(h.decision). Set size $(length(h.rows)), failures $(h.failures), non-failures $(h.nonfailures), underpowered: $(h.underpowered); dropped replicates $(h.dropped).", ""]
    append!(L, mdt(["statistic", "value"], [["AUC_I", h.auc_i], ["AUC_P", h.auc_p], ["Delta_AUC", h.delta], ["lower bound", h.lower], ["upper bound", h.upper]]))
    append!(L, ["", "Power at the achieved failures: " * (pw.hc === nothing ? "not recomputed (target reached)" : "$(pw.hc.successes) of $(pw.hc.runs), power $(_fmt(pw.hc.power)), standard error $(_fmt(pw.hc.se)), planned $(_fmt(pw.hc.planned))"), ""])
    append!(L, mdt(["orbit", "rho", "n_C", "n_R", "f", "s_I", "s_P"], [[r.id, r.rho, r.nC, r.nR, r.f, r.sI, r.sP] for r in h.rows]))
    return L
end

function _p3t_t9(A)
    R = A.rab
    names = ("chaotic_supported", "regular_supported", "unresolved")
    rows = [[replace(a, "_" => " "), replace(b, "_" => " "), get(R.pairs, "$a|$b", 0)] for a in names for b in names]
    dec = [r for r in A.O.rows if haskey(R.out, r.id)]
    err(f) = (d = [r for r in dec if f(r) !== :unresolved]; isempty(d) ? nothing : P3T_R(count(r -> (f(r) === :chaotic_supported && r.rho === :R) || (f(r) === :regular_supported && r.rho === :C), d), length(d)))
    L = ["R_a against R_ab on D ($(length(dec)) orbits): the nine pairs of outputs; ladder overrides $(R.overrides).", ""]
    append!(L, mdt(["R_a", "R_ab", "orbits"], rows))
    append!(L, ["", "Error of the decided orbits: R_a $(_c(err(r -> r.out))), R_ab $(_c(err(r -> R.out[r.id])))."])
    return L
end

function _p3t_t10(A)
    q, t = A.qt, A.c.t
    vals, _ = p3t_measured_costs(A.man)
    rows = [[spec, get(vals, Symbol("cost_", lowercase(spec)), nothing), getfield(q, Symbol("cost_", lowercase(spec))), getfield(t, Symbol("pilot_cost_", lowercase(spec)))] for spec in ("L1", "tg", "pa1", "bk1", "fq1")]
    return vcat(["Cost per orbit-step by spec (thread-ns): measured in the test run (manifest), the cost model keys of the 71-key table, and the pilot keys.", ""],
                mdt(["spec", "measured", "cost key", "pilot key"], rows))
end

function _p3t_t11(A)
    cnt = Dict{Tuple{String, String, String}, Int}()
    for it in A.items, f in it.flags
        k = (it.arms[1], it.id, f)
        cnt[k] = get(cnt, k, 0) + 1
    end
    for f in ("drift", "excluded", "diagnostics exceeded"), a in ("ref",)
        haskey(cnt, (a, "L5", f)) || (cnt[(a, "L5", f)] = 0)
    end
    failed = [(c.file, c.cell) for c in A.c.chunks if c.status == "failed"]
    L = ["Stored flags by arm and item kind (counts of items):", ""]
    append!(L, mdt(["arm", "item", "flag", "count"], [[a, i, f, n] for ((a, i, f), n) in sort!(collect(cnt))]))
    append!(L, ["", "Chunks failed twice and their orbits:", ""], mdt(["chunk file", "cell", "orbits"], [[f, c, join(_p3t_cell_orbits(A, c), ", ")] for (f, c) in failed]))
    return L
end

# The orbits of a cell of E74 of a chunk that failed twice, by the panel (E80).
_p3t_cell_orbits(A, cell) = sort!(Int[Int(o["id"]) for o in A.c.panel_orbits if p3t_cell_of(o) == cell])

"The tables file of E73; `sha` is the SHA-256 of the machine file."
function p3t_tables(A, files, sha, rec)
    D = A.decisions
    L = String["# Study p3-abstention: analysis of the test run", "", "Study slug: $(A.c.slug); preregistration header hash: $(A.c.prereg_hash)", "", "results SHA-256: $sha", "",
               "## Decisions", "", "- HA: $(D.HA)", "- HB: $(D.HB)", "- HC: $(D.HC)", ""]
    for (n, title, lines) in ((1, "calibration record", _p3t_t1(A, rec)), (2, "test draw and sets", _p3t_t2(A)), (3, "reference check on D_B", _p3t_t3(A, files)), (4, "outputs of the rule", _p3t_t4(A)),
                              (5, "HA", _p3t_t5(A)), (6, "HB", _p3t_t6(A)), (7, "HB sensitivities and power", _p3t_t7(A)), (8, "HC", _p3t_t8(A)), (9, "R_a against R_ab", _p3t_t9(A)),
                              (10, "cost", _p3t_t10(A)), (11, "flags and health", _p3t_t11(A)))
        append!(L, ["## Table T$n: $title", ""], lines, [""])
    end
    for (id, f) in (("F2", A.f2), ("F3", A.f3))
        append!(L, ["## $id data", ""], md_table(f.columns, f.rows), ["", "digest: $(f.digest)", ""])
    end
    append!(L, ["## Readings chosen", ""], ["- $x" for x in P3T_READINGS], ["", "## Deviations", ""], ["- $x" for x in P3T_DEVIATIONS], ["", "## Exploratory items", "", "none", ""])
    return join(L, "\n")
end
