# Study p3-abstention (Amendment D93): the readable tables file of E58: the header (hashes, versions, the SHA-256 of the record), table T1 (the content of the record in readable form) and the data
# block of figure F1 with its digest. No clock, no absolute path.

_p3_t(x::Rational) = "$(numerator(x))/$(denominator(x))"
_p3_t(x::Nothing) = "undefined"
_p3_t(x) = string(x)
_p3_tr(x::Rational) = "$(_p3_t(x)) ($(round(Float64(x); sigdigits = 6)))"
_p3_tr(x) = _p3_t(x)
_p3_theta_text(t) = "V = " * join(t.V, ", ") * "; qC = $(t.qC), qR = $(t.qR), w = $(t.w), v = $(t.v)"

"""
    p3_f1_table(sel) -> (columns, rows, digest)

The data block of figure F1 (E58): one row per element of F in the order of E52 with `index, theta, C_R, FR(W2; theta), FR(T_s; theta), J, in Theta_tol, theta*`, and the SHA-256 of
the table written as markdown (cells as text, LF line ends, the rule of D89 E37).
"""
function p3_f1_table(sel)
    cols = ["index", "theta", "C_R", "FR(W2; theta)", "FR(T_s; theta)", "J", "in Theta_tol", "theta*"]
    tol = Set(sel.tol)
    rows = [Any[i, _p3_theta_text(sel.thetas[i]), _p3_t(sel.stats[i].cost), _p3_t(sel.stats[i].fr_w2), _p3_t(sel.stats[i].fr_ts), _p3_t(get(sel.J, i, nothing)), i in tol, i == sel.star]
            for i in sel.feasible]
    return (columns = cols, rows = rows, digest = bytes2hex(sha256(codeunits(join(md_table(cols, rows), "\n")::String))))
end

function _p3_header_lines(fp, program)
    L = String["- preregistration header SHA-256 (recomputed under the rule of the SELF cell): $(fp.preregistration)", "- p2 manifest SHA-256: $(fp.p2_manifest)", "- P1 manifest SHA-256: $(fp.p1_manifest)",
               "- panel file SHA-256: $(fp.panel)", "- p3 parameter table SHA-256: $(fp.p3_parameters)", "- p2 parameter table SHA-256: $(fp.p2_parameters)", "- P1 parameter table SHA-256: $(fp.p1_parameters)",
               "- p2 preregistration SHA-256: $(fp.p2_preregistration)", "- P1 preregistration SHA-256: $(fp.p1_preregistration)",
               "- module version: $(program.module_version)", "- git commit: $(fp.git_commit)", "- git working tree clean: $(fp.git_clean)" * (fp.git_clean ? "" : " (not a freeze run)"),
               "- generator: $(program.generator)", "- quantile function: $(program.quantile)"]
    return L
end

function _p3_t1_lines(res, q, rec)
    L = String[]
    sets = rec.sets
    push!(L, "### Relabelling C0 and the drops", "")
    append!(L, md_table(["K index", "old", "new", "count"], [Any[k, o, n, c] for (k, o, n, c) in [(r.k_index, r.old, r.new, r.count) for r in rec.relabelling.changes_panel]]))
    push!(L, "", "Panel: $(rec.relabelling.panel_changed) orbits changed class; grid: $(rec.relabelling.grid_changed). Dropped: $(rec.relabelling.dropped.total) ($(join(["$k $v" for (k, v) in sort(collect(rec.relabelling.dropped.by_reason))], ", "))); relabelling data missing: $(length(rec.relabelling.missing_data)).", "")
    push!(L, "### Sets of Cal (pooled, panel, grid)", "")
    append!(L, md_table(["set", "pooled", "panel", "grid"], [Any[string(n), getfield(sets, n).pooled, getfield(sets, n).panel, getfield(sets, n).grid] for n in (:cal, :S2, :T_s, :W2, :R_W, :regular, :gate)]))
    sel = res.selection
    push!(L, "", "### Feasible set F and the selection", "", "- elements of Theta: $(length(sel.thetas)); cost cap of F1: $(_p3_tr(sel.cost_cap))", "- elements satisfying F1 to F6 alone: $(join(sel.alone, ", "))",
         "- |F| = $(length(sel.feasible)); outcome: $(rec.outcome.selection)", "- J_min: $(_p3_tr(sel.J_min)); |Theta_tol| = $(length(sel.tol))")
    sel.star === nothing && return L
    st, t = sel.stats[sel.star], sel.thetas[sel.star]
    push!(L, "- theta* (index $(sel.star)): $(_p3_theta_text(t)); C_R = $(_p3_tr(st.cost))",
         "- calibration error e = $(_p3_tr(st.e)); unresolved share on Cal $(_p3_tr(st.u_cal)), on S2 $(_p3_tr(st.u_s2)); gate count $(st.g); J = $(_p3_tr(sel.J[sel.star]))", "")
    return L
end

function _p3_c7_lines(res)
    c = res.c7
    L = String["", "### HB analogue on Cal at theta* (C7)", ""]
    for A in (:T_s, :W2)
        push!(L, "Set $A:", "")
        append!(L, md_table(["X", "b", "cost", "capped", "cost-unmatched", "FR_X", "FR_rule", "n10", "n01", "p", "margin", "coverage"],
                            [Any[r.X, r.b, _p3_t(r.cost), r.capped, r.cost_unmatched, _p3_t(r.fr_x), _p3_t(r.fr_rule), r.n10, r.n01, _p3_tr(r.p), r.margin, _p3_t(r.coverage)] for r in c.hb[A]]), [""])
    end
    push!(L, "Risk and coverage of the rule family around theta* and of the single indicators:", "")
    append!(L, md_table(["qC", "qR", "coverage", "e"], [Any[r.qC, r.qR, _p3_t(r.coverage), _p3_t(r.e)] for r in c.risk_coverage.family]), [""])
    append!(L, md_table(["X", "b", "coverage", "e"], [Any[r.X, r.b, _p3_t(r.coverage), _p3_t(r.e)] for r in c.risk_coverage.single]), [""])
    push!(L, "First $(length(c.top)) elements of Theta_tol in the order of C6: " * join(c.top, ", ")::String, "")
    return L
end

function _p3_pp_lines(res)
    L = String["", "### PP for HB and HC", ""]
    if res.hb !== nothing
        h = res.hb
        push!(L, "- optimism: used $(res.optimism.used), skipped $(res.optimism.skipped); r_corr T_s $(round(h.rates.r_corr[:T_s]; sigdigits = 12)), W2 $(round(h.rates.r_corr[:W2]; sigdigits = 12)), R_W $(round(h.rates.r_corr[:R_W]; sigdigits = 12))",
             "- rule rates: T_s $(round(h.rates.p_rule[:T_s]; sigdigits = 12)), R_W $(round(h.rates.p_rule[:R_W]; sigdigits = 12)); targets: $(h.targets.outcome), cap $(h.targets.cap), n_T $(h.targets.n_T), n_W $(h.targets.n_W)", "")
        append!(L, md_table(["j", "rho", "nT", "nRW", "nW", "allowed", "successes", "runs", "power", "se"],
                            [Any[x.j, _p3_t(x.rho), x.nT, x.nRW, x.nW, x.allowed, _p3_t(x.successes), x.runs, _p3_t(x.power), _p3_t(x.se)] for x in h.path]), [""])
        push!(L, "Type-I check (joint and single rejection rates): maximum joint $(h.type1.max_joint)", "")
    end
    hc = res.hc
    hc === nothing || push!(L, hc.undefined ? "- HC: the AUC is undefined (failures $(hc.n_fail), non-failures $(hc.n_nonfail)); no PP for HC" :
         "- HC: AUC_I_cal $(hc.auc_i_cal) (failures $(hc.n_fail), non-failures $(hc.n_nonfail)), AUC_I_used $(hc.auc_i_used), target n_f $(_p3_t(hc.target)), underpowered $(hc.underpowered)", "")
    return L
end

"The tables file as text (E58); `record_sha` is the SHA-256 of the record file (a file cannot state its own hash)."
function p3_render(res, q, rec, fp, program, record_sha, f1)
    L = String["# Calibration tables: p3-abstention", "", "record SHA-256: $record_sha", "", "## Header", ""]
    append!(L, _p3_header_lines(fp, program), ["", "## Table T1", ""], _p3_t1_lines(res, q, rec))
    res.c7 === nothing || append!(L, _p3_c7_lines(res))
    append!(L, _p3_pp_lines(res), ["", "## F1 data", ""], md_table(f1.columns, f1.rows), ["", "digest: $(f1.digest)", ""])
    return join(L, "\n")
end
