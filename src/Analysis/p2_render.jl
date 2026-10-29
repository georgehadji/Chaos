# Study p2-indicators (Amendment D89): the tables file of E37. No clock, no absolute path (I13). The data tables of F1 to F6 carry a `digest:` line, the SHA-256 of the
# table written as markdown with 12 significant digits and LF line ends.

const P2_T_TITLES = Dict("T1" => "Orbits by source and K index", "T2" => "Verdict counts by indicator, budget and set", "T3" => "H1 table", "T4" => "H2 table", "T5" => "H3 table",
    "T6" => "H4 table", "T7" => "H5 table", "T8" => "H6 table", "T9" => "AUC of the chaos scores", "T10" => "Statistics and decisions under table R2", "T11" => "Sensitivities",
    "T12" => "Rule-validity gate on the integrable control", "T13" => "I_8 extension and ROBUST flags", "T14" => "Regular range and sticky orbits", "T15" => "Packages",
    "T16" => "Flags, omissions, controls and the K = 0 line", "T17" => "Cost")
const P2_F_TITLES = Dict("F1" => "agreement of the pairs of I_8", "F2" => "detection times and outcomes", "F3" => "fraction chaotic per indicator and budget", "F4" => "AUC rows",
    "F5" => "lattice of the grid", "F6" => "stretching distance by class")

p2_digest(t) = bytes2hex(sha256(codeunits(join(md_table(t.columns, t.rows), "\n")::String)))

function _decision_line(h, d)
    v = _vstr(d.verdict)
    h === :H1 && return "- H1: $v; orbits $(d.n), decided for all indicators $(d.used), U $(_fmt(d.U)) [$(_fmt(d.U_lo)), $(_fmt(d.U_hi))], A_min $(_fmt(d.A_min)), U on S2 $(_fmt(d.U_S2))"
    h === :H2 && return "- H2: $v; orbits $(d.n), Q $(_fmt(d.Q)) [$(_fmt(d.Q_lo)), $(_fmt(d.Q_hi))]"
    h === :H3 && return "- H3: $v; orbits W2 $(d.n_w), S2 $(d.n_s), indicators passing $(d.passing === nothing ? "undefined" : d.passing)"
    h === :H4 && return "- H4: $v; orbits $(d.n), Delta_phi $(_fmt(d.dphi)) [$(_fmt(d.dphi_lo)), $(_fmt(d.dphi_hi))]; design ensemble: $(_vstr(d.design.verdict)), orbits $(d.design.n)"
    h === :H5 && return "- H5: $v; orbits $(d.n), F_med $(_fmt(d.F_med)) [$(_fmt(d.F_med_lo)), $(_fmt(d.F_med_hi))]"
    return "- H6: $v; orbits S_D* $(d.dir.n), S_L* $(d.prec.n), largest delta_dir $(_fmt(_maxdelta(d.dir))), largest delta_prec $(_fmt(_maxdelta(d.prec)))"
end

_maxdelta(part) = (v = filter(!isnothing, collect(values(part.delta))); isempty(v) ? nothing : maximum(v))

function _decision_lines(res)
    L = String["Primary decisions (indicators I_7, sets of E_all, verdicts of table R1), in the order H1 to H6:", ""]
    append!(L, _decision_line(h, res.decisions[h]) for h in (:H1, :H2, :H3, :H4, :H5, :H6))
    append!(L, ["", "Extension (indicators I_8, sets restricted to E8) and the ROBUST flag of each decision:", ""])
    append!(L, "- $h: $(_vstr(res.robust[h])) (primary $(_vstr(res.decisions[h].verdict)), extension $(_vstr(res.extension[h].verdict)))" for h in (:H1, :H2, :H3, :H4, :H5, :H6))
    append!(L, ["", "Table R2 (sensitivity, not decisive; the primary decisions stay those of R1):", ""])
    append!(L, "- $h: $(_vstr(res.r2[h].verdict))" for h in (:H1, :H2, :H3, :H4, :H5, :H6))
    append!(L, ["", "Rule-validity gate on the integrable control (share chaotic at the last budget, limit $(res.params.gate_max_share)):", ""])
    append!(L, "- $(g.X): decided $(g.decided), chaotic $(g.chaotic), share $(_fmt(g.share))$(g.violated ? ", VIOLATED" : "")" for g in res.gate)
    for g in res.gate
        g.violated && push!(L, "", "rule invalid on the integrable control: $(g.X) (share $(_fmt(g.share)) above $(res.params.gate_max_share)); the decisions stay as computed and a second version without $(g.X) is added")
    end
    any(g -> g.violated, res.gate) || push!(L, "", "No indicator violates the gate.")
    return L
end

function _header_lines(res)
    h = res.header
    L = String["- preregistration header SHA-256: $(h.prereg_header_sha256)", "- manifest SHA-256: $(h.manifest_sha256)", "- panel file SHA-256: $(h.panel_sha256)",
               "- parameter table SHA-256: $(h.params_sha256)", "- P1 preregistration header SHA-256: $(h.p1_prereg_header_sha256)", "- P1 manifest SHA-256: $(h.p1_manifest_sha256)",
               "- P1 parameter table SHA-256: $(h.p1_params_sha256)", "- plan hash: $(h.plan_sha256)", "- chunks: $(h.chunks)", "- items: $(h.items)", "- bootstrap seed: $(h.seed)",
               "- module version: $(h.module_version)", "- git commit: $(h.git_sha)", "- parameters (58 keys):"]
    for k in P2_PARAM_KEYS
        v = getfield(res.params, k)
        push!(L, "  - $k = $(v isa AbstractVector ? "[" * join(v, ", ") * "]" : v)")
    end
    return L
end

"The tables file as text."
function p2_render(res)
    L = String["# Analysis tables: p2-indicators", "", "## Header", ""]
    append!(L, _header_lines(res), ["", "## Decisions", ""], _decision_lines(res), [""])
    section(id, title, t; extra = t.extra) = append!(L, ["## $id $title (confirmatory)", ""], md_table(t.columns, t.rows), extra, [""])
    data(id, t) = append!(L, ["## $id data (confirmatory): $(P2_F_TITLES[id])", ""], md_table(t.columns, t.rows), ["", "digest: $(res.digests[id])", ""])
    it = res.items
    for i in 1:17
        id = "T$i"
        section(id, P2_T_TITLES[id], getproperty(it, Symbol(id)))
        for (X, v) in sort!(collect(res.second_items); by = kv -> string(kv[1])), (sid, t) in v
            id == sid && section(id * " without $X", P2_T_TITLES[id], t)
        end
        i == 4 && data("F2", it.F2)
        i == 3 && data("F1", it.F1)
        i == 6 && data("F3", it.F3)
        i == 9 && data("F4", it.F4)
        i == 14 && (data("F5", it.F5); data("F6", it.F6))
    end
    append!(L, ["## Readings", ""], ["- " * r for r in res.readings], ["", "## Deviations", ""], ["- " * d for d in res.deviations], ["", "## Exploratory items", "", "none", ""])
    return join(L, "\n")
end

function p2_write_tables(out::AbstractString, res)
    mkpath(out)
    path = joinpath(out, P2_SLUG * "-tables.md")
    write(path, p2_render(res))
    return path
end
