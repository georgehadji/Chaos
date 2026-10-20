# The tables file (E19): T1 to T7, the figure data with digests, the readings, the deviations. No clock, no absolute path (I3).

_fmt(x::AbstractFloat) = isfinite(x) ? string(round(x; sigdigits = 12)) : string(x)
_fmt(x) = string(x)

md_table(cols, rows) = vcat("| " * join(cols, " | ") * " |", "|" * repeat("---|", length(cols)), ["| " * join(_fmt.(r), " | ") * " |" for r in rows])

figure_lines(t) = vcat(md_table(t.columns, t.rows), hasproperty(t, :extra) ? t.extra : String[])
figure_digests(fig) = Dict{String, String}(k => bytes2hex(sha256(codeunits(join(figure_lines(t), "\n")::String))) for (k, t) in fig)

_verdict(v) = replace(string(v), "_" => " ")
_cell(o) = o["control"] ? "controls" : "K index $(o["k_index"]), class $(o["class"])"

function _t1(res)
    lines = ["Orbits per cell of the panel (controls in their own cell), orbits dropped from the sets (flags drift or excluded).", ""]
    append!(lines, md_table(["K index", "class", "orbits", "dropped"], [[c.k_index, c.class, c.orbits, c.dropped] for c in res.cells]))
    append!(lines, ["", "Cells of the panel manifest:", ""])
    append!(lines, md_table(["K index", "class", "drawn", "qualifying", "status", "reason"],
                            [[c["k_index"], c["class"], c["drawn"], c["qualifying"], c["status"], c["reason"]] for c in res.panel_cells]))
    append!(lines, ["", "Orbits dropped as drift or excluded: " * (isempty(res.sets.dropped) ? "none" : join(["$(d.id) ($(d.reason))" for d in res.sets.dropped], ", ")::String)])
    return lines
end

function _t2(res)
    h1 = res.decisions.H1
    byid = Dict(r.id => r for r in res.records)
    lines = ["Verdict: $(_verdict(h1.verdict)); eligible orbits $(h1.n); pass fractions H1a $(_fmt(h1.fractions.a)), H1b $(_fmt(h1.fractions.b)), H1c $(_fmt(h1.fractions.c)) " *
             "(passing: $(h1.passes.a), $(h1.passes.b), $(h1.passes.c)).", ""]
    ks = Dict(Int(o["id"]) => (o["control"] ? 0 : Int(o["k_index"])) for o in res.orbits)
    rows = [[k == 0 ? "controls" : "K index $k", length(ix), count(i -> byid[i].h1.pass_a, ix), count(i -> byid[i].h1.pass_b, ix), count(i -> byid[i].h1.pass_c, ix)]
            for (k, ix) in sort!(collect(Dict(k => [i for i in res.sets.E_ch if ks[i] == k] for k in unique(ks[i] for i in res.sets.E_ch))))]
    append!(lines, md_table(["cell", "orbits", "H1a pass", "H1b pass", "H1c pass"], rows))
    fails = [[i, byid[i].h1.delta_qb, byid[i].h1.diff_l9, byid[i].h1.diff_l2] for i in res.sets.E_ch if !(byid[i].h1.pass_a && byid[i].h1.pass_b && byid[i].h1.pass_c)]
    append!(lines, ["", "Failing orbits:", ""])
    append!(lines, isempty(fails) ? ["none"] : md_table(["orbit", "Delta_QB", "|L9 - L5|", "|L2 - L5|"], fails))
    append!(lines, ["", "Exact-exponent controls (deviation of lambda_1(N_1) from the analytic exponent, largest over V):", ""])
    append!(lines, md_table(["control", "scope", "largest deviation", "L9 deviation"], [[c.id, c.scope, maximum(values(c.deviations)), c.deviations["L9"]] for c in res.controls]))
    return lines
end

# Median and bootstrap interval of a value per panel cell for the orbits of one set.
function _per_cell(res, ids, vals, tag)
    p = res.params
    idset = Set(ids)
    rows = Vector{Any}[]
    for ((cls, k), ix) in _cell_groups(res.orbits)
        sel = [vals[Int(res.orbits[i]["id"])] for i in ix if Int(res.orbits[i]["id"]) in idset]
        isempty(sel) && continue
        b = bootstrap_interval(sel, "$tag/$cls/$k"; seed = res.header.seed, resamples = p.bootstrap_resamples, level = p.bootstrap_level)
        push!(rows, Any[cls, k, length(sel), b.point, b.lo, b.hi, length(sel) < p.min_set ? "not decidable" : ""])
    end
    return rows
end

function _t3(res)
    h = res.decisions.H2
    byid = Dict(r.id => r for r in res.records)
    lines = ["Verdict: $(_verdict(h.verdict)); W orbits $(h.n); median R(n_H2) $(_fmt(h.median)), 95 percent bootstrap interval [$(_fmt(h.lo)), $(_fmt(h.hi))].", "", "By cell (W orbits):", ""]
    return append!(lines, md_table(["class", "K index", "orbits", "median", "lo", "hi", "note"], _per_cell(res, res.sets.W, Dict(i => byid[i].R_nH2 for i in res.sets.W), "T3")))
end

function _t4(res)
    h = res.decisions.H3
    byid = Dict(r.id => r for r in res.records)
    n(ids) = isempty(ids) ? NaN : median_of([byid[i].conv.n for i in ids])
    lines = ["Verdict: $(_verdict(h.verdict)); median n_conv over W $(_fmt(n(res.sets.W))), over S $(_fmt(n(res.sets.S))); ratio $(_fmt(h.ratio)), 95 percent bootstrap interval " *
             "[$(_fmt(h.lo)), $(_fmt(h.hi))]; orbits W $(h.n_w), S $(h.n_s); share of undefined n_conv $(_fmt(h.undefined_share)); S orbits without positive growth $(res.sets.S_nonpositive).",
             "", "By cell (W and S orbits):", ""]
    ids = vcat(res.sets.W, res.sets.S)
    return append!(lines, md_table(["class", "K index", "orbits", "median", "lo", "hi", "note"], _per_cell(res, ids, Dict(i => Float64(byid[i].conv.n) for i in ids), "T4")))
end

function _t5(res)
    h = res.decisions.H4
    rows = [["K = pi/2", h.halfpi.L1, h.halfpi.L5, res.params.h4_ref_halfpi, res.params.h4_tol_halfpi, h.halfpi.pass],
            ["K = 2 pi", h.twopi.L1, h.twopi.L5, res.params.h4_ref_2pi, res.params.h4_tol_2pi, h.twopi.pass]]
    lines = ["Verdict: $(_verdict(h.verdict)).", ""]
    append!(lines, md_table(["start", "L1", "L5", "reference", "tolerance", "pass"], rows))
    append!(lines, ["", "Cross-check (no verdict): L5 at the cross-check step, Float64 against Double64:", ""])
    append!(lines, md_table(["K", "Float64", "Double64"], [[x.K, x.float64, x.double64] for x in h.crosscheck]))
    append!(lines, ["", "Comparators (no verdict): ensemble mean of lambda_L5(N_1) over kept chaotic-sea orbits, the Shevchenko2014 fit and the line ln(K/2):", ""])
    return append!(lines, md_table(["K index", "K", "orbits", "mean", "Shevchenko2014", "ln(K/2)"], [[c.k_index, c.K, c.n, c.mean, c.shevchenko === nothing ? "no fit" : c.shevchenko, c.chirikov] for c in h.comparators]))
end

function _t6(res)
    rows = [[c.k_index, c.class, l[1], l[2], k] for c in res.cells for (l, k) in sort!(collect(c.labels))]
    lines = ["Two-axis label counts per cell:", ""]
    append!(lines, md_table(["K index", "class", "instability", "structure", "orbits"], rows))
    append!(lines, ["", "Labels forced to unresolved, by test (could be applied to / forced):", ""])
    tests = [:ladder, :direction, :estimator, :diagnostics, :excluded]
    return append!(lines, md_table(["test", "applicable", "forced"], [[string(t), count(r -> t in r.applicable || t === :excluded, res.records), count(r -> t in r.label.forced, res.records)] for t in tests]))
end

function _t7(res)
    h = res.decisions.H5
    lines = ["Verdict: $(_verdict(h.verdict)).", ""]
    rows = [[nm, s.n, s.fraction, s.pass === nothing ? "not decidable" : (s.pass ? "pass" : "fail")] for (nm, s) in (("H5a", h.h5a), ("H5b", h.h5b), ("H5c", h.h5c))]
    append!(lines, md_table(["sub-rule", "orbits", "fraction", "result"], rows))
    h.h5b.pass === false || h.h5c.pass === false ? push!(lines, "", "The label scheme is not informative; decisions D36 or D37 go back to the owner.") : nothing
    return lines
end

const FIG_TITLES = Dict("F1" => "lambda_m(n) curves", "F2" => "R(n) median and band", "F3" => "n_conv distributions", "F4" => "label stability",
                        "F5" => "symplectic residuals", "F6" => "cost", "F7" => "sensitivity")
const T_TITLES = Dict("T1" => "Panel table", "T2" => "H1 table", "T3" => "H2 table", "T4" => "H3 table", "T5" => "H4 table", "T6" => "Label table", "T7" => "H5 table")

"The tables file as text."
function render_tables(study::AbstractString, res)
    L = String["# Analysis tables: $study", "", "## Header", ""]
    h = res.header
    append!(L, ["- preregistration header SHA-256: $(h.prereg_header_sha256)", "- panel file SHA-256: $(h.panel_sha256)", "- manifest SHA-256: $(h.manifest_sha256)",
                "- parameter file SHA-256: $(h.params_sha256)", "- bootstrap seed: $(h.seed)", "- module version: $(h.module_version)", "- git commit: $(h.git_sha)", "- parameters: " * join(["$k = $(getfield(res.params, k))" for k in sort!(collect(keys(res.params)))], "; ")::String, ""])
    section(id, body) = append!(L, ["## $id $(T_TITLES[id]) (confirmatory)", ""], body, [""])
    fsection(id) = append!(L, ["## $id data (confirmatory): $(FIG_TITLES[id])", ""], figure_lines(res.figures[id]), ["", "digest: $(res.digests[id])", ""])
    section("T1", _t1(res)); section("T2", _t2(res)); fsection("F1"); fsection("F2"); section("T3", _t3(res)); fsection("F3"); section("T4", _t4(res))
    section("T5", _t5(res)); section("T6", _t6(res)); fsection("F4"); section("T7", _t7(res)); fsection("F5"); fsection("F6"); fsection("F7")
    append!(L, ["## Readings", ""], ["- " * r for r in res.readings], ["", "## Deviations", ""], ["- " * d for d in res.deviations], ["", "## Exploratory items", "", "none", ""])
    return join(L, "\n")
end

function write_tables(out::AbstractString, study::AbstractString, res)
    mkpath(out)
    path = joinpath(out, study * "-tables.md")
    write(path, render_tables(study, res))
    return path
end
