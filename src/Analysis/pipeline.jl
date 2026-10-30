# The analysis pipeline (Algorithm of the spec): checks, items, per-orbit quantities, controls, sets, decisions H1 to H5, outputs.

default_seed() = Int(TOML.parsefile(joinpath(@__DIR__, "..", "..", "_config", "numerics-defaults.toml"))["seed"])

"The function that draws the figure files from the result; set by the plotting extension when it loads. Without it the figures are skipped."
const FIGURE_DRAWER = Ref{Any}(nothing)
figures_available() = FIGURE_DRAWER[] !== nothing

# Union of the flags of the items of each orbit.
function _flag_map(items)
    m = Dict{Int, Set{String}}()
    for it in items
        isempty(it.flags) || union!(get!(m, it.orbit, Set{String}()), it.flags)
    end
    return m
end

# Directions of the L5 items at the primary settings and horizon N_1 of each orbit.
function _dirs_map(items, p)
    m = Dict{Int, Vector{String}}()
    for it in items
        it.method == "L5" && it.precision == "Float64" && it.mode == "lifted" && it.tau == 1 && it.eps === nothing && it.horizon == p.N_1 &&
            push!(get!(m, it.orbit, String[]), it.direction)
    end
    return m
end

_bf_label(items) = (i = findfirst(it -> startswith(it.precision, "BigFloat"), items); i === nothing ? nothing : items[i].precision)

# Controls of the main arm with their stored scope and deviation (E9).
function _control_rows(idx, orbits, p)
    rows = NamedTuple[]
    for o in orbits
        o["control"] || continue
        id = Int(o["id"])
        its = Dict(m => _primary(idx, id, m, p.N_1, p) for m in V_METHODS)
        scope = get(its["L5"].attrs, "control_scope", nothing)
        if scope === nothing
            Float64(o["exponent"]) > 0 && throw(ArgumentError("control $id has a positive exponent but no control_scope attribute in its main items"))
            continue
        end
        stored = Dict{String, Float64}(m => Float64(it.attrs["control_deviation"]) for (m, it) in its if haskey(it.attrs, "control_deviation"))
        push!(rows, (id = id, exponent = Float64(o["exponent"]), scope = String(scope), lam = Dict(m => _val(it, p.N_1) for (m, it) in its),
                     stored = isempty(stored) ? nothing : stored))
    end
    return rows
end

# The two H4 starts (extra ids above the panel) and their cross-check values (E13).
function _h4_data(idx, items, nP, p)
    rows = NamedTuple[]
    for id in sort!(unique(it.orbit for it in items if it.orbit > nP))
        K = parse(Float64, first(it for it in items if it.orbit == id).attrs["params"])
        f1, f5 = _need(idx, id, "L1", "Float64", p.h4_T), _need(idx, id, "L5", "Float64", p.h4_T)
        d5 = _need(idx, id, "L5", "Double64", p.h4_T_cross)
        push!(rows, (id = id, K = K, L1 = _val(f1, p.h4_T), L5 = _val(f5, p.h4_T), float64 = _val(f5, p.h4_T_cross), double64 = _val(d5, p.h4_T_cross)))
    end
    pick(K) = (i = findfirst(r -> abs(r.K - K) <= 1e-9, rows); i === nothing ? throw(ArgumentError("the store has no H4 start at K = $K")) : rows[i])
    return (halfpi = pick(pi / 2), twopi = pick(2 * pi))
end

# Ensemble mean of lambda_L5(N_1) over the kept chaotic-sea orbits of each panel K, next to the comparators (E13); no verdict.
function _comparators(orbits, recs, sets, p)
    keep = Set(sets.E_ch)
    rows = NamedTuple[]
    for k in sort!(unique(Int(o["k_index"]) for o in orbits if !o["control"]))
        os = [(o, r) for (o, r) in zip(orbits, recs) if !o["control"] && Int(o["k_index"]) == k && o["class"] == "C" && Int(o["id"]) in keep]
        isempty(os) && continue
        K = Float64(first(os)[1]["K"])
        push!(rows, (k_index = k, K = K, n = length(os), mean = sum(r.lam1["L5"] for (_, r) in os) / length(os), shevchenko = shevchenko_fit(K), chirikov = chirikov_line(K)))
    end
    return rows
end

function _cells(orbits, recs, sets)
    dropped = Set(d.id for d in sets.dropped)
    groups = Dict{Tuple{Int, String}, Vector{Int}}()
    for (i, o) in enumerate(orbits)
        push!(get!(groups, o["control"] ? (0, "controls") : (Int(o["k_index"]), String(o["class"])), Int[]), i)
    end
    return [(k_index = k, class = c, orbits = length(ix), dropped = count(i -> Int(orbits[i]["id"]) in dropped, ix),
             labels = Dict{Tuple{String, String}, Int}(l => count(i -> (recs[i].label.instability, recs[i].label.structure) == l, ix) for l in
                         unique((recs[i].label.instability, recs[i].label.structure) for i in ix))) for ((k, c), ix) in sort!(collect(groups))]
end

function _deviations(p, ctrl, man, recs)
    bound = sort!([r.id for r in ctrl.rows if r.scope == "bound" && r.exponent > 0])
    info = sort!([r.id for r in ctrl.rows if r.scope == "info"])
    ov = get(man, "overrides", Dict())
    nbr = count(r -> "breakdown_tau" in r.flags, recs)
    return String["L2 runs in wrapped coordinates (A5 item 1); H1c compares the wrapped L2 value with the lifted L5 value of the same orbit",
                  "The tau sweep records an orthonormalisation breakdown as a result, not as a failure (A5 item 2); orbits with the flag: $nbr",
                  "Control bound scope (A5 item 3): under the bound $(join(bound, ", ")); reported as info $(join(info, ", "))",
                  "The control bound of L9 is $(p.control_bound_L9), of the other methods $(p.control_bound) (A6)",
                  isempty(ov) ? "Overrides recorded in the manifest: none" : "Overrides recorded in the manifest: $(join(sort!(collect(keys(ov))), ", "))"]
end

const READINGS = ["R(n) and n_conv use the horizon-N_1 series of every orbit; lambda_ref and the labels use the longest horizon of the orbit (P3)",
    "A failing label test forces both axes to unresolved (P4)",
    "The second-half rotation number is (x_N - x_h) / (2 pi (N - h)); h is the recorded step nearest N/2, N/2 itself where it is recorded (P5)",
    "Labels across the three precisions are compared at N_B (P5)",
    "H5a is read from the per-precision labels before the ladder test acts; H5 is supported only if H5a, H5b and H5c all pass (P6)",
    "Percentile bootstrap, orbits resampled with replacement, W and S resampled independently, streams keyed by seed and label (P7)",
    "F1 shows the first drawn orbit of each cell, F3 empirical distribution functions, F4 shares per test, F5 quantiles, F7 paired differences (P8)",
    "S is the class C orbits at K = 2 pi and K = 10 without a label condition; W takes class C at 0 < K <= K_c; E_ch takes class C at K > 0",
    "Periodic points are the controls of kind fixed or newton; kind y5 is not periodic",
    "The ladder compares the three precisions with one common half point: the recorded step nearest N_B/2 that all three items have",
    "Diagnostics exceeded is read on the L5 item at N_1 and on the L5 item at the label horizon",
    "F6 is the unweighted mean of thread-nanoseconds per orbit-step over all items of a method and precision (all arms, all sweep values)",
    "F2 has a pooled series over W next to the series by class and K; orbits with lambda_ref <= 0 are left out of F2",
    "Flags drift and excluded are the union over the items of an orbit; diagnostics exceeded is read on the primary L5 item of arm main"]

function _write_defect(out, study::String, ctrl, hashes)
    mkpath(out)
    path = joinpath(out, study * "-defect.md")::String
    open(path, "w") do io
        println(io, "# Control defect: $study\n\nA gated exact-exponent control exceeds its bound. This is a code defect, not a scientific result; the analysis stops (spec E9).\n")
        for d in ctrl.defects
            println(io, "- control $(d.id), method $(d.method): deviation $(d.deviation) above the bound $(d.bound)")
        end
    end
    return path
end

# The checks, the items and the control rows that the labelling and the full analysis of study p1-le-benchmark share (E17).
function _p1_prepare(prereg, manifest, panel, params, store)
    c = check_inputs(prereg, manifest, panel, params, store)
    items = read_store(manifest, store)
    check_items(items, c.man, c.p)
    orbits = TOML.parsefile(panel)["orbit"]
    idx = index_items(items)
    ctrl = _controls(idx, orbits, c.p)
    return (c = c, items = items, orbits = orbits, idx = idx, ctrl = ctrl)
end

"""
    analyze_p1(; prereg, manifest, panel, params, store, out, figures = false, threads, seed) -> NamedTuple

Run the analysis of the spec for study `p1-le-benchmark`. Every input is checked before anything is written (an `ArgumentError` names the cause). A gated control above
its bound stops the run with status `:control_defect` and a defect report; otherwise the tables file `<study>-tables.md` is written and, when
`figures` is true and the plotting extension is loaded, the figures. Status `:ok`, `:figures_skipped` or `:control_defect`.
"""
function analyze_p1(; prereg, manifest, panel, params, store, out, figures::Bool = false, threads::Integer = Threads.nthreads(), seed::Union{Nothing, Integer} = nothing)
    pre = _p1_prepare(prereg, manifest, panel, params, store)
    c, items, orbits, idx, ctrl = pre.c, pre.items, pre.orbits, pre.idx, pre.ctrl
    p, man = c.p, c.man
    sd = seed === nothing ? default_seed() : Int(seed)
    study = String(man["study"])
    ctrl.ok || return merge(_empty_result(:control_defect, c.hashes), (defects = ctrl.defects, tables = _write_defect(out, study, ctrl, c.hashes)))
    flags, dirs, bf = _flag_map(items), _dirs_map(items, p), _bf_label(items)
    recs = [orbit_record(idx, o, p, flags, dirs, bf) for o in orbits]
    res = _decide(orbits, recs, idx, items, ctrl, c, sd, threads)
    recheck_inputs(c.raw, prereg, manifest, panel, params)
    path = write_tables(out, study, res)
    return merge(res, (status = _figures(res, out, figures), tables = path))
end

# Draw the figures when asked: `:ok`, `:figures_skipped` (no plotting package loaded) or `:figures_failed` (the drawer raised; the tables are written).
_figures(res, out, wanted::Bool) = _draw_figures(res, out, wanted).status

# The same with the message of the exception, which the command of study p2-indicators prints (E45); `error` is empty unless the drawer raised.
function _draw_figures(res, out, wanted::Bool)
    wanted || return (status = :ok, error = "")
    figures_available() || return (status = :figures_skipped, error = "")
    try
        FIGURE_DRAWER[](res, out)
        return (status = :ok, error = "")
    catch e
        @warn "drawing the figures failed" exception = (e, catch_backtrace())
        return (status = :figures_failed, error = sprint(showerror, e))
    end
end

_empty_result(status, hashes) = (status = status, tables = nothing, decisions = nothing, sets = nothing, labels = nothing, controls = nothing, cells = nothing,
                                 digests = nothing, deviations = String[], readings = READINGS, header = hashes, defects = NamedTuple[])

# Control rows, the check and the refusal on stored deviations that do not match the series (E9, I4).
function _controls(idx, orbits, p)
    rows = _control_rows(idx, orbits, p)
    cc = check_controls(rows, p)
    isempty(cc.mismatches) || throw(ArgumentError("control_deviation of control $(join(cc.mismatches, ", ")) differs from the deviation recomputed from the stored series"))
    return merge(cc, (rows = rows,))
end

# Sets, the five decisions in the fixed order (I5), the labels, cells, controls, deviations and the figure data.
function _decide(orbits, recs, idx, items, ctrl, c, sd, threads)
    p = c.p
    sets_in = [(id = r.id, K = Float64(get(o, "K", NaN)), class = String(o["class"]), control = Bool(o["control"]), exponent = o["control"] ? Float64(o["exponent"]) : 0.0,
            flags = r.flags, instability = r.label.instability) for (o, r) in zip(orbits, recs)]
    sets = orbit_sets(sets_in, p)
    byid = Dict(r.id => r for r in recs)
    h1 = decide_h1([(id = i, pass_a = byid[i].h1.pass_a, pass_b = byid[i].h1.pass_b, pass_c = byid[i].h1.pass_c) for i in sets.E_ch], p)
    h2 = decide_h2([byid[i].R_nH2 for i in sets.W], p; seed = sd)
    h3 = decide_h3([byid[i].conv.n for i in sets.W], [byid[i].conv.n for i in sets.S], p; seed = sd)
    wsn = vcat(sets.W, sets.S)
    h3 = merge(h3, (undefined_share = isempty(wsn) ? 0.0 : count(i -> !byid[i].conv.defined, wsn) / length(wsn),))
    h4d = _h4_data(idx, items, maximum(Int(o["id"]) for o in orbits), p)
    h4 = merge(decide_h4((halfpi = (L1 = h4d.halfpi.L1, L5 = h4d.halfpi.L5), twopi = (L1 = h4d.twopi.L1, L5 = h4d.twopi.L5)), p),
               (crosscheck = [(K = x.K, float64 = x.float64, double64 = x.double64) for x in (h4d.halfpi, h4d.twopi)], comparators = _comparators(orbits, recs, sets, p)))
    h5 = decide_h5([(id = r.id, in_SL = r.in_SL, ladder_same = r.label.ladder_same, instability = r.label.instability, structure = r.label.structure,
                     excluded = :excluded in r.label.forced) for r in recs], p)
    labels = Dict{Int, Any}(r.id => (instability = r.label.instability, structure = r.label.structure, forced = r.label.forced, raw = r.label.raw,
                                      ladder_same = r.label.ladder_same, applicable = r.applicable) for r in recs)
    fig = figure_data(orbits, recs, items, idx, sets, p, sd, threads)
    return (status = :ok, tables = nothing, decisions = (H1 = h1, H2 = h2, H3 = h3, H4 = h4, H5 = h5), sets = sets, labels = labels,
            controls = [(id = r.id, scope = r.scope, deviations = control_deviations(r.lam, r.exponent)) for r in ctrl.rows],
            cells = _cells(orbits, recs, sets), digests = figure_digests(fig), deviations = _deviations(p, ctrl, c.man, recs), readings = READINGS,
            header = merge(c.hashes, (seed = sd, module_version = string(Experiments.package_version()),
                                               git_sha = Experiments.git_state(pkgdir(parentmodule(@__MODULE__))::String).sha)), defects = NamedTuple[], study = String(c.man["study"]), figures = fig, records = recs, orbits = orbits, params = p, h4 = h4d,
            panel_cells = get(get(TOML.parsefile(c.panel_path), "manifest", Dict()), "cells", []))
end
