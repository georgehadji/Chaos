# Study p3-abstention (Amendment D93): the pipeline of the Algorithm of the amendment, in the order of its steps. Nothing is written before every check has passed.

# E48 g for the P1 items read: the stored start text (and K, if the item has it) equals the panel file.
function _p3_check_p1_identity(idx1, panel_orbits, p2)
    N, d1 = p2.N_1, p2.directions[1]
    for o in panel_orbits
        o["control"] === true && continue
        id = Int(o["id"])
        for m in ("L5", "L1")
            it = get(idx1, item_key(id, m, "Float64", N, 1, nothing, d1, "lifted"), nothing)
            it === nothing && continue
            a = it.attrs
            lbl = "P1 chunk $(_chunk_name(it.chunk)), orbit $id, method $m"
            haskey(a, "ic") && a["ic"] != "$(o["x0_dec"]),$(o["p0_dec"])" && _fail("$lbl: stored ic \"$(a["ic"])\" differs from the start of the panel file \"$(o["x0_dec"]),$(o["p0_dec"])\"")
            haskey(a, "params") && haskey(o, "K") && tryparse(Float64, String(a["params"])) != Float64(o["K"]) &&
                _fail("$lbl: stored params \"$(a["params"])\" differs from the K $(o["K"]) of the panel file")
        end
    end
end

# The flags of an orbit: the union over the items read (p2 items, and the P1 items of C0 for a panel orbit).
function _p3_flags_of(items2, idx1, panel_orbits, p2)
    fl = Dict{Int, Set{String}}()
    for it in items2
        union!(get!(fl, it.orbit, Set{String}()), it.flags)
    end
    for o in panel_orbits, m in ("L5", "L1")
        it = get(idx1, item_key(Int(o["id"]), m, "Float64", p2.N_1, 1, nothing, p2.directions[1], "lifted"), nothing)
        it === nothing || isempty(it.flags) || union!(get!(fl, it.orbit, Set{String}()), it.flags)
    end
    return fl
end

# The stages up to the verdicts: items, relabelling, Cal and the sets, the verdicts of Cal and of the gate set (E48 to E51).
function _p3_prepare(i, q, p2, c2, p1)
    panel_orbits = TOML.parsefile(i.panel)["orbit"]
    grid = p2_grid_starts(i.p2_prereg; panel = i.panel)
    _check_grid(p2, grid)
    items2 = _p3_read_p2(c2.man, i.p2_store, p2, panel_orbits)
    p2_check_items(items2, c2, panel_orbits, grid)
    idx2 = p2_index(items2)
    _p3_check_p1_identity(p1.idx, panel_orbits, p2)
    nonctl = vcat(Int[Int(o["id"]) for o in panel_orbits if o["control"] !== true], Int[g.id for g in grid])
    p2_check_coverage(idx2, c2.man, [(id = "L1", precision = "Float64", horizon = p2.N_1, direction = p2.directions[1], orbits = nonctl)])
    ctx = (p = p2, idx = idx2)
    rel = p3_relabel_orbits(panel_orbits, grid, p1.idx, idx2, p2, p1.old)
    orbits = p3_orbit_rows(panel_orbits, grid, ctx, rel, _p3_flags_of(items2, p1.idx, panel_orbits, p2))
    S = p3_sets(orbits, q, p2)
    need = sort!(union(S.cal, S.gate))
    p2_check_coverage(idx2, c2.man, [(id = id, precision = pr, horizon = h, direction = d, orbits = os) for (id, pr, h, d, os) in _p3_verdict_requirements(p2, need)])
    verdicts = p3_verdicts(ctx, need)
    return (orbits = orbits, S = S, rel = rel, ctx = ctx, verdicts = verdicts, rows = p3_rows(orbits, S, verdicts, p2))
end

# C7 and PP (E54 C7, E56, E57) for theta*; `nothing` fields for the outcomes without a selection.
function _p3_pp(prep, q, sel, threads)
    none = (star = nothing, matched = nothing, c7 = nothing, optimism = nothing, hb = nothing, hc = nothing, rab = NamedTuple[])
    sel.star === nothing && return none
    th, rows = sel.thetas[sel.star], prep.rows
    m = p3_matched(q, th.V)
    any(X -> m.X[X].capped, Symbol.(q.voter_bits)) && _fail("defect: the feasible element $(sel.star) has a capped comparator (F1 forbids it)")
    c7 = p3_c7(rows, q, sel)
    opt = p3_optimism(rows, q; threads = threads)
    rates = p3_hb_rates(rows, q, th, opt.optimism)
    model = p3_hb_model(rates)
    path = p3_hb_path(model, q; runs = q.pp_runs, threads = threads)
    tg = p3_hb_targets(path, q)
    pt = tg.j === nothing ? nothing : path[tg.j + 1]
    type1 = pt === nothing ? nothing : p3_hb_type1(model, pt.nT, pt.nRW, q; point = pt.j, runs = q.pp_runs, threads = threads)
    pop = p3_hc_population(prep.rows, q, th)
    hc = pop.auc === nothing ? (undefined = true, auc_i_cal = nothing, n_fail = pop.n_fail, n_nonfail = pop.n_nonfail) :
         merge(p3_pp_hc(q, pop.auc; runs = q.pp_runs, threads = threads), (auc_i_cal = pop.auc, n_fail = pop.n_fail, n_nonfail = pop.n_nonfail, undefined = false))
    return (star = th, matched = m, c7 = c7, optimism = opt, hb = (rates = rates, model = model, path = path, targets = tg, type1 = type1), hc = hc, rab = p3_rab(prep.ctx, rows, th, _p3_kP(q)))
end

"""
    calibrate(; prereg, params, out, threads, p2_store, p2_manifest, p2_prereg, p2_params, panel, p1_store, p1_manifest, p1_prereg, p1_params) -> NamedTuple

The calibration of study `p3-abstention` (Amendment D93): checks every input (an `ArgumentError` names the cause; nothing is written, the output folder is not created), relabels the
orbits of the p2 panel and grid (C0), builds Cal and the sets, computes the verdicts, runs C3 to C8, the optimism bootstrap and PP for HB and HC, and writes the record
`<out>/p3-abstention-calibration-record.json` and the tables `<out>/p3-abstention-calibration-tables.md` (T1 and the data of figure F1). Files not given default to the study layout.
"""
function calibrate(; prereg, params, out = nothing, threads::Integer = Threads.nthreads(), kw...)
    i = _p3_inputs(prereg, params, out, kw)
    _p3_check_out(i.out)
    slug = study_slug(prereg)
    slug == P3_SLUG || _fail("calibrate is the command of study $P3_SLUG; this preregistration is of study \"$slug\"")
    ph = _p2_check_prereg(prereg)
    qp = _p3_read(params, prereg, i.p2_params, i.p2_prereg, i.p2_manifest, i.p1_params, nothing)
    q, p2 = qp.q, qp.p2
    c2 = p2_check_inputs(i.p2_prereg, i.p2_manifest, i.panel, i.p2_params, i.p2_store, i.p1_params)
    c1 = _p3_check_p1(i, c2.man)
    p1 = _p3_p1_old(c1, i)
    files = [i.prereg, i.params, i.p2_prereg, i.p2_manifest, i.p2_params, i.panel, i.p1_prereg, i.p1_manifest, i.p1_params]
    raw = _p3_hashes(files)
    prep = _p3_prepare(i, q, p2, c2, p1)
    sel = p3_select(prep.rows, q)
    pp = _p3_pp(prep, q, sel, threads)
    return _p3_finish(i, q, c2, c1, ph, raw, files, prep, sel, pp)
end

# The fingerprints of C1 (E55): the hashes of the inputs, the commit of the library and whether its working tree was clean.
function _p3_fingerprints(i, c2, c1, ph, raw)
    g = git_state(pkgdir(parentmodule(@__MODULE__))::String)
    return (preregistration = ph, p2_manifest = c2.hashes.manifest_sha256, p1_manifest = c1.hashes.manifest_sha256, panel = c2.hashes.panel_sha256, p3_parameters = raw[i.params],
            p2_parameters = c2.hashes.params_sha256, p1_parameters = c1.hashes.params_sha256, p2_preregistration = c2.hashes.prereg_header_sha256,
            p1_preregistration = c1.hashes.prereg_header_sha256, git_commit = g.sha, git_clean = !g.dirty)
end

const P3_PROGRAM = (name = "Chaoslab.Analysis.calibrate", generator = "xoshiro256++ seeded by SplitMix64 from the first 8 bytes, big-endian, of SHA-256 over chaoslab-analysis-1|<test_seed>|p3|<purpose>|<unit>",
                    quantile = "AS241 (PPND16)")

# The result, the record and the tables; the input hashes are checked again before the two files are written (E48 j).
function _p3_finish(i, q, c2, c1, ph, raw, files, prep, sel, pp)
    src = Dict(o.id => o for o in prep.orbits)
    dropped = [(id = d.id, source = src[d.id].source, kindex = src[d.id].kindex, reason = d.reason) for d in prep.S.dropped]
    res = (status = :ok, outcome = sel.outcome, record = nothing, tables = nothing, rows = prep.rows, dropped = dropped, sets = prep.S,
           relabel = (rows = prep.rel.rows, changes = prep.rel.changes, missing = prep.rel.missing), verdicts = prep.verdicts, selection = sel,
           star = pp.star, matched = pp.matched, c7 = pp.c7, optimism = pp.optimism, hb = pp.hb, hc = pp.hc, rab = pp.rab, f1 = p3_f1_table(sel))
    fp = _p3_fingerprints(i, c2, c1, ph, raw)
    program = merge(P3_PROGRAM, (module_version = string(package_version()),))
    rec = p3_record(res, q, fp, program)
    text = p3_json(rec)
    sha = bytes2hex(sha256(codeunits(text)))
    tables = p3_render(res, q, rec, fp, program, sha, res.f1)
    for f in files
        isfile(f) && bytes2hex(sha256(read(f))) == raw[f] || _fail("input file $f changed during the run")
    end
    mkpath(i.out)
    rp, tp = joinpath(i.out, P3_SLUG * "-calibration-record.json"), joinpath(i.out, P3_SLUG * "-calibration-tables.md")
    write(rp, text); write(tp, tables)
    return merge(res, (record = rp, tables = tp))
end
