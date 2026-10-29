# D89c input refusals (E36 a to p), D89d reader and identity, D89r read-only after every refusal. Spec: analysis-spec.md, Amendment D89, E24, E25, E36, I12.
# Each malformed input is made by editing a copy of the default study in place (and restoring it afterwards): the cause is named, the output folder is not created,
# every input has the same bytes after the refusal.

# the preregistration with its hash rewritten (the SELF rule) and the manifest and chunk files following, so that the edit is the only defect
function p2_edit_prereg(c, f)
    text = f(read(c.prereg, String))
    prep(t) = replace(replace(t, "\r\n" => "\n"), r"(?m)^(\| SHA-256 of this file \(current version\) \|)[^\n|]*\|" => s"\1 SELF |")
    h = an_sha256(prep(text))
    text = replace(text, r"(?m)^(\| SHA-256 of this file \(current version\) \| )[0-9a-f]{64}" => SubstitutionString("\\g<1>" * h))
    write(c.prereg, text)
    p2_edit_manifest(c, m -> m["prereg"] = Dict("header" => h, "sha256_start" => h, "sha256_end" => h))
    for fl in filter(f -> endswith(f, ".h5"), readdir(c.store; join = true))
        p2_set_attr(fl, "", "prereg_sha256", h)
    end
    return h
end

# rewrite one chunk file with the real writer after a change of its orbit records (the old file is moved aside by the writer: the caller cleans the store)
function p2_rewrite(c, i, f)
    path = p2_chunk(c, i)
    ch = AX.read_chunk(path)
    AX.write_chunk(path, Dict{String, Any}(ch.attrs), AX.OrbitRecord[f(o) for o in ch.orbits if f(o) !== nothing])
end
# a record without the rows `drop` (indices into t): t, state, winding, exponents, diagnostics and the series of the indicators
function p2_rows_off(o, drop)
    keep = setdiff(1:length(o.t), drop)
    ind = Dict{String, NamedTuple}(id => merge(t, t.series === nothing ? (;) : (series = t.series[keep, :],)) for (id, t) in o.indicators)
    return AX.OrbitRecord(o.id, o.attrs, o.t[keep], o.state[keep, :], o.winding === nothing ? nothing : o.winding[keep, :], Dict(m => v[keep, :] for (m, v) in o.lyap),
                          Dict(m => v[keep] for (m, v) in o.diagnostics), o.thresholds, o.first_exceedance, o.cost, o.failure, ind)
end
function p2_rows_off_final(o, id, rowlabel)
    t = o.indicators[id]
    keep = [r for r in eachindex(t.rows) if t.rows[r] != rowlabel]
    ind = copy(o.indicators)
    ind[id] = merge(t, (final = t.final[keep, :], rows = t.rows[keep], flags = t.flags[keep]))
    return AX.OrbitRecord(o.id, o.attrs, o.t, o.state, o.winding, o.lyap, o.diagnostics, o.thresholds, o.first_exceedance, o.cost, o.failure, ind)
end

# the malformed inputs: (tag, words the message must contain, mutation) for a copy `st` of the default study
function p2_refusal_cases(st)
    chunk(i) = p2_chunk(st, i)
    h64(c) = string(c)^64
    param_edit(f) = (d = TOML.parsefile(st.params); f(d); open(io -> TOML.print(io, d), st.params, "w"))
    p1m(f) = (m = AX.read_json(st.p1.manifest); f(m); AX._write_json(st.p1.manifest, m))
    return Tuple{String, Vector{String}, Function}[
        # a, b: the study slug and the preregistration
        ("slug-absent", ["slug"], () -> p2_edit_prereg(st, t -> replace(t, r"(?m)^\| Study slug \|[^\n]*\n" => ""))),
        ("slug-several", ["slug"], () -> p2_edit_prereg(st, t -> replace(t, "| Study slug | p2-indicators |\n" => "| Study slug | p2-indicators |\n| Study slug | p2-indicators |\n"))),
        ("slug-invalid", ["slug"], () -> p2_edit_prereg(st, t -> replace(t, "| Study slug | p2-indicators |" => "| Study slug | P2 Indicators |"))),
        ("hash-row-absent", ["header"], () -> write(st.prereg, replace(read(st.prereg, String), r"(?m)^\| SHA-256 of this file[^\n]*\n" => ""))),
        ("hash-row-repeated", ["header"], () -> write(st.prereg, replace(read(st.prereg, String), "| Study slug | p2-indicators |\n" => "| Study slug | p2-indicators |\n| SHA-256 of this file (current version) | " * "a"^64 * " |\n"))),
        ("hash-wrong", ["preregistration"], () -> write(st.prereg, read(st.prereg, String) * "x")),
        ("manifest-prereg-header", ["prereg"], () -> p2_edit_manifest(st, m -> m["prereg"]["header"] = h64('f'))),
        ("manifest-prereg-start", ["prereg"], () -> p2_edit_manifest(st, m -> m["prereg"]["sha256_start"] = h64('f'))),
        ("manifest-prereg-end", ["prereg"], () -> p2_edit_manifest(st, m -> m["prereg"]["sha256_end"] = h64('f'))),
        ("chunk-prereg", ["prereg"], () -> p2_set_attr(chunk(5), "", "prereg_sha256", h64('d'))),
        # c: the manifest
        ("manifest-schema", ["schema_version"], () -> p2_edit_manifest(st, m -> m["schema_version"] = 3)),
        ("status", ["status"], () -> p2_edit_manifest(st, m -> m["status"] = "running")),
        ("chunk-failed", ["status", "chunk"], () -> p2_edit_manifest(st, m -> m["chunks"][3]["status"] = "failed")),
        ("chunk-pending", ["status", "chunk"], () -> p2_edit_manifest(st, m -> m["chunks"][3]["status"] = "pending")),
        ("git-dirty", ["clean"], () -> p2_edit_manifest(st, m -> m["git"]["clean"] = false)),
        ("session-dirty", ["session"], () -> p2_edit_manifest(st, m -> m["sessions"][1]["clean"] = false)),
        ("control-stops", ["control"], () -> p2_edit_manifest(st, m -> m["control_stops"] = Any[Dict("chunk" => 1)])),
        ("chunk-count", ["chunk"], () -> p2_edit_manifest(st, m -> m["plan_counts"]["chunks"] = 34)),
        ("item-count", ["items"], () -> p2_edit_manifest(st, m -> m["plan_counts"]["items"] = 465)),
        # d, e: the plan and the panel
        ("plan-manifest", ["plan"], () -> p2_edit_manifest(st, m -> m["plan_sha256"] = h64('f'))),
        ("plan-chunk", ["plan"], () -> p2_set_attr(chunk(5), "", "plan_sha256", h64('f'))),
        ("panel-file", ["panel"], () -> write(st.panel, read(st.panel, String) * "\n")),
        ("panel-start", ["panel"], () -> p2_edit_manifest(st, m -> m["panel"]["sha256_start"] = h64('f'))),
        ("panel-end", ["panel"], () -> p2_edit_manifest(st, m -> m["panel"]["sha256_end"] = h64('f'))),
        ("panel-chunk", ["panel"], () -> p2_set_attr(chunk(5), "", "panel_sha256", h64('f'))),
        ("panel-build", ["build"], () -> p2_edit_manifest(st, m -> m["panel"]["build_hash"] = h64('f'))),
        # f: the chunk files
        ("chunk-missing", ["chunk-", "missing"], () -> rm(chunk(5))),
        ("pilot", ["pilot"], () -> p2_set_attr(chunk(5), "", "pilot", true)),
        ("incomplete", ["complete"], () -> p2_set_attr(chunk(5), "", "complete", false)),
        ("chunk-schema", ["schema_version"], () -> p2_set_attr(chunk(5), "", "schema_version", 1)),
        ("chunk-hash", ["chunk_hash"], () -> p2_set_attr(chunk(5), "", "chunk_hash", h64('d'))),
        ("extra-file", ["chunk-09999"], () -> cp(chunk(5), joinpath(st.store, "chunk-09999.h5"))),
        # g: items and omissions
        ("omission-not-listed", ["not_applicable"], () -> p2_edit_manifest(st, m -> delete!(m["chunks"][10], "not_applicable"))),
        ("omission-has-item", ["not_applicable"], () -> p2_edit_manifest(st, m -> m["chunks"][10]["not_applicable"] = Dict("fq1" => [20, 5]))),
        ("item-missing", ["pa1", "3"], () -> (p2_rewrite(st, 14, o -> o.id == 3 ? nothing : o); p2_edit_manifest(st, m -> m["plan_counts"]["items"] -= 1))),
        # h: identity
        ("grid-params", ["params", "21"], () -> p2_set_attr(chunk(3), "orbits/21", "params", "0.75")),
        ("grid-ic", ["ic", "22"], () -> p2_set_attr(chunk(3), "orbits/22", "ic", "0.1,0.2")),
        ("panel-ic", ["ic", "3"], () -> p2_set_attr(chunk(1), "orbits/3", "ic", "0.5,0.9")),
        ("panel-params", ["params", "4"], () -> p2_set_attr(chunk(1), "orbits/4", "params", "0.6")),
        # i: recorded steps and window rows
        ("step-missing", ["step"], () -> p2_rewrite(st, 2, o -> o.id == 7 ? p2_rows_off(o, [searchsortedfirst(P2_GA, 704)]) : o)),
        ("window-missing", ["window"], () -> p2_rewrite(st, 10, o -> o.id == 7 ? p2_rows_off_final(o, "fq1", "W=1280") : o)),
        # j: non-finite values without their flag
        ("nan-frequency", ["non-finite"], () -> p2_set_value(chunk(10), "orbits/3/indicators/fq1/final", (1, 11), NaN)),
        ("nan-series", ["non-finite"], () -> p2_set_value(chunk(2), "orbits/3/indicators/tg/series", (3, 1), NaN)),
        ("incomplete-phase", ["complete"], () -> p2_set_value(chunk(14), "orbits/3/indicators/pa1/scalars", (3,), 0.0)),
        # k: the drift flag against the stored points
        ("drift-unflagged", ["drift"], () -> p2_set_value(chunk(1), "orbits/6/state", (3, 2), 2e6)),
        ("drift-spurious", ["drift"], () -> p2_set_attr(chunk(1), "orbits/6", "flags", ["drift"])),
        # l, n, p: the table, the control deviation, the flag counts
        ("param-missing", ["budgets"], () -> param_edit(d -> delete!(d, "budgets"))),
        ("control-deviation", ["control_deviation"], () -> p2_set_attr(chunk(1), "orbits/17", "control_deviation", 0.5)),
        ("flag-counts", ["flag"], () -> p2_edit_manifest(st, m -> m["chunks"][10]["indicator_flags"]["below_floor"] += 1)),
        # m: the P1 reference inputs
        ("p1-prereg", ["preregistration"], () -> write(st.p1.prereg, read(st.p1.prereg, String) * "x")),
        ("p1-status", ["status"], () -> p1m(m -> m["status"] = "running")),
        ("p1-panel", ["panel"], () -> p1m(m -> m["panel"]["sha256_end"] = h64('f'))),
        ("p1-chunk-missing", ["chunk"], () -> rm(joinpath(st.p1.store, "chunk-00004.h5"))),
        ("p1-control-defect", ["control"],
         () -> p2_set_value(joinpath(st.p1.store, "chunk-00008.h5"), "orbits/17/lyap/L9", (length(an_grid(st.p1.p.N_1)), 1), 0.69 + 6e-8)),
    ]
end

@testset "D89c more than 40 malformed inputs are refused before any output, naming the cause" begin
    base = p2_copy(p2_cached("default"), "refuse-base")
    cases = p2_refusal_cases(base)
    @test length(cases) >= 40 && allunique(c[1] for c in cases)
    failed = String[]
    for (i, (tag, words, _)) in enumerate(cases)
        st = p2_copy(base, "refuse-$tag")                                    # a fresh copy of the study for every case
        (_, _, mutate) = p2_refusal_cases(st)[i]
        pre = p2_fingerprint(st)
        mutate()
        mid = p2_fingerprint(st)
        mid != pre || @warn "refusal case $tag: the mutation changed no file"            # a case that edits nothing would test nothing
        ok = mid != pre && an_rejects(() -> p2_run(st), words)
        ok && !ispath(st.out) && p2_fingerprint(st) == mid || push!(failed, tag)          # D89r: nothing written, no input changed by the refused run
        GC.gc()
        try rm(st.dir; recursive = true, force = true) catch end
    end
    @test isempty(failed)
    isempty(failed) || println("refusal cases not refused as expected: ", failed)
    # nothing under data/ was written, and the unmutated study is accepted (also after the refusals)
    @test p2_datatree_hash() == p2_datatree_hash()
    r = p2_run(base; threads = 2)
    @test r.status === :ok
end


@testset "D89c things that are not refusals: overrides in the manifest are listed, the stale area is ignored" begin
    st = p2_copy(p2_cached("default"), "notrefused")
    p2_edit_manifest(st, m -> m["overrides"] = Dict{String, Any}("backend" => "cpu", "dirty" => "no"))
    mkpath(joinpath(st.store, "stale")); cp(p2_chunk(st, 5), joinpath(st.store, "stale", "chunk-00005.h5"))
    r = p2_run(st)
    @test r.status === :ok && any(d -> occursin("backend", d) && occursin("overrides", lowercase(d)), r.deviations)
    text = read(r.tables, String)
    @test occursin("backend", first(split(text[findfirst("## Deviations", text).start:end], "\n## ")))
end

@testset "D89d reader and identity: every series, table, flag word and attribute of the 33 chunks, both axis orders" begin
    st = p2_cached("default")
    items = AN.p2_read_store(st.manifest, st.store)
    @test length(items) == st.items == 466
    bad = String[]
    for (i, sp) in enumerate(st.specs)
        ch = AX.read_chunk(p2_chunk(st, i))
        dir = sp.role in (:L1, :L5, :tg, :st) ? sp.direction : nothing
        for o in ch.orbits
            it = AN.p2_find_item(items; orbit = o.id, id = sp.id, precision = sp.precision, horizon = sp.horizon, tau = 1, direction = dir)
            if it === nothing
                push!(bad, "chunk $i orbit $(o.id): not found"); continue
            end
            it.chunk == i && it.t == o.t && it.state == o.state && it.winding == o.winding && it.flags == o.attrs["flags"] && it.arms == ch.attrs["arms"] &&
                it.attrs["design_class"] == o.attrs["design_class"] && it.attrs["ic"] == o.attrs["ic"] && it.mode == "wrapped" || push!(bad, "chunk $i orbit $(o.id): common fields")
            get(it.attrs, "control_scope", nothing) == get(o.attrs, "control_scope", nothing) && get(it.attrs, "control_deviation", nothing) == get(o.attrs, "control_deviation", nothing) ||
                push!(bad, "chunk $i orbit $(o.id): control attributes")
            if haskey(o.lyap, sp.id)
                it.lam == o.lyap[sp.id][:, 1] || push!(bad, "chunk $i orbit $(o.id): first exponent")
            else
                it.lam === nothing || push!(bad, "chunk $i orbit $(o.id): lam of an indicator item")
                t = o.indicators[sp.id]
                same(a, b) = a === nothing ? b === nothing : (b !== nothing && isequal(a, b))
                (same(it.series, t.series) && it.columns == t.columns && same(it.final, t.final) && it.rows == t.rows && it.fcolumns == t.fcolumns && same(it.flagwords, t.flags) &&
                 it.bits == t.bits && same(it.scalars, t.scalars) && it.scolumns == t.scolumns && same(it.range, t.range)) || push!(bad, "chunk $i orbit $(o.id): tables")
                # the Julia shapes: recorded step first for a series, the rows first for a table; n differs from the number of columns
                t.series === nothing || size(it.series) == (length(it.t), length(it.columns)) || push!(bad, "chunk $i orbit $(o.id): series shape")
                t.final === nothing || size(it.final) == (length(it.rows), length(it.fcolumns)) || push!(bad, "chunk $i orbit $(o.id): table shape")
            end
        end
    end
    @test isempty(bad)
    # the Julia shapes of the schema ("Axis order"): the recorded step first, the rows first; the number of steps differs from the number of columns
    n = length(P2_GA)
    s1 = AN.p2_find_item(items; orbit = 1, id = "tg", precision = "Float64", horizon = st.p.N_1, tau = 1, direction = "pi/4")
    s2 = AN.p2_find_item(items; orbit = 1, id = "bk1", precision = "Float64", horizon = st.p.N_c, tau = 1, direction = nothing)
    @test size(s1.series) == (n, 10) && n != 10 && size(s1.state) == (n, 2) && size(s2.final) == (40, 8) && s2.rows[1] == "W=20 cos(1 x[1])"
    # the flag words: the bit of the j-th name has the value 2^j (j = 0 for the first name)
    fq = AN.p2_find_item(items; orbit = 3, id = "fq1", precision = "Float64", horizon = st.p.N_c, tau = 1, direction = nothing)
    @test fq.bits[1] == "below_floor" && fq.flagwords[1] == 513 && AN.p2_flag(fq, 1, "below_floor") && AN.p2_flag(fq, 1, "refined1") && !AN.p2_flag(fq, 1, "degenerate") && !AN.p2_flag(fq, 1, "mirror1")
    @test all(r -> all(j -> AN.p2_flag(fq, r, fq.bits[j]) == ((fq.flagwords[r] & (Int64(1) << (j - 1))) != 0), eachindex(fq.bits)), eachindex(fq.flagwords))
    d15 = AN.p2_find_item(items; orbit = 15, id = "fq1", precision = "Float64", horizon = st.p.N_c, tau = 1, direction = nothing)
    @test d15.flagwords[1] == 2 + 1024 && AN.p2_flag(d15, 1, "degenerate") && AN.p2_flag(d15, 1, "refined2") && !AN.p2_flag(d15, 1, "below_floor")
    # the grid: label, K index, i, j, K and start texts from the [[run.h4]] tables, in the order of the tables (ids 21 to 44)
    gs = AN.p2_grid_starts(st.prereg)
    @test length(gs) == 24 && [g.id for g in gs] == collect(21:44) && all(((g, h),) -> (g.label, g.kindex, g.i, g.j, g.K, g.x0, g.p0) == (h.label, h.kindex, h.i, h.j, h.K, h.x0, h.p0), zip(gs, st.h4))
    @test gs[1].label == "grid-k1-i00-j00" && gs[13].kindex == 2 && gs[13].label == "grid-k2-i00-j00" && gs[24].label == "grid-k2-i02-j03"
    # the orbit attributes of every grid item hold the K text and the start texts of its table (neither the K index nor the label is stored)
    g = AN.p2_find_item(items; orbit = 33, id = "L1", precision = "Float64", horizon = st.p.N_1, tau = 1, direction = "pi/4")
    @test g.attrs["params"] == gs[13].K && g.attrs["ic"] == "$(gs[13].x0),$(gs[13].p0)" && !haskey(g.attrs, "label") && !haskey(g.attrs, "k_index")
    # an item that does not exist is not found
    @test AN.p2_find_item(items; orbit = 1, id = "L1", precision = "Float64", horizon = 7, tau = 1, direction = "pi/4") === nothing
end
