# D93c input refusals (E48 a to k), D93p read-only after every refusal. Spec: analysis-spec.md, Amendment D93, E46 to E48, I30, criteria D93c, D93p (S).
# Each malformed input is made by editing a fresh copy of the default study: the cause is named, the output folder is not created, every input has the same bytes after the refusal.
# E48 d comes before E48 e: a mutation of the p2 manifest is followed by a resynchronisation of the SHA-256 in the p3 table (`p3_resync!`), so that the mutation is the only defect.
# Items the command does not read (st, pr, the direction items, the sensitivity items) are not checked: the second testset pins that nothing is refused there.

function p3_refusal_cases(st)
    chunk(i) = p3_chunk(st, i)
    h64(c) = string(c)^64
    ptab(f) = (d = TOML.parsefile(st.params); f(d); open(io -> TOML.print(io, d), st.params, "w"))
    p2tab(f) = (d = TOML.parsefile(st.s2.params); f(d); open(io -> TOML.print(io, d), st.s2.params, "w"))
    p1m(f) = (m = AX.read_json(st.p1.manifest); f(m); AX._write_json(st.p1.manifest, m))
    man(f) = (p3_edit_manifest(st, f); p3_resync!(st))
    return Tuple{String, Vector{String}, Function}[
        # a, b: the slug and the preregistration of the study
        ("slug-absent", ["slug"], () -> p3_edit_prereg(st, t -> replace(t, r"(?m)^\| Study slug \|[^\n]*\n" => ""))),
        ("slug-several", ["slug"], () -> p3_edit_prereg(st, t -> replace(t, "| Study slug | p3-abstention |\n" => "| Study slug | p3-abstention |\n| Study slug | p3-abstention |\n"))),
        ("slug-invalid", ["slug"], () -> p3_edit_prereg(st, t -> replace(t, "| Study slug | p3-abstention |" => "| Study slug | P3 Abstention |"))),
        ("hash-row-absent", ["header"], () -> write(st.prereg, replace(read(st.prereg, String), r"(?m)^\| SHA-256 of this file[^\n]*\n" => ""))),
        ("hash-row-repeated", ["header"], () -> write(st.prereg, replace(read(st.prereg, String), "| Study slug | p3-abstention |\n" => "| Study slug | p3-abstention |\n| SHA-256 of this file (current version) | " * "a"^64 * " |\n"))),
        ("hash-wrong", ["preregistration"], () -> write(st.prereg, read(st.prereg, String) * "x")),
        # c: the parameter table
        ("param-missing", ["budgets"], () -> ptab(d -> delete!(d, "budgets"))),
        ("param-unknown", ["bogus_key"], () -> ptab(d -> d["bogus_key"] = Dict("value" => 1, "section" => "x", "quote" => "1", "text" => "1"))),
        ("param-quote", ["min_set"], () -> ptab(d -> d["min_set"]["quote"] = "text that is not in the section")),
        ("param-text", ["cost_tg"], () -> ptab(d -> d["cost_tg"]["text"] = "text that is not in the quote")),
        ("param-mirror", ["n_s"], () -> p3_set_q!(st, merge(st.q, (; n_s = 801)))),
        ("param-range", ["e_cal"], () -> p3_set_q!(st, merge(st.q, (; e_cal = 1.5)))),
        ("param-voters", ["voter_bits"], () -> p3_set_q!(st, merge(st.q, (; voter_bits = reverse(st.q.voter_bits))))),
        ("p2-param-missing", ["ratio_threshold"], () -> p2tab(d -> delete!(d, "ratio_threshold"))),
        # d: the p2 manifest against p2_manifest_sha256 (the p2 manifest edited, the p3 table not resynchronised)
        ("p2-manifest-hash", ["p2_manifest_sha256"], () -> p3_edit_manifest(st, m -> m["reasons"] = String["edited after the table"])),
        # e: the p2 inputs (D89 E36 b to f and h to k, for the items read)
        ("p2-study", ["study"], () -> man(m -> m["study"] = "p3-abstention")),
        ("p2-prereg-header", ["prereg"], () -> man(m -> m["prereg"]["header"] = h64('f'))),
        ("p2-prereg-start", ["prereg"], () -> man(m -> m["prereg"]["sha256_start"] = h64('f'))),
        ("p2-prereg-end", ["prereg"], () -> man(m -> m["prereg"]["sha256_end"] = h64('f'))),
        ("p2-chunk-prereg", ["prereg"], () -> p2_set_attr(chunk(5), "", "prereg_sha256", h64('d'))),
        ("p2-manifest-schema", ["schema_version"], () -> man(m -> m["schema_version"] = 3)),
        ("p2-status", ["status"], () -> man(m -> m["status"] = "running")),
        ("p2-chunk-failed", ["status", "chunk"], () -> man(m -> m["chunks"][3]["status"] = "failed")),
        ("p2-git-dirty", ["clean"], () -> man(m -> m["git"]["clean"] = false)),
        ("p2-control-stops", ["control"], () -> man(m -> m["control_stops"] = Any[Dict("chunk" => 1)])),
        ("p2-plan-manifest", ["plan"], () -> man(m -> m["plan_sha256"] = h64('f'))),
        ("p2-panel-file", ["panel"], () -> (write(st.s2.panel, read(st.s2.panel, String) * "\n"))),
        ("p2-panel-manifest", ["panel"], () -> man(m -> m["panel"]["sha256_end"] = h64('f'))),
        ("p2-chunk-missing", ["chunk-", "missing"], () -> rm(chunk(5))),
        ("p2-pilot", ["pilot"], () -> p2_set_attr(chunk(5), "", "pilot", true)),
        ("p2-incomplete", ["complete"], () -> p2_set_attr(chunk(5), "", "complete", false)),
        ("p2-chunk-hash", ["chunk_hash"], () -> p2_set_attr(chunk(5), "", "chunk_hash", h64('d'))),
        ("p2-nan-series", ["non-finite"], () -> p2_set_value(chunk(2), "orbits/3/indicators/tg/series", (3, 1), NaN)),
        ("p2-nan-frequency", ["non-finite"], () -> p2_set_value(chunk(10), "orbits/3/indicators/fq1/final", (1, 11), NaN)),
        ("p2-drift-unflagged", ["drift"], () -> p2_set_value(chunk(1), "orbits/6/state", (3, 2), 2e6)),
        ("p2-drift-spurious", ["drift"], () -> p2_set_attr(chunk(1), "orbits/6", "flags", ["drift"])),
        # f: the P1 reference inputs (E48 f), also when only the change counts need them (P46)
        ("p1-prereg", ["preregistration"], () -> write(st.p1.prereg, read(st.p1.prereg, String) * "x")),
        ("p1-status", ["status"], () -> p1m(m -> m["status"] = "running")),
        ("p1-chunk-failed", ["status", "chunk"], () -> p1m(m -> m["chunks"][5]["status"] = "failed")),
        ("p1-chunk-missing", ["chunk"], () -> rm(p3_p1chunk(st, 5))),
        ("p1-chunk-hash", ["hash"], () -> p2_set_attr(p3_p1chunk(st, 5), "", "chunk_hash", h64('d'))),
        ("p1-study", ["study"], () -> p1m(m -> m["study"] = "p2-indicators")),
        ("p1-panel", ["panel"], () -> p1m(m -> m["panel"]["sha256_end"] = h64('f'))),
        # g: identity of the stored texts with the panel file and the grid tables (p2 items and P1 items read)
        ("id-grid-params", ["params", "21"], () -> p2_set_attr(chunk(3), "orbits/21", "params", "0.75")),
        ("id-grid-ic", ["ic", "22"], () -> p2_set_attr(chunk(3), "orbits/22", "ic", "0.1,0.2")),
        ("id-panel-ic", ["ic", "3"], () -> p2_set_attr(chunk(1), "orbits/3", "ic", "0.5,0.9")),
        ("id-panel-params", ["params", "4"], () -> p2_set_attr(chunk(1), "orbits/4", "params", "0.6")),
        ("id-p1-ic", ["ic", "3"], () -> p2_set_attr(p3_p1chunk(st, 5), "orbits/3", "ic", "0.5,0.9")),
        # h: an item of E51 missing for an orbit of Cal (or of the gate) that the manifest does not list in not_applicable; a listed item that exists
        ("item-missing-cal", ["pa1", "3"], () -> (p3_rewrite(st, 14, o -> o.id == 3 ? nothing : o); p3_edit_manifest(st, m -> m["plan_counts"]["items"] -= 1); p3_resync!(st))),
        ("item-missing-gate", ["pa1", "9"], () -> (p3_rewrite(st, 14, o -> o.id == 9 ? nothing : o); p3_edit_manifest(st, m -> m["plan_counts"]["items"] -= 1); p3_resync!(st))),
        ("item-missing-grid", ["L1", "22"], () -> (p3_rewrite(st, 3, o -> o.id == 22 ? nothing : o); p3_edit_manifest(st, m -> m["plan_counts"]["items"] -= 1); p3_resync!(st))),
        ("item-listed-exists", ["not_applicable"], () -> man(m -> m["chunks"][10]["not_applicable"] = Dict("fq1" => [20, 5]))),
        # i: a recorded step or a window row that E51 needs is missing
        ("step-missing", ["step"], () -> p3_rewrite(st, 2, o -> o.id == 7 ? p2_rows_off(o, [searchsortedfirst(P2_GA, 704)]) : o)),
        ("window-missing", ["window"], () -> p3_rewrite(st, 10, o -> o.id == 7 ? p2_rows_off_final(o, "fq1", "W=1280") : o)),
    ]
end

@testset "D93c more than 30 malformed inputs are refused before any output, naming the cause; D93p inputs unchanged, output folder not created" begin
    base = p3_copy(p3_cached("default"), "refuse-base")
    cases = p3_refusal_cases(base)
    @test length(cases) >= 30 && allunique(c[1] for c in cases)
    failed = String[]
    for (i, (tag, words, _)) in enumerate(cases)
        (P3_FULL || tag in P3_CASE_SAMPLE) || continue                          # the other cases: verify/analysis-p3.jl (D38)
        st = p3_copy(base, "refuse-$tag")                                    # a fresh copy of the study for every case
        (_, _, mutate) = p3_refusal_cases(st)[i]
        pre = p3_fingerprint(st)
        mutate()
        mid = p3_fingerprint(st)
        mid != pre || @warn "refusal case $tag: the mutation changed no file"            # a case that edits nothing would test nothing
        ok = mid != pre && an_rejects(() -> p3_run(st), words)
        ok && !ispath(st.out) && p3_fingerprint(st) == mid || push!(failed, tag)          # I30: nothing written, no input changed by the refused run
        GC.gc()
        try rm(st.dir; recursive = true, force = true) catch end
    end
    @test isempty(failed)
    isempty(failed) || println("refusal cases not refused as expected: ", failed)
    # nothing under data/ was written; the unmutated study is accepted (also after the refusals)
    h0 = p2_datatree_hash()
    r = p3_run(base; threads = 2)
    @test r.status === :ok && p2_datatree_hash() == h0
end

@testset "D93c things that are not refusals: items the command does not read, overrides in the manifest, a missing item that only the relabelling needs" begin
    st0 = p3_cached("default")
    ref = p3_run(p3_copy(st0, "notrefused-ref"); threads = 1)
    # (a) the items that E49 and E51 do not read: pr (RLI), st (stretching), the direction and sensitivity items, the Double64 items of SALI and so on of orbits outside S_L:
    # values written as NaN in them change nothing (no refusal, the same record)
    st = p3_copy(st0, "notrefused-unread")
    p2_set_value(p3_chunk(st, 9), "orbits/21/indicators/pr/series", (3, 2), NaN)           # pr, grid
    p2_set_value(p3_chunk(st, 6), "orbits/1/indicators/st/series", (3, 1), NaN)            # st, panel
    p2_set_value(p3_chunk(st, 7), "orbits/21/indicators/st/series", (3, 1), NaN)           # st, grid
    p2_set_value(p3_chunk(st, 18), "orbits/7/lyap/L1", (3, 1), NaN)                        # a direction item (dindex 3) of S_D
    p2_set_value(p3_chunk(st, 28), "orbits/1/indicators/pa_planar/final", (1, 1), NaN)     # a sensitivity item of S_L
    r = p3_run(st; threads = 1)
    @test r.status === :ok && read(r.record, String) == read(ref.record, String) && read(r.tables, String) == read(ref.tables, String)
    # (b) the overrides of the manifest are recorded, not refused (D89 E36)
    st2 = p3_copy(st0, "notrefused-overrides")
    p3_edit_manifest(st2, m -> m["overrides"] = Dict{String, Any}("backend" => "cpu", "dirty" => "no")); p3_resync!(st2)
    @test p3_run(st2).status === :ok
    # (c) every cause of missing data of C0 (E49 e): the orbit leaves Cal, is listed as dropped with the reason "missing", is counted in relabel.missing by source, K index and item, and nothing
    # is refused (a missing item that only the relabelling needs is not a refusal, D93d)
    rewrite_p1(c, i, f) = (path = p3_p1chunk(c, i); ch = AX.read_chunk(path); AX.write_chunk(path, Dict{String, Any}(ch.attrs), AX.OrbitRecord[f(o) for o in ch.orbits if f(o) !== nothing]))
    step(o, n) = p2_rows_off(o, [searchsortedfirst(o.t, n)])
    causes = [("L5 item", 5, o -> o.id == 4 ? nothing : o, 4, "L5"),
              ("L5 step N_1/10", 5, o -> o.id == 6 ? step(o, P2_N1 ÷ 10) : o, 6, "L5 step N_1/10"),
              ("L5 step N_1/2", 5, o -> o.id == 6 ? step(o, P2_N1 ÷ 2) : o, 6, "L5 step N_1/2"),
              ("L5 step N_1", 5, o -> o.id == 6 ? step(o, P2_N1) : o, 6, "L5 step N_1"),
              ("L1 item", 1, o -> o.id == 7 ? nothing : o, 7, "L1")]
    for (tag, chunk, f, id, item) in causes
        c = p3_copy(st0, "notrefused-missing-" * replace(tag, r"[^A-Za-z0-9]" => "-"))
        rewrite_p1(c, chunk, f)
        r = p3_run(c)
        @test r.status === :ok && id in [d.id for d in r.dropped] && first(d.reason for d in r.dropped if d.id == id) == "missing" && !(id in r.sets.cal) && length(r.sets.cal) == length(ref.sets.cal) - 1
        @test any(m -> m.id == id && m.source === :panel && m.kindex == Int(AN_ORBITS[id]["k_index"]) && m.item == item, r.relabel.missing) && length(r.relabel.missing) == 1
        GC.gc(); try rm(c.dir; recursive = true, force = true) catch end
    end
end
