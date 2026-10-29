# Pipeline on synthetic studies: AC1 (refusals), AC2 (reader), AC6 (control defect), AC14 to AC17 (parameters, content, determinism, read-only),
# AC18 (figures skipped), and the scenario values of E7 to E15 end to end. Spec: analysis-spec.md. Fixture: an_helpers.jl.
# The expected values follow from the closed forms of the fixture (lam = lref (1 + c / n^q)); see the test plan, "Scenarios".

const AN_RES = Dict{String, Any}()
an_res(tag, st; kw...) = get!(() -> an_run(st; kw...), AN_RES, tag)
an_chunk(st, i) = joinpath(st.store, "chunk-" * lpad(string(i), 5, '0') * ".h5")
an_find(st; method, horizon = st.p.N_1, mode = "lifted", tau = 1, precision = "Float64", dir = "pi/4", eps = nothing) =
    findfirst(s -> s.method == method && s.horizon == horizon && s.mode == mode && s.tau == tau && s.precision == precision && s.direction == dir && s.eps == eps, st.specs)
function an_set_root_attr(file, name, val)
    AX.HDF5.h5open(file, "r+") do h
        AX.HDF5.delete_attribute(h, name)
        AX.HDF5.write_attribute(h, name, val)
    end
end
function an_set_attr(file, path, name, val)
    AX.HDF5.h5open(file, "r+") do h
        g = h[path]
        AX.HDF5.delete_attribute(g, name)
        AX.HDF5.write_attribute(g, name, val)
    end
end
function an_set_value(file, path, idx, val)          # HDF5.jl indexes in the Julia shape (n, k): recorded row j, vector k is [j, k]
    AX.HDF5.h5open(file, "r+") do h
        d = h[path]
        d[idx...] = val
    end
end
function an_edit_manifest(st, f)
    m = AX.read_json(st.manifest)
    f(m)
    AX._write_json(st.manifest, m)
end
an_fingerprint(st) = [hash(read(f)) for f in vcat([st.manifest, st.prereg, st.panel, st.params], sort(filter(f -> endswith(f, ".h5"), readdir(st.store; join = true))))]
an_dec(res, k) = repr(getfield(res.decisions, k))

@testset "AC2 reader: every series, point, flag and attribute is what the writer wrote" begin
    st = an_cached("default")
    items = AN.read_store(st.manifest, st.store)
    @test length(items) == sum(length(sp.orbits) for sp in st.specs)
    bad = String[]
    for (i, sp) in enumerate(st.specs)
        ch = AX.read_chunk(an_chunk(st, i))
        for o in ch.orbits
            it = AN.find_item(items; orbit = o.id, method = sp.method, precision = sp.precision, horizon = sp.horizon, tau = sp.tau, eps = sp.eps,
                              direction = sp.direction, mode = sp.mode)
            if it === nothing
                push!(bad, "chunk $i orbit $(o.id): item not found")
                continue
            end
            it.t == o.t || push!(bad, "chunk $i orbit $(o.id): t")
            it.lam == o.lyap[sp.method][:, 1] || push!(bad, "chunk $i orbit $(o.id): lambda")
            it.state == o.state || push!(bad, "chunk $i orbit $(o.id): state")
            it.flags == o.attrs["flags"] || push!(bad, "chunk $i orbit $(o.id): flags")
            it.arms == ch.attrs["arms"] || push!(bad, "chunk $i orbit $(o.id): arms")
            it.attrs["design_class"] == o.attrs["design_class"] || push!(bad, "chunk $i orbit $(o.id): attributes")
            get(it.attrs, "control_scope", nothing) == get(o.attrs, "control_scope", nothing) || push!(bad, "chunk $i orbit $(o.id): control_scope")
            get(it.attrs, "control_deviation", nothing) == get(o.attrs, "control_deviation", nothing) || push!(bad, "chunk $i orbit $(o.id): control_deviation")
            Set(keys(it.diag)) == Set(keys(o.diagnostics)) || push!(bad, "chunk $i orbit $(o.id): diagnostics")
        end
    end
    @test isempty(bad)
    # the winding of a wrapped item is read too; a wrapped item has no lifted state
    w = AN.find_item(items; orbit = 1, method = "L2", precision = "Float64", horizon = st.p.N_1, tau = 1, eps = 1e-8, direction = "pi/4", mode = "wrapped")
    @test w !== nothing && w.eps == 1e-8 && size(w.winding) == (length(w.t), 2)
    @test AN.find_item(items; orbit = 1, method = "L5", precision = "Float64", horizon = st.p.N_1, tau = 99, eps = nothing, direction = "pi/4", mode = "lifted") === nothing
end

@testset "AC15 and the default scenario: content, sets, decisions, labels" begin
    st = an_cached("default")
    res = an_res("default", st)
    @test res.status === :ok && isfile(res.tables) && basename(res.tables) == "p1-le-benchmark-tables.md"
    @test an_files(st.out) == ["p1-le-benchmark-tables.md"]                                  # figures not requested: tables only
    @test collect(keys(res.decisions)) == [:H1, :H2, :H3, :H4, :H5]                      # decided in this order (I5)
    # sets
    @test sort(res.sets.E_ch) == [1, 2, 3, 4, 11, 12, 13, 14, 17, 18, 20] && sort(res.sets.W) == collect(1:8) && sort(res.sets.S) == [11, 12, 13, 14]
    @test isempty(res.sets.dropped) && res.sets.S_nonpositive == 0
    # H1: every E_ch orbit passes
    h1 = res.decisions.H1
    @test h1.verdict === :supported && h1.n == 11 && h1.fractions.a == 1.0 && h1.fractions.b == 1.0 && h1.fractions.c == 1.0
    # H2: R(n_H2) = 0.02 / 100 for every W orbit (L9 is the largest, L1 and the QR family the smallest), not supported
    h2 = res.decisions.H2
    @test h2.verdict === :not_supported && h2.n == 8 && isapprox(h2.median, 2e-4; rtol = 1e-9) && isapprox(h2.lo, 2e-4; rtol = 1e-9)
    # H3: every V method is inside the 5 percent band from the first step, ratio 1
    h3 = res.decisions.H3
    @test h3.verdict === :not_supported && h3.n_w == 8 && h3.n_s == 4 && h3.ratio ≈ 1.0 && h3.undefined_share == 0.0
    # H4: the h4 starts reproduce the references exactly
    h4 = res.decisions.H4
    @test h4.verdict === :supported && h4.halfpi.pass && h4.twopi.pass
    @test length(h4.crosscheck) == 2 && all(isapprox(r.float64, ref) && isapprox(r.double64, ref) for (r, ref) in zip(h4.crosscheck, (0.298, 1.172)))      # no verdict, both K
    # H5
    h5 = res.decisions.H5
    @test h5.verdict === :supported && h5.h5a.n == 4 && h5.h5a.fraction == 1.0 && h5.h5b.n == 20 && h5.h5b.fraction == 0.0 && h5.h5c.n == 5
    # labels of every panel orbit, once
    @test sort(collect(keys(res.labels))) == collect(1:20)
    lab(i) = (res.labels[i].instability, res.labels[i].structure)
    for i in vcat(1:8, 11:14); @test lab(i) == (AN_POS, "irregular"); end
    for i in (9, 10, 15, 16); @test lab(i) == (AN_NONE, "quasiperiodic"); end
    @test lab(17) == (AN_POS, "periodic") && lab(19) == (AN_NONE, "periodic")
    @test lab(18) == (AN_POS, "irregular")                        # the info control: its start left the orbit, the return residual is 2e-6
    @test lab(20) == (AN_POS, "irregular")                        # Y5 is not a periodic point of the standard map
    @test all(isempty(res.labels[i].forced) for i in 1:20)
    @test res.labels[1].ladder_same === true && res.labels[2].ladder_same === nothing        # S_L = {1, 5, 9, 11}
    # cells: every panel orbit in one cell, label counts add up (I7, I10)
    @test sum(c.orbits for c in res.cells) == 20 && all(sum(values(c.labels)) == c.orbits for c in res.cells)
    # controls table: gated controls and the info control with its deviation
    rows = Dict(r.id => r for r in res.controls)
    @test Set(keys(rows)) == Set([17, 18, 20]) && rows[17].scope == "bound" && rows[18].scope == "info"
    @test isapprox(rows[18].deviations["L5"], 0.3; atol = 1e-12) && rows[17].deviations["L9"] == 0.0
    # header: hashes of the four inputs and no clock or path
    text = read(res.tables, String)
    @test res.header.prereg_header_sha256 == st.prereg_hash && res.header.panel_sha256 == st.panel_hash
    @test res.header.manifest_sha256 == an_sha256(read(st.manifest)) && res.header.params_sha256 == an_sha256(read(st.params))
    @test !occursin(st.dir, text) && !occursin(AN_DIR, text) && !occursin(r"20\d\d-\d\d-\d\dT\d\d:\d\d", text)
    # the sections of the analysis plan, each marked confirmatory; the figure data with their digests
    for id in ("T1", "T2", "T3", "T4", "T5", "T6", "T7", "F1", "F2", "F3", "F4", "F5", "F6", "F7")
        @test occursin(Regex("(?m)^## $id[^\\n]*confirmatory"), text)
    end
    for i in 1:7
        j = findfirst("## F$i data", text)
        @test j !== nothing
        m = match(r"digest: ([0-9a-f]{64})", text[j.start:end])
        @test m !== nothing && m.captures[1] == res.digests["F$i"]
    end
    @test Set(keys(res.digests)) == Set("F$i" for i in 1:7)
    for h in ("## Header", "## Readings", "## Deviations", "## Exploratory items")
        @test occursin(h, text)
    end
    ex = text[findfirst("## Exploratory items", text).start:end]
    @test occursin("none", lowercase(first(split(ex, "\n## "))))                        # I8: the module produces no exploratory item
    rd = first(split(text[findfirst("## Readings", text).start:end], "\n## ")[1:1])
    @test count(l -> startswith(l, "- "), split(rd, '\n')) >= 8
    @test any(d -> occursin("L2", d) && occursin("wrapped", d), res.deviations) && any(d -> occursin("5e-8", d) || occursin("5.0e-8", d), res.deviations)
    @test any(d -> occursin("overrides", lowercase(d)), res.deviations)
end

@testset "AC14 weak scenario: H2 and H3 supported, H1 untouched" begin
    cover = Dict(i => Dict("L9" => (3e5, 3.0)) for i in 1:8)      # L9 differs by 3e5 / n^3: R(100) = 0.3, Delta at N_1 about 1e-7
    st = an_cached("weak"; cover = cover)
    res = an_res("weak", st)
    @test res.decisions.H1.verdict === :supported
    h2 = res.decisions.H2
    @test h2.verdict === :supported && h2.n == 8 && isapprox(h2.median, 0.3; rtol = 1e-9) && isapprox(h2.lo, 0.3; rtol = 1e-9)
    h3 = res.decisions.H3
    @test h3.verdict === :supported && h3.ratio ≈ 200.0 && h3.lo ≈ 200.0           # n_conv 200 on W (first grid step with 3e5 / n^3 <= 0.05), 1 on S
end

@testset "AC12 and AC13 disagreement tests end to end" begin
    cover = Dict(2 => Dict("L9" => 100.0))
    lamh(id, sp, steps, vals) = id == 9 && sp.precision == "BigFloat(256)" ? vals .+ 0.05 :
                                id == 5 && sp.precision == "BigFloat(256)" ? vals .+ 0.01 :
                                id == 7 && sp.method == "L5" && sp.direction == "0.3" ? fill(1e-5, length(vals)) : vals
    flagh(id, sp) = (id == 13 && sp.method == "L5" && "main" in sp.arms) ? ["diagnostics exceeded"] : String[]
    st = an_cached("force"; cover = cover, lam = lamh, flags = flagh)
    res = an_res("force", st)
    L = res.labels
    @test L[9].forced == [:ladder] && L[5].forced == [:ladder] && L[7].forced == [:direction] && L[2].forced == [:estimator] && L[13].forced == [:diagnostics]
    for i in (2, 5, 7, 9, 13)
        @test (L[i].instability, L[i].structure) == (AN_UNR, AN_UNR)
    end
    @test L[1].forced == Symbol[] && L[1].instability == AN_POS
    @test L[9].raw.instability == AN_NONE && L[9].ladder_same === false && L[5].ladder_same === true      # the raw labels differ for 9 only
    # forced orbits leave W; E_ch and S do not depend on the labels
    @test sort(res.sets.W) == [1, 3, 4, 6, 8] && 2 ∈ res.sets.E_ch && 13 ∈ res.sets.S && res.sets.S_nonpositive == 1
    # H1: orbit 2 fails H1b, 10 of 11 pass
    @test res.decisions.H1.verdict === :falsified && res.decisions.H1.passes.b == 10 && res.decisions.H1.passes.a == 11
    # H5a from the raw labels: 3 of 4, not 1.0 by construction
    h5 = res.decisions.H5
    @test h5.h5a.n == 4 && h5.h5a.fraction == 0.75 && h5.h5a.pass === false && h5.verdict === :not_supported
    @test h5.h5b.fraction == 0.25 && h5.h5b.pass === true && h5.h5c.n == 4 && h5.h5c.fraction == 0.0
    # the label table says which tests could be applied: S_L orbits carry the ladder test, S_D orbits the direction test
    @test :ladder in L[9].applicable && :direction in L[9].applicable && :ladder ∉ L[2].applicable && :direction in L[2].applicable
    @test :ladder ∉ L[3].applicable && :direction ∉ L[3].applicable && :estimator in L[3].applicable && :diagnostics in L[3].applicable
end

@testset "AC7 and AC12 flags: drift, excluded, diagnostics" begin
    flagh(id, sp) = id == 3 && sp.mode == "lifted" ? ["drift"] : id == 4 ? ["excluded"] : (id == 5 && sp.method == "L5" && "main" in sp.arms) ? ["diagnostics exceeded"] : String[]
    stateh(id, st) = id == 3 ? (s = copy(st); s[end, 2] = 2e6; s) : st
    st = an_cached("flags"; flags = flagh, state = stateh)
    res = an_res("flags", st)
    @test Dict(d.id => d.reason for d in res.sets.dropped) == Dict(3 => "drift", 4 => "excluded")
    @test sort(res.sets.E_ch) == [1, 2, 11, 12, 13, 14, 17, 18, 20] && sort(res.sets.W) == [1, 2, 6, 7, 8]
    @test res.labels[4].forced == [:excluded] && (res.labels[4].instability, res.labels[4].structure) == (AN_UNR, AN_UNR)
    @test res.labels[5].forced == [:diagnostics]
    @test res.labels[3].forced == Symbol[]                                       # drift removes from the sets, it does not change the label
    @test res.decisions.H1.n == 9 && res.decisions.H5.h5b.n == 19                # excluded orbits are left out of every fraction
    @test sum(c.orbits for c in res.cells) == 20                                  # but they stay in the cell counts (I7)
end

@testset "AC6 control defect stops the analysis" begin
    lamh(id, sp, steps, vals) = id == 17 && sp.method == "L9" && "main" in sp.arms ? vals .+ 6e-8 : vals
    st = an_cached("defect"; lam = lamh)
    res = an_res("defect", st)
    @test res.status === :control_defect && res.decisions === nothing
    @test length(res.defects) == 1 && res.defects[1].id == 17 && res.defects[1].method == "L9"
    @test an_files(st.out) == ["p1-le-benchmark-defect.md"]                             # a defect report and nothing else
    rep = read(joinpath(st.out, "p1-le-benchmark-defect.md"), String)
    @test occursin("17", rep) && occursin("L9", rep) && !occursin("## T2", rep)
end

@testset "AC1 refusals: malformed inputs are refused before any output" begin
    st = an_cached("default")
    iL5 = an_find(st; method = "L5")
    iL1 = an_find(st; method = "L1")
    rewrite_params(st2, f) = (d = TOML.parsefile(st2.params); f(d); open(io -> TOML.print(io, d), st2.params, "w"))
    cases = [
        ("status", ["status"], s -> an_edit_manifest(s, m -> m["status"] = "running")),
        ("chunk-status", ["status", "chunk"], s -> an_edit_manifest(s, m -> m["chunks"][3]["status"] = "pending")),
        ("prereg-text", ["preregistration"], s -> write(s.prereg, read(s.prereg, String) * "x")),
        ("prereg-manifest", ["preregistration"], s -> an_edit_manifest(s, m -> m["prereg"]["header"] = "f"^64)),
        ("panel", ["panel"], s -> write(s.panel, read(s.panel, String) * "\n")),
        ("chunk-missing", ["chunk-", "missing"], s -> rm(an_chunk(s, 5))),
        ("pilot", ["pilot"], s -> an_set_root_attr(an_chunk(s, 5), "pilot", true)),
        ("incomplete", ["complete"], s -> an_set_root_attr(an_chunk(s, 5), "complete", false)),
        ("schema", ["schema_version"], s -> an_set_root_attr(an_chunk(s, 5), "schema_version", 1)),
        ("manifest-schema", ["schema_version"], s -> an_edit_manifest(s, m -> m["schema_version"] = 3)),
        ("chunk-hash", ["chunk_hash"], s -> an_set_root_attr(an_chunk(s, 5), "chunk_hash", "d"^64)),
        ("extra-file", ["chunk-09999"], s -> cp(an_chunk(s, 5), joinpath(s.store, "chunk-09999.h5"))),
        ("nan", ["non-finite"], s -> an_set_value(an_chunk(s, iL5), "orbits/1/lyap/L5", (3, 1), NaN)),
        ("drift-unflagged", ["drift"], s -> an_set_value(an_chunk(s, iL5), "orbits/6/state", (length(an_grid(s.p.N_1)), 2), 2e6)),
        ("control-deviation", ["control_deviation"], s -> an_set_attr(an_chunk(s, iL1), "orbits/17", "control_deviation", 0.5)),
        ("param-missing", ["n_H2"], s -> rewrite_params(s, d -> delete!(d, "n_H2"))),
        ("param-unknown", ["bogus_key"], s -> rewrite_params(s, d -> d["bogus_key"] = 1)),
    ]
    @test length(cases) >= 15
    for (tag, words, mutate) in cases
        s = an_copy(st, "refuse-$tag")
        mutate(s)
        ok = an_rejects(() -> an_run(s), words)
        ok || println("refusal case not refused as expected: ", tag)
        @test ok && !ispath(s.out)
    end
    # a store file under stale/ is ignored
    s = an_copy(st, "stale")
    mkpath(joinpath(s.store, "stale"))
    cp(an_chunk(s, 5), joinpath(s.store, "stale", "chunk-00005.h5"))
    r = an_run(s)
    @test r.status === :ok && an_dec(r, :H1) == an_dec(an_res("default", st), :H1)
end

@testset "AC17 read-only: inputs unchanged after a run and after a refusal" begin
    st = an_copy(an_cached("default"), "readonly")
    before = an_fingerprint(st)
    r = an_run(st)
    @test r.status === :ok && an_fingerprint(st) == before
    an_edit_manifest(st, m -> m["status"] = "running")
    before = an_fingerprint(st)
    @test_throws ArgumentError an_run(st; out = an_path("readonly-out2"))
    @test an_fingerprint(st) == before
end

@testset "AC16 determinism: byte-identical tables at any thread count and output place" begin
    st = an_cached("default")
    a = an_res("default", st)
    b = an_run(st; out = an_path("det-out-1"), threads = 1)
    c = an_run(st; out = an_path("det-out-4"), threads = 4)
    @test read(a.tables) == read(b.tables) == read(c.tables)
    @test a.digests == b.digests == c.digests
end

@testset "AC14 a parameter changes only what depends on it" begin
    st = an_cached("default")
    base = an_res("default", st)
    ks = (:H1, :H2, :H3, :H4, :H5)
    function changed(tag, kv...)
        pf = an_write_params(joinpath(st.dir, "params-$tag.toml"), merge(st.p, NamedTuple(kv)))
        r = an_run(st; params = pf, out = an_path("pert-$tag"))
        return Set(k for k in ks if an_dec(r, k) != an_dec(base, k))
    end
    @test changed("h4", :h4_ref_halfpi => 0.5) == Set([:H4])
    @test changed("h5a", :h5a_min => 1.1) == Set([:H5])
    @test changed("h2", :h2_median => 1e-4, :h2_interval => 1e-4) == Set([:H2])
    @test changed("h3", :h3_ratio => 0.5, :h3_interval => 0.5) == Set([:H3])
    @test changed("h1c", :h1c_fraction => 1.01) == Set([:H1])
    @test changed("minset", :min_set => 12) == Set([:H1, :H2, :H3, :H5])
end

@testset "AC18 figures: skipped without the plotting extension" begin
    st = an_cached("default")
    r = an_run(st; out = an_path("fig-out"), figures = true)
    @test r.status === :figures_skipped && isfile(r.tables)
    @test !any(f -> endswith(f, ".pdf"), an_files(an_path("fig-out")))
    @test Set(keys(r.digests)) == Set("F$i" for i in 1:7) && all(occursin(r"^[0-9a-f]{64}$", d) for d in values(r.digests))
    text = read(r.tables, String)
    @test all(occursin(d, text) for d in values(r.digests))
end
