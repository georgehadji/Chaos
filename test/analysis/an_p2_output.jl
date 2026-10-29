# D89p output content, D89q determinism, D89r read-only, D89s command, D89u figure data. Spec: analysis-spec.md, Amendment D89, E37, I12, I13, I17,
# criteria D89p, D89q, D89r, D89s, D89u. The tables are res.items.T1 .. T17 and F1 .. F6, each (columns, rows); the file holds the same tables as markdown.

const P2_DETRUNS = Dict{Symbol, Any}()
@testset "D89p T2: counts add up to the set sizes (I17), by source and by K, the step used" begin
    st = p2_cached("default"); res = p2_res("default", st); (; ref, S) = p2_refs("default", st)
    T2 = res.items.T2
    kc = p2_kfield(T2)
    cols = ["chaotic", "regular", "undecided", "not applicable"]
    sets = Dict("W2" => S.W2, "S2" => S.S2, "T_s" => S.T_s, "regular" => S.regular)
    bad = String[]
    for X in P2_X8, k in 0:9, (nm, ids) in sets
        rows = p2_rowsof(T2, "indicator" => String(X), "k" => k, "set" => nm, "source" => "all", kc => "all")
        length(rows) == 1 || (push!(bad, "$X $k $nm rows"); continue)
        sum(rows[1][p2_tcol(T2, c)] for c in cols) == length(ids) || push!(bad, "$X $k $nm sum")
        panel = p2_rowsof(T2, "indicator" => String(X), "k" => k, "set" => nm, "source" => "panel", kc => "all")[1]
        grid = p2_rowsof(T2, "indicator" => String(X), "k" => k, "set" => nm, "source" => "grid", kc => "all")[1]
        k1 = p2_rowsof(T2, "indicator" => String(X), "k" => k, "set" => nm, "source" => "all", kc => p2_klabel(T2, 0.5, 1))[1]
        k2 = p2_rowsof(T2, "indicator" => String(X), "k" => k, "set" => nm, "source" => "all", kc => p2_klabel(T2, 2 * pi, 2))[1]
        for c in cols
            j = p2_tcol(T2, c)
            (panel[j] + grid[j] == rows[1][j] && k1[j] + k2[j] == rows[1][j]) || push!(bad, "$X $k $nm by-source or by-K $c")
        end
    end
    @test isempty(bad)
    # hand rows: LE at k = 6 on W2: seven chaotic (W and W2 profiles), four regular (the sticky orbits); the step n(b_6) = 2816; RLI exists for 8 of the 11 orbits
    r = p2_rowsof(T2, "indicator" => "LE", "k" => 6, "set" => "W2", "source" => "all", kc => "all")[1]
    @test Tuple(r[p2_tcol(T2, c)] for c in cols) == (7, 4, 0, 0) && r[p2_tcol(T2, "step")] == P2_NK[7] == 2816 && r[p2_tcol(T2, "budget")] == 2560
    r = p2_rowsof(T2, "indicator" => "RLI", "k" => 6, "set" => "W2", "source" => "all", kc => "all")[1]
    @test Tuple(r[p2_tcol(T2, c)] for c in cols) == (8, 0, 0, 3) && r[p2_tcol(T2, "step")] == P2_NKB[7]
    r = p2_rowsof(T2, "indicator" => "E1", "k" => 6, "set" => "W2", "source" => "all", kc => "all")[1]
    @test r[p2_tcol(T2, "step")] == 1280 && r[p2_tcol(T2, "regular")] == 11                      # a window row: the window W = b / 2
    r = p2_rowsof(T2, "indicator" => "WBA", "k" => 3, "set" => "regular", "source" => "all", kc => "all")[1]
    @test r[p2_tcol(T2, "undecided")] == 1                                                         # orbit 15 has the flag zero at every window
end

@testset "D89p the data tables of F1 to F6 and their content" begin
    st = p2_cached("default"); res = p2_res("default", st); (; ref, S) = p2_refs("default", st); p = st.p
    # F1: A_XY, kappa_XY and the 2 x 2 counts of every pair of I_8 at b_P and N_c on W2 and S2
    F1 = res.items.F1
    @test length(F1.rows) == 2 * 2 * 28
    function f1row(set, k, X, Y)
        rs = p2_rowsof(F1, "set" => set, "budget" => p.budgets[k + 1], "X" => String(X), "Y" => String(Y))
        length(rs) == 1 ? rs[1] : nothing
    end
    badf = String[]
    for (set, ids) in (("W2", S.W2), ("S2", S.S2)), k in (6, 9), a in 1:8, b in (a + 1):8
        X, Y = P2_X8[a], P2_X8[b]
        r = f1row(set, k, X, Y)
        r === nothing && (push!(badf, "$set $k $X $Y missing"); continue)
        have = [i for i in ids if haskey(ref.V1, (i, X)) && haskey(ref.V1, (i, Y))]
        t = p2ref_table([ref.V1[(i, X)][k + 1] for i in have], [ref.V1[(i, Y)][k + 1] for i in have])
        got = (r[p2_tcol(F1, "cc")], r[p2_tcol(F1, "cr")], r[p2_tcol(F1, "rc")], r[p2_tcol(F1, "rr")], r[p2_tcol(F1, "N")])
        (got == (t.cc, t.cr, t.rc, t.rr, t.n) && p2_same(r[p2_tcol(F1, "A")], t.A) && p2_same(r[p2_tcol(F1, "kappa")], t.kappa)) || push!(badf, "$set $k $X $Y")
    end
    @test isempty(badf)
    # hand row: W2 at b_P, (LE, E1): LE chaotic on 7 orbits and regular on 4, E1 regular on all 11: cr = 7, rr = 4, A = 4 / 11, kappa = 0 (p_e = p_o)
    r = f1row("W2", 6, :LE, :E1)
    @test (r[p2_tcol(F1, "cc")], r[p2_tcol(F1, "cr")], r[p2_tcol(F1, "rc")], r[p2_tcol(F1, "rr")]) == (0, 7, 0, 4) && isapprox(r[p2_tcol(F1, "A")], 4 / 11) && abs(r[p2_tcol(F1, "kappa")]) < 1e-15
    # F2: the distribution of the censored detection time and of the outcomes per indicator for S2, W2, T_s: the counts add up to the orbits that have the item
    F2 = res.items.F2
    sets3 = Dict("S2" => S.S2, "W2" => S.W2, "T_s" => S.T_s)
    badg = String[]
    for (nm, ids) in sets3, X in P2_X8
        n = count(i -> haskey(ref.Det, (i, X)), ids)
        tt = p2_rowsof(F2, "set" => nm, "indicator" => String(X), "kind" => "t_censored")
        oc = p2_rowsof(F2, "set" => nm, "indicator" => String(X), "kind" => "outcome")
        sum(r[p2_tcol(F2, "count")] for r in tt; init = 0) == n || push!(badg, "$nm $X t")
        sum(r[p2_tcol(F2, "count")] for r in oc; init = 0) == n || push!(badg, "$nm $X outcome")
        for r in tt
            r[p2_tcol(F2, "value")] in vcat(P2_BUDGETS[1:8], [40960]) || push!(badg, "$nm $X value $(r[p2_tcol(F2, "value")])")
            r[p2_tcol(F2, "count")] == count(i -> haskey(ref.Det, (i, X)) && ref.Det[(i, X)].censored == r[p2_tcol(F2, "value")], ids) || push!(badg, "$nm $X count")
        end
        for r in oc
            r[p2_tcol(F2, "count")] == count(i -> haskey(ref.Det, (i, X)) && String(ref.Det[(i, X)].outcome) == r[p2_tcol(F2, "value")], ids) || push!(badg, "$nm $X outcome count")
        end
    end
    @test isempty(badg)
    # F3: the fraction of chaotic verdicts per indicator and budget for W2, S2, T_s and the regular set (decided orbits in the denominator)
    F3 = res.items.F3
    badh = String[]
    for (nm, ids) in (("W2", S.W2), ("S2", S.S2), ("T_s", S.T_s), ("regular", S.regular)), X in P2_X8, k in 0:9
        rs = p2_rowsof(F3, "set" => nm, "indicator" => String(X), "k" => k)
        length(rs) == 1 || (push!(badh, "$nm $X $k"); continue)
        vs = [ref.V1[(i, X)][k + 1] for i in ids if haskey(ref.V1, (i, X))]
        dec = count(in((c_, r_)), vs); ch = count(==(c_), vs)
        (rs[1][p2_tcol(F3, "decided")] == dec && p2_same(rs[1][p2_tcol(F3, "fraction")], dec == 0 ? nothing : ch / dec)) || push!(badh, "$nm $X $k value")
    end
    @test isempty(badh)
    # F4 holds the AUC rows (T9 and F4 share the data): one row per (indicator, budget, K, source)
    F4 = res.items.F4
    @test length(F4.rows) == length(res.auc) && Set(r[p2_tcol(F4, "indicator")] for r in F4.rows) == Set(String.(vcat(P2_X8, :STR)))
    r = p2_rowsof(F4, "indicator" => "LE", "k" => 6, p2_kfield(F4) => p2_klabel(F4, 0.5, 1), "source" => "pooled")[1]
    @test (r[p2_tcol(F4, "n_C")], r[p2_tcol(F4, "n_R")]) == (13, 7) && p2_same(r[p2_tcol(F4, "AUC")], p2_kget(res.auc, K2 -> (:LE, 6, K2, "pooled"), 0.5, 1).auc)
    # F5: the lattice of each K index: one row per grid orbit with rho and the verdicts of I_8 at b_P and N_c
    F5 = res.items.F5
    @test length(F5.rows) == 24 && [r[p2_tcol(F5, "id")] for r in F5.rows] == collect(21:44)
    r = F5.rows[1]
    @test r[p2_tcol(F5, "label")] == "grid-k1-i00-j00" && r[p2_tcol(F5, "rho")] == "C" && r[p2_tcol(F5, "K index")] == 1 && r[p2_tcol(F5, "i")] == 0 && r[p2_tcol(F5, "j")] == 0
    @test all(r -> all(X -> r[p2_tcol(F5, "$X b_P")] == String(ref.V1[(r[p2_tcol(F5, "id")], X)][7]) && r[p2_tcol(F5, "$X N_c")] == String(ref.V1[(r[p2_tcol(F5, "id")], X)][10]), P2_X8), F5.rows)
    # F6: the count and the 5th, 25th, 50th, 75th and 95th percentile of the stretching D at each recorded step for the classes C and R and for T_s
    F6 = res.items.F6
    cls = Dict("C" => [i for i in S.E_all if res.rho[i] === :C], "R" => S.regular, "T_s" => S.T_s)
    badi = String[]
    for (nm, ids) in cls, t in (P2_GB[1], 44, 2816, 20480)
        rs = p2_rowsof(F6, "class" => nm, "step" => t)
        length(rs) == 1 || (push!(badi, "$nm $t"); continue)
        sc = [ref.items[(i, "st", "Float64", p.N_c, 1)].indicators["st"].series[searchsortedfirst(P2_GB, t), 1] for i in ids]
        want = (length(sc), (p2ref_percentile(sc, q) for q in (0.05, 0.25, 0.5, 0.75, 0.95))...)
        got = (rs[1][p2_tcol(F6, "count")], (rs[1][p2_tcol(F6, c)] for c in ("p5", "p25", "p50", "p75", "p95"))...)
        p2_same(collect(got), collect(want)) || push!(badi, "$nm $t value")
    end
    @test isempty(badi) && length(F6.rows) == 3 * length(P2_GB)
end

@testset "D89p the tables file: header, sections, digests, readings, deviations (E37)" begin
    st = p2_cached("default"); res = p2_res("default", st)
    text = read(res.tables, String)
    h = res.header
    @test h.prereg_header_sha256 == st.prereg_hash && h.panel_sha256 == st.panel_hash && h.manifest_sha256 == an_sha256(read(st.manifest)) && h.params_sha256 == an_sha256(read(st.params))
    @test h.p1_prereg_header_sha256 == st.p1.prereg_hash && h.p1_manifest_sha256 == an_sha256(read(st.p1.manifest)) && h.p1_params_sha256 == an_sha256(read(st.p1.params))
    @test h.plan_sha256 == st.p.plan_sha256 && h.chunks == 33 && h.items == 466 && h.seed == AN_SEED
    head = first(split(text[findfirst("## Header", text).start:end], "\n## "; limit = 2))
    for hs in (st.prereg_hash, st.panel_hash, h.manifest_sha256, h.params_sha256, st.p1.prereg_hash, h.p1_manifest_sha256, h.p1_params_sha256, st.p.plan_sha256)
        @test occursin(hs, head)
    end
    @test occursin("module version", head) && occursin("git commit", head) && occursin(string(AN_SEED), head) && all(k -> occursin(string(k), head), P2_PARAM_KEYS)        # the 58 parameters
    @test !occursin(st.dir, text) && !occursin(AN_DIR, text) && !occursin(r"20\d\d-\d\d-\d\dT\d\d:\d\d", text)                    # no clock, no absolute path
    # sections: the decisions, T1 to T17 and F1 to F6 (data tables with digests), readings, deviations, exploratory items
    @test occursin("## Decisions", text)
    for id in vcat(["T$i" for i in 1:17])
        @test occursin(Regex("(?m)^## $id[^\\n]*confirmatory"), text)
    end
    for i in 1:6
        @test occursin(Regex("(?m)^## F$i[^\\n]*confirmatory"), text) || occursin(Regex("(?m)^## F$i data"), text)
        j = findfirst("## F$i data", text)
        @test j !== nothing
        m = match(r"digest: ([0-9a-f]{64})", text[j.start:end])
        @test m !== nothing && m.captures[1] == res.digests["F$i"]
        # the digest is the SHA-256 of the table written with 12 significant digits and LF line ends, recomputed here from the numbers
        t = getproperty(res.items, Symbol("F$i"))
        @test p2_digest(t.columns, t.rows) == res.digests["F$i"]
    end
    @test Set(keys(res.digests)) == Set("F$i" for i in 1:6) && Set(propertynames(res.items)) == Set(Symbol.(vcat(["T$i" for i in 1:17], ["F$i" for i in 1:6])))
    for sec in ("## Readings", "## Deviations", "## Exploratory items")
        @test occursin(sec, text)
    end
    ex = text[findfirst("## Exploratory items", text).start:end]
    @test occursin("none", lowercase(first(split(ex, "\n## "))))                                             # I8: no exploratory item
    rd = first(split(text[findfirst("## Readings", text).start:end], "\n## "))
    @test count(l -> startswith(l, "- "), split(rd, '\n')) >= 15               # every reading chosen is listed, one line each
    dv = join(res.deviations, "\n")
    low = lowercase(dv)
    @test occursin("overrides", low) && (occursin("reduced", low) || occursin("grid", low)) && occursin("lambda_ref", low) && occursin("sali", low) && occursin("undecided", low) &&
          (occursin("n(b)", low) || occursin("first recorded step", low) || occursin("budget step", low)) && all(d -> d isa String, res.deviations)
    dtext = first(split(text[findfirst("## Deviations", text).start:end], "\n## "))
    @test all(d -> occursin(d, dtext), res.deviations)
    # the decisions section prints the verdicts, the ROBUST flags and the gate result
    dec = first(split(text[findfirst("## Decisions", text).start:end], "\n## "))
    @test all(h -> occursin(h, dec), ("H1", "H2", "H3", "H4", "H5", "H6")) && occursin("ROBUST", dec) && occursin("NOT COMPARABLE", dec) && !occursin("rule invalid on the integrable control", dec)
    @test occursin("not decidable", lowercase(dec)) && occursin("supported", lowercase(dec))
end

@testset "D89r read-only: inputs, store and data tree unchanged; D89q determinism at 1 and 4 threads and another folder" begin
    st = p2_copy(p2_cached("default"), "det")
    before = p2_fingerprint(st); tree = p2_datatree_hash()
    a = p2_run(st; threads = 1, figures = true)
    ta = read(a.tables)
    @test a.status === :figures_skipped && isempty(filter(f -> endswith(f, ".pdf"), an_files(st.out)))        # D89u: no drawer for the study, tables only
    @test p2_fingerprint(st) == before && p2_datatree_hash() == tree                                    # I12: nothing but --out changes, nothing under data/
    @test all(f -> startswith(f, "p2-indicators"), an_files(st.out)) && an_files(st.out) == ["p2-indicators-tables.md"]
    b = p2_run(st; threads = 4, out = joinpath(st.dir, "out-other"))
    P2_DETRUNS[:a] = a; P2_DETRUNS[:b] = b                                    # I13: any thread count, any output folder
    @test read(b.tables) == ta && a.digests == b.digests
    # one command regenerates everything: delete the output folder and run again
    rm(st.out; recursive = true)
    d = p2_run(st; threads = 3)
    P2_DETRUNS[:d] = d
    @test read(d.tables) == ta
    @test p2_datatree_hash() == tree
    # the in-memory P1 labelling writes no P1 output: the P1 folder holds the inputs and the store only
    @test sort(readdir(st.p1.dir)) == sort(["p1-le-benchmark-manifest.json", "p1-le-benchmark-panel.toml", "p1-le-benchmark-parameters.toml", "p1-le-benchmark-preregistration.md", "store"])
end

@testset "D89s command: usage errors, exit code 2 with the figures skipped, one command regenerates" begin
    st = p2_copy(p2_cached("default"), "cli")
    io = IOBuffer()
    @test AN.cli(p2_cli_args(st; extra = ["--threads", "2"]); io = io) == 2
    msg = String(take!(io))
    @test occursin("figures skipped", msg) && isfile(joinpath(st.out, "p2-indicators-tables.md"))
    ta = isfile(joinpath(st.out, "p2-indicators-tables.md")) ? read(joinpath(st.out, "p2-indicators-tables.md")) : UInt8[]
    ispath(st.out) && rm(st.out; recursive = true)
    io = IOBuffer(); @test AN.cli(p2_cli_args(st); io = io) == 2
    @test isfile(joinpath(st.out, "p2-indicators-tables.md")) && !isempty(ta) && read(joinpath(st.out, "p2-indicators-tables.md")) == ta
    # the options --p1-* are accepted for p2-indicators; usage errors write nothing
    a = p2_cli_args(st)
    for (tag, args) in (("value", vcat(a[1:5], ["--store"])), ("unknown", vcat(a[1:5], ["--bogus", "x", "--out", an_path("p2-cli-o1")])), ("few", a[1:3]), ("none", ["analyze"]),
                        ("p1value", vcat(a[1:5], ["--out", an_path("p2-cli-o3"), "--p1-store"])), ("threads", vcat(a[1:5], ["--out", an_path("p2-cli-o4"), "--threads", "many"])))
        io = IOBuffer()
        @test AN.cli(args; io = io) == 1
        m = String(take!(io))
        @test (occursin("analyze", m) || occursin("usage", lowercase(m))) && !ispath(an_path("p2-cli-o$tag"))
    end
    io = IOBuffer(); @test AN.cli(["analyze", "--help"]; io = io) == 0
    h = String(take!(io))
    @test all(o -> occursin(o, h), ("--p1-store", "--p1-preregistration", "--p1-manifest", "--p1-parameters", "--store", "--out", "--threads"))
    # a refusal: exit code 1, the cause named, nothing written
    bad = p2_copy(p2_cached("default"), "cli-bad"); write(bad.prereg, read(bad.prereg, String) * "x")
    io = IOBuffer(); @test AN.cli(p2_cli_args(bad); io = io) == 1
    @test occursin("preregistration", lowercase(String(take!(io)))) && !ispath(bad.out)
end

@testset "D89u figure data: six digests, equal in every run, and the figure step without a drawer" begin
    st = p2_cached("default"); res = p2_res("default", st)
    @test all(occursin(r"^[0-9a-f]{64}$", d) for d in values(res.digests)) && Set(keys(res.digests)) == Set("F$i" for i in 1:6)
    @test res.digests == P2_DETRUNS[:a].digests == P2_DETRUNS[:b].digests == P2_DETRUNS[:d].digests              # the runs of the determinism testset
end
