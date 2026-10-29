# Review findings of stage 04 for Amendment D89 (tests written from the spec, not from the code). One testset per finding; the finding numbers are those of the review
# (8, 10 and 11 became readings and are not tested). Spec: analysis-spec.md, Amendment D89, E22, E24, E31, E36, T16.
#
# The study "review" is the default study with the grid K list swapped: grid K values ["0.75", "0.5"] (K index 1 = 0.75, K index 2 = 0.5), so that a panel orbit and a grid orbit
# with the same K VALUE have different K indexes (panel index 1 / grid index 2 at K = 0.5) and orbits with the same index have different K values (panel and grid index 1: 0.5 and
# 0.75; index 2: 2 pi and 0.5). Design variant :review adds: FLI of the regular orbits 41 and 42 undecided at k = 0 .. 2 and FLI of the sticky orbit 39 chaotic from k = 3 (finding 5),
# nu1 = NaN of the gate orbit 10 at window 320 (finding 6).

const P2_REVIEW_K = [0.5, 0.75, 2 * pi]
p2_review_study() = p2_cached("review"; variant = :review, p = p2_params(grid_k_values = ["0.75", "0.5"]))
# a K label of a table or dictionary key: a number equal to one of the K values (a K index 1, 2 is not)
p2_isK(x) = x isa Real && any(K -> isapprox(Float64(x), K; atol = 1e-12), P2_REVIEW_K)
function p2_Kof(x)
    i = x isa Real ? findfirst(K -> isapprox(Float64(x), K; atol = 1e-12), P2_REVIEW_K) : nothing
    return i === nothing ? NaN : P2_REVIEW_K[i]                    # NaN for a label that is not a K value (a K index, "all")
end
# the entry of a dictionary whose key satisfies pred, or `default`
function p2_getK(d, pred, default)
    for (k, v) in d
        pred(k) && return v
    end
    return default
end
# true if f() throws an ArgumentError whose message contains all the words of one alternative
function p2_rejects_any(f, alternatives)
    try
        f()
        return false
    catch e
        e isa ArgumentError || return false                      # a KeyError, BoundsError ... is not a named refusal
        msg = lowercase(sprint(showerror, e))
        return any(ws -> all(w -> occursin(lowercase(w), msg), ws), alternatives)
    end
end

@testset "review 1: pooling by K value, not by K index (AUC, regular range, leaving budget, share below the median, T9, F4, T14)" begin
    st = p2_review_study(); res = p2_res("review", st); (; ref, S) = p2_refs("review", st); p = st.p
    orbs = p2_ref_orbits(st)
    @test p2_ref_Kvals(orbs, S) == P2_REVIEW_K && [g.K for g in st.h4[[1, 13]]] == ["0.75", "0.5"]
    # the study has the intended collisions: panel index 1 = 0.5 and grid index 2 = 0.5 (same value, other index); panel index 2 = 2 pi and grid index 2 = 0.5 (same index, other value)
    byid = Dict(o.id => o for o in orbs)
    @test byid[1].kindex == 1 && byid[1].K == 0.5 && byid[33].kindex == 2 && byid[33].K == 0.5 && byid[21].kindex == 1 && byid[21].K == 0.75 && byid[11].kindex == 2 && byid[11].K == 2 * pi
    # AUC: rows labelled by the K value; counts and values equal to the reference grouped by the K value; hand counts of the plan
    auc = p2_ref_auc_K(st, ref, S)
    @test Set(Float64(k[3]) for k in keys(res.auc) if k[3] isa Real) == Set(P2_REVIEW_K) && all(k -> p2_isK(k[3]), keys(res.auc))
    @test Set((X, k, p2_Kof(K), src) for (X, k, K, src) in keys(res.auc)) == Set(keys(auc))
    A(X, k, K, src) = p2_getK(res.auc, key -> key[1] == X && key[2] == k && key[4] == src && p2_Kof(key[3]) == K, (n_C = -1, n_R = -1, auc = :missing))
    @test isempty([key for (key, want) in auc if !p2_same(A(key...), want)])
    # K = 0.5: chaotic orbits of E_all 7 (panel 1, 3 to 8) + 8 (grid 33 to 40), regular 2 (9, 10) + 3 (41 to 43); K = 0.75: grid only, 6 and 5; K = 2 pi: panel only, 3 and 2
    @test (a = A(:LE, 6, 0.5, "pooled"); (a.n_C, a.n_R) == (15, 5) && a.auc !== nothing) && (a = A(:LE, 6, 0.75, "pooled"); (a.n_C, a.n_R) == (6, 5) && a.auc !== nothing)
    @test (a = A(:LE, 6, 2 * pi, "pooled"); (a.n_C, a.n_R) == (3, 2) && a.auc === nothing)
    @test (a = A(:LE, 6, 0.5, "panel"); (a.n_C, a.n_R) == (7, 2) && a.auc === nothing) && (a = A(:LE, 6, 0.5, "grid"); (a.n_C, a.n_R) == (8, 3) && a.auc === nothing)
    @test (a = A(:LE, 6, 0.75, "grid"); (a.n_C, a.n_R) == (6, 5) && a.auc !== nothing) && (a = A(:LE, 6, 0.75, "panel"); (a.n_C, a.n_R) == (0, 0) && a.auc === nothing)
    # regular range at the K value (at least min_set = 5 regular orbits): K = 0.5 holds {9, 10, 41, 42, 43} and K = 0.75 holds {27 .. 31}; K = 2 pi holds two: no range
    rng = p2_ref_ranges_K(st, ref, S)
    R(X, K, k) = p2_getK(res.range, key -> key[1] == X && key[3] == k && p2_Kof(key[2]) == K, :missing)
    @test Set((X, p2_Kof(K), k) for (X, K, k) in keys(res.range)) == Set(keys(rng)) && all(key -> p2_isK(key[2]), keys(res.range))
    @test isempty([key for (key, want) in rng if !(let got = R(key...); want === nothing ? got === nothing : (got !== nothing && got !== :missing && p2_same(got, want)) end)])
    nR(X, K, k) = (r = R(X, K, k); r isa NamedTuple ? r.n : r)
    @test nR(:LE, 0.5, 3) == 5 && nR(:LE, 0.75, 3) == 5 && nR(:LE, 2 * pi, 3) === nothing
    # leaving budget of the sticky orbits (5, 39 at K = 0.5; 23, 24 at K = 0.75) against the range of their own K; the share below the regular median at k = 9
    lv = p2_ref_leaving_K(st, ref, S, rng)
    @test isempty([key for (key, want) in lv if get(res.leaving, key, :missing) !== want])
    bl = p2_ref_below_K(st, ref, S, rng)
    Bm(X, K) = p2_getK(res.below_median, key -> key[1] == X && p2_Kof(key[2]) == K, :missing)
    @test Set((X, p2_Kof(K)) for (X, K) in keys(res.below_median)) == Set(keys(bl)) && isempty([key for (key, want) in bl if !p2_same(Bm(key...), want)])
    # the tables carry the same rows, labelled by the K value; the Y5 control (K absent) is in no row
    for T in (res.items.T9, res.items.F4)
        kc = p2_kfield(T)
        labels = [r[p2_tcol(T, kc)] for r in T.rows]
        @test all(p2_isK, labels) && Set(p2_Kof(x) for x in labels) == Set(P2_REVIEW_K)
        r = p2_rowsof(T, "indicator" => "LE", "k" => 6, kc => 0.5, "source" => "pooled")
        @test length(r) == 1 && (r[1][p2_tcol(T, "n_C")], r[1][p2_tcol(T, "n_R")]) == (15, 5)
    end
    T14 = res.items.T14
    labels = [r[p2_tcol(T14, p2_kfield(T14))] for r in T14.rows]
    @test all(p2_isK, labels) && !any(x -> x isa Real && isnan(x), labels)
    r = [r for r in T14.rows if r[p2_tcol(T14, "part")] == "range" && r[p2_tcol(T14, "indicator")] == "LE" && r[p2_tcol(T14, "k")] == 3 && r[p2_tcol(T14, p2_kfield(T14))] == 0.5]
    @test length(r) == 1 && r[1][p2_tcol(T14, "n")] == 5
end

@testset "review 1: the by-K rows of T2, T3 to T8, T10 and T13 pool by K value and are labelled by it; the pooled decisions are unchanged" begin
    st = p2_review_study(); res = p2_res("review", st); (; ref, S) = p2_refs("review", st)
    orbs = p2_ref_orbits(st); byid = Dict(o.id => o for o in orbs)
    # T2: the counts of the by-K rows are those of the orbits of the set with that K value (source all), whatever their panel or grid index
    T2 = res.items.T2; kc = p2_kfield(T2)
    sets = Dict("W2" => S.W2, "S2" => S.S2, "T_s" => S.T_s, "regular" => S.regular)
    cols = ["chaotic", "regular", "undecided", "not applicable"]
    bad = String[]; seen = Set{Float64}()
    for r in T2.rows
        lab = r[p2_tcol(T2, kc)]
        (lab == "all" || r[p2_tcol(T2, "source")] != "all") && continue
        p2_isK(lab) || (push!(bad, "label $lab"); continue)
        push!(seen, p2_Kof(lab))
        X, k, nm = Symbol(r[p2_tcol(T2, "indicator")]), r[p2_tcol(T2, "k")], r[p2_tcol(T2, "set")]
        ids = [id for id in sets[nm] if Float64(byid[id].K) == p2_Kof(lab)]
        vs = [get(ref.V1, (id, X), nothing) for id in ids]
        want = (count(v -> v !== nothing && v[k + 1] === c_, vs), count(v -> v !== nothing && v[k + 1] === r_, vs), count(v -> v !== nothing && v[k + 1] === u_, vs), count(isnothing, vs))
        Tuple(r[p2_tcol(T2, c)] for c in cols) == want || push!(bad, "$X $k $nm $lab")
    end
    @test isempty(bad) && seen == Set(P2_REVIEW_K)
    # T3 to T8: the statistic "orbits" of the by-K rows (hand counts: W2 = {1, 3, 5, 7, 8, 39} at 0.5 and {21, 22, 23, 24, 26} at 0.75; S2 = {33 .. 38} at 0.5 and {11, 12} at 2 pi;
    # T_s = {5, 39} at 0.5 and {23, 24} at 0.75)
    expect = Dict(:T3 => ("orbits", Dict(0.5 => 6, 0.75 => 5)), :T4 => ("orbits", Dict(0.5 => 6, 2 * pi => 2)), :T5 => ("orbits W2", Dict(0.5 => 6, 0.75 => 5)),
                  :T6 => ("orbits", Dict(0.5 => 2, 0.75 => 2)), :T7 => ("orbits", Dict(0.5 => 6, 0.75 => 5)))
    for (nm, (stat, counts)) in expect
        T = getproperty(res.items, nm); kc = p2_kfield(T)
        rows = [r for r in T.rows if r[p2_tcol(T, "source")] == "all" && r[p2_tcol(T, kc)] != "all" && r[p2_tcol(T, "statistic")] == stat]
        @test all(r -> p2_isK(r[p2_tcol(T, kc)]), rows) && !isempty(rows)
        got = Dict(p2_Kof(r[p2_tcol(T, kc)]) => r[p2_tcol(T, "value")] for r in rows if p2_isK(r[p2_tcol(T, kc)]))
        @test all(((K, n),) -> get(got, K, 0) == n, counts) && all(((K, n),) -> n == 0 || haskey(counts, K), got)
    end
    for nm in (:T8, :T10, :T13)                              # tables whose by-K rows (if any) carry K values and no other label
        T = getproperty(res.items, nm); kc = p2_kfield(T)
        @test all(r -> r[p2_tcol(T, kc)] == "all" || p2_isK(r[p2_tcol(T, kc)]), T.rows)
    end
    # the pooled decisions do not depend on the K labelling
    rf = p2_ref_hyp(st, ref, S, P2_X7)
    @test p2_decision_diff(res.decisions, rf) == String[] && p2_decision_diff(res.r2, p2_ref_hyp(st, ref, S, P2_X7; rule = :R2)) == String[]
    @test p2_decision_diff(res.extension, p2_ref_hyp(st, ref, S, P2_X8; suffix = "_8", srule = "E8")) == String[]
    @test res.sets.W2 == S.W2 && res.sets.S2 == S.S2 && res.sets.T_s == S.T_s
end

@testset "review 2: type 7 percentile with infinite neighbours is never NaN (E31)" begin
    p = p2_params()
    # 100 scores, ten of them -Inf (RLI = 0 gives log10 = -Inf): the 5th percentile has h = 99 x 0.05 + 1 = 5.95, both neighbours (the 5th and 6th values) are -Inf: the value is -Inf
    s1 = vcat(fill(-Inf, 10), collect(1.0:90.0))
    rr = AN.p2_regular_range(s1, p)
    @test rr.lo == -Inf && !isnan(rr.lo)
    @test isapprox(rr.hi, 85.05; atol = 1e-12) && rr.median == 40.5 && rr.n == 100                # h = 95.05: the 95th and 96th values are 85 and 86
    # one -Inf next to a finite value: h = 4 x 0.05 + 1 = 1.2 lies between -Inf and 1: no NaN; the value is the infinite one or the order statistic nearer to h (both are -Inf)
    s2 = AN.p2_regular_range([-Inf, 1.0, 2.0, 3.0, 4.0], p)
    @test !isnan(s2.lo) && s2.lo == -Inf && !isnan(s2.hi) && isapprox(s2.hi, 3.8; atol = 1e-12) && s2.median == 2.0       # h = 4 x 0.95 + 1 = 4.8 between 3 and 4
    # an interpolation weight nearer to the finite neighbour: percentile 15 gives h = 4 x 0.15 + 1 = 1.6 between -Inf and 1: the finite order statistic or the infinite one, never NaN
    s3 = AN.p2_regular_range([-Inf, 1.0, 2.0, 3.0, 4.0], merge(p, (; range_lower_percentile = 15.0)))
    @test !isnan(s3.lo) && s3.lo in (-Inf, 1.0)
    # all finite: unchanged (type 7: 1.95 on 1 .. 20)
    @test isapprox(AN.p2_regular_range(collect(1.0:20.0), p).lo, 1.95; atol = 1e-12)
    # the median with infinities
    @test AN.p2_regular_range([-Inf, -Inf, -Inf, 1.0, 2.0, 3.0], p).median == -Inf && !isnan(AN.p2_regular_range([-Inf, -Inf, -Inf, 1.0, 2.0, 3.0], p).lo)
end

@testset "review 3: a chunk with a stripped attribute is refused in P2, naming it (E36 b, d, e, f)" begin
    st = p2_cached("default")
    failed = String[]
    for name in ("plan_sha256", "prereg_sha256", "panel_sha256", "status", "study")
        c = p2_copy(st, "rev3-$name")
        p2_del_attr(p2_chunk(c, 5), "", name)
        ok = an_rejects(() -> p2_run(c), [name]) && !ispath(c.out)
        ok || push!(failed, name)
        GC.gc(); try rm(c.dir; recursive = true, force = true) catch end
    end
    @test isempty(failed)
    # the P1 behaviour stays: a P1 chunk without panel_sha256 or prereg_sha256 is accepted (the P1 helper does not require them)
    q = p1_copy(p1_gold(), "rev3-p1")
    f = joinpath(q.store, "chunk-00005.h5")
    p2_del_attr(f, "", "panel_sha256"); p2_del_attr(f, "", "prereg_sha256")
    r = AN.analyze(; p1_args(q)..., out = joinpath(q.dir, "o"))
    @test r.status === :ok
end

@testset "review 4: a numeric text must equal its value (E22, I14, A2)" begin
    st = p2_cached("default")
    base = merge(p2_params(), (; panel_sha256 = st.p.panel_sha256))
    # the text and the quote agree with each other (the quote follows the text) but not with the value: refused naming the key
    cases = (("persist", :persist_divisor, base, "8"), ("minset", :min_set, base, "50"), ("e1", :e1_d_threshold, base, "-5"), ("q", :h2_q_min, merge(base, (; h2_q_min = 100.0)), "10"),
             ("budgets", :budgets, base, "40, 80, 160, 320, 640, 1280, 2560, 5120, 10240, 20481"), ("level", :bootstrap_level, base, "0.5"), ("resamples", :bootstrap_resamples, base, "1e4"))
    failed = String[]
    for (tag, key, p, text) in cases
        pre, par = p2_pair("rev4-$tag", p; texts = Dict(key => text))
        ok = an_rejects(() -> AN.p2_read_params(par; prereg = pre, p1_params = st.p1.params), [string(key)])
        ok || push!(failed, tag)
    end
    @test isempty(failed)
    # numeric texts that equal the value are accepted whatever their format; texts that are not numbers are not compared (95 percent, 5th, 10^3 are in the default table)
    for (tag, key, text) in (("fmt1", :h1_u_min, "0.2000"), ("fmt2", :bootstrap_resamples, "1e3"), ("fmt3", :le_floor, "0.0001"), ("fmt4", :e1_d_threshold, "-4.0"), ("fmt5", :min_set, "5"))
        pre, par = p2_pair("rev4-$tag", base; texts = Dict(key => text))
        got = AN.p2_read_params(par; prereg = pre, p1_params = st.p1.params)
        @test isequal(getfield(got, key), getfield(base, key))
    end
    tab = TOML.parsefile(st.params)
    @test tab["bootstrap_level"]["text"] == "95 percent" && tab["range_lower_percentile"]["text"] == "5th" && tab["bootstrap_resamples"]["text"] == "10^3"
    @test AN.p2_read_params(st.params; prereg = st.prereg, p1_params = st.p1.params).bootstrap_level == 0.95
end

@testset "review 5: the leaving budget when the regular range does not exist at the first budgets (E31, P25)" begin
    p = p2_params()
    sc = Union{Nothing, Float64}[9, 9, 9, 9, 9, 9, 9, 9, 9, 9]
    # no range at k = 0 .. 2 (nothing): a budget without a range counts as not above, the result is the budget of k = 3, b_3 = 320, not undefined
    up = Union{Nothing, Float64}[nothing, nothing, nothing, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0]
    @test AN.p2_leaving_budget(sc, up, p) == 320
    @test AN.p2_leaving_budget(sc, Union{Nothing, Float64}[3.0, 3.0, 3.0, 3.0, 3.0, nothing, 3.0, 3.0, 3.0, 3.0], p) == 2560           # no range at k = 5: from k = 6
    @test AN.p2_leaving_budget(sc, Union{Nothing, Float64}[3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, nothing], p) === nothing          # no range at k = 9: never above
    @test AN.p2_leaving_budget(sc, fill(nothing, 10), p) === nothing
    # end to end (design variant :review on the default grids, so that the K index and the K value agree and the finding is not mixed up with the pooling of finding 1): the regular set at
    # K = 2 pi is {15, 16, 41, 42, 43}, five orbits = min_set; FLI is undecided for 41 and 42 at k = 0 .. 2, so the range of FLI exists from k = 3; the sticky orbit 39 (K = 2 pi) has
    # FLI chaotic from k = 3 (score r >= 8 against the regular r of 1.2 to 1.5): its leaving budget is b_3 = 320
    st = p2_cached("review-dg"; variant = :review); res = p2_res("review-dg", st); (; ref, S) = p2_refs("review-dg", st)
    rng = p2_ref_ranges_K(st, ref, S)
    @test all(k -> rng[(:FLI, 2 * pi, k)] === nothing, 0:2) && all(k -> rng[(:FLI, 2 * pi, k)] !== nothing && rng[(:FLI, 2 * pi, k)].n == 5, 3:9)
    @test p2_str(ref.V1[(39, :FLI)]) == "rrrccccccc" && S.T_s == [5, 23, 24, 39]
    @test res.leaving[(39, :FLI)] == 320
    @test p2_ref_leaving_K(st, ref, S, rng)[(39, :FLI)] == 320
end

@testset "review 6: non-finite nu1 distances of the K = 0 line are filtered and counted (T16)" begin
    st = p2_review_study(); res = p2_res("review", st); (; ref, S) = p2_refs("review", st); p = st.p
    kl = p2_ref_k0(st, ref, S)
    # orbit 10 (gate) has nu1 = NaN at window 320 (E1 and E2 undecided there, the flags degenerate and undefined window allow it): the distance of that orbit is left out and counted,
    # the largest distance is that of orbit 9: k x 1e-10 with k = 4 for fq1 and k x 3e-10 for pa1
    w = 320
    @test !isnan(res.k0line[w].fq1) && !isnan(res.k0line[w].pa1) && isapprox(res.k0line[w].fq1, 4e-10; atol = 1e-14) && isapprox(res.k0line[w].pa1, 1.2e-9; atol = 1e-14)
    @test res.k0line[w].fq1_nonfinite == 1 && res.k0line[w].pa1_nonfinite == 1
    @test all(W -> W == w || (res.k0line[W].fq1_nonfinite == 0 && res.k0line[W].pa1_nonfinite == 0), p.windows)
    @test all(W -> p2_same((fq1 = res.k0line[W].fq1, pa1 = res.k0line[W].pa1), (fq1 = kl[W].fq1, pa1 = kl[W].pa1); rtol = 1e-9), p.windows)
    @test all(W -> isfinite(res.k0line[W].fq1) && isfinite(res.k0line[W].pa1), p.windows)
end

@testset "review 7: an unknown indicator kind and wrong bit names are named refusals (E36)" begin
    st = p2_cached("default")
    chunk(c, i) = p2_chunk(c, i)
    failed = String[]
    function case(tag, alternatives, mutate)
        c = p2_copy(st, "rev7-$tag")
        mutate(c)
        ok = p2_rejects_any(() -> p2_run(c), alternatives) && !ispath(c.out)
        ok || push!(failed, tag)
        GC.gc(); try rm(c.dir; recursive = true, force = true) catch end
    end
    ids(c, i) = [o.id for o in AX.read_chunk(chunk(c, i)).orbits]
    rename(bits, j, nm) = (b = copy(bits); b[j] = nm; b)
    # an unknown kind in the root attributes of a chunk of an indicator spec: an ArgumentError naming the kind (not a KeyError)
    case("kind", [["kind", "bogus"]], c -> p2_set_attr(chunk(c, 10), "", "kind", "bogus"))
    # the bit names are verified at every position that the reader uses: the second name (degenerate / undefined window / zero), the third of the frequency table (failed) and of the
    # Birkhoff table (reldig_undefined)
    case("bit2", [["degenerate"], ["bogus_bit"]], c -> foreach(id -> p2_set_attr(chunk(c, 10), "orbits/$id/indicators/fq1/flags", "bits", rename(AX.IND_FREQ_BITS, 2, "bogus_bit")), ids(c, 10)))
    case("bit-failed", [["failed"], ["bogus_bit"]], c -> foreach(id -> p2_set_attr(chunk(c, 10), "orbits/$id/indicators/fq1/flags", "bits", rename(AX.IND_FREQ_BITS, 3, "bogus_bit")), ids(c, 10)))
    # the flag is set in this chunk: the count of the manifest is renamed too, so that the recomputation of the counts (a check that is not the one under test) does not refuse the chunk
    case("bit-reldig", [["reldig_undefined"], ["bogus_bit"]], c -> (foreach(id -> p2_set_attr(chunk(c, 11), "orbits/$id/indicators/bk1/flags", "bits", rename(AX.IND_BK_BITS, 3, "bogus_bit")), ids(c, 11));
                                                                     p2_edit_manifest(c, m -> (fl = m["chunks"][11]["indicator_flags"]; fl["bogus_bit"] = pop!(fl, "reldig_undefined")))))
    case("bit-zero", [["zero"], ["bogus_bit"]], c -> foreach(id -> p2_set_attr(chunk(c, 11), "orbits/$id/indicators/bk1/flags", "bits", rename(AX.IND_BK_BITS, 2, "bogus_bit")), ids(c, 11)))
    case("bit-undefined-window", [["undefined window"], ["bogus_bit"]], c -> foreach(id -> p2_set_attr(chunk(c, 14), "orbits/$id/indicators/pa1/flags", "bits", rename(AX.IND_PA_BITS, 2, "bogus_bit")), ids(c, 14)))
    @test isempty(failed)
    isempty(failed) || println("review 7 cases not refused as expected: ", failed)
end

@testset "review 9: a control item without control_deviation is refused (E36 n)" begin
    st = p2_cached("default")
    failed = String[]
    for (tag, id) in (("17", 17), ("20", 20))                       # controls with a positive exponent (bound scope)
        c = p2_copy(st, "rev9-$tag")
        p2_del_attr(p2_chunk(c, 1), "orbits/$id", "control_deviation")
        ok = an_rejects(() -> p2_run(c), ["control_deviation"]) && !ispath(c.out)
        ok || push!(failed, tag)
        GC.gc(); try rm(c.dir; recursive = true, force = true) catch end
    end
    @test isempty(failed)
end

@testset "review 12: a growth ratio that overflows to Inf gives the verdict chaotic and no score (E26, E31)" begin
    p = p2_params()
    r = AN.p2_ratio(1e308, 1e-10)                                  # both growth values are finite and the denominator is positive; the quotient overflows
    @test r == Inf
    for (X, v) in ((:LE, (n = 100, lambda = 0.3, r = r)), (:SALI, (r = r, sali = 0.1)), (:FLI, (n = 100, r = r, fli = 3.0)), (:MEGNO, (r = r, ybar = 2.0)))
        @test AN.p2_verdict(X, :R1, v, p) === :chaotic
        @test AN.p2_score(X, v) === nothing                       # undefined, not Inf (only RLI = 0 keeps the infinity, as minus infinity)
    end
    r = AN.p2_ratio(-1e308, 1e-10)                                 # the other sign: regular, still no score
    @test r == -Inf && AN.p2_verdict(:LE, :R1, (n = 100, lambda = 0.3, r = r), p) === :regular && AN.p2_score(:LE, (n = 100, lambda = 0.3, r = r)) === nothing
    @test AN.p2_score(:RLI, (rli = 0.0,)) == -Inf && AN.p2_score(:LE, (n = 100, lambda = 0.3, r = 9.5)) == 9.5
end
