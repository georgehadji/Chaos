# Amendment D94, panel level: D94c (points in the panel), D94d (class bits), D94e (draw: hand tables and the builder), D94f (prefix and
# independence), D94g (gates and ids), D94h (draw independence of D), D94j (controls, counts, no screening value), D94p (determinism),
# D94t (a second fixture). Spec: experiments-spec.md, "Amendment D94": E32 to E38, R22 to R24, I28 to I36. API assumed: the test plan.
# Fixtures: 64 random points per K (48 in fixture B), r_max 3 (2), sizes at most 5, the K list of the base fixture (and another).

# ------------------------------------------------------------------------------------------------------------------ D94d
@testset "D94d class bits equal the per-point evaluation with the verified screening (tolerance 0)" begin
    for b in (d94_build("A", D94_FA), d94_build("B", D94_FB))
        fx = b.fx
        for i in eachindex(fx.K)
            k = b.pools["k"][string(i)]
            bits = d94_side_bits(b.pools, i)
            @test k["points"] == fx.M && length(bits) == fx.M
            @test bits == d94_ref_bits(b, i)
            @test k["bits_sha256"] == ex_sha256(bits)
            for c in ("T", "C", "S")                                           # run-length lists: maximal runs of consecutive ids, ascending
                runs = k[c]
                @test all(r -> r[1] <= r[2], runs) && all(j -> runs[j][1] > runs[j - 1][2] + 1, 2:length(runs))
            end
            @test sort(collect(keys(k))) == ["C", "S", "T", "bits_sha256", "points"]                          # I34: ids, bits and counts only
        end
    end
    # K = 0 (index 1 of the fixtures): the closed form of O2 holds for every random point, and the same rule screens K = 0 and K > 0
    b = d94_build("A", D94_FA)
    blk = ex_block_default()
    N, al = blk["N_scr"], ex_f64(blk["alpha"])
    closed = log(sqrt((cos(al) + N * sin(al))^2 + sin(al)^2)) / N
    @test all(p -> abs(EX.screen_point(0.0, p[1], p[2]; N_scr = N, w = blk["w"], alpha = al, tau = 1).lambda - closed) <= 8 * EX_EPS, d94_points(D94_FA.seed, 1, D94_FA.M))
    bits0 = d94_side_bits(b.pools, 1)
    @test count(x -> x & 0x01 != 0, bits0) > 0 && all(x -> x & 0x06 == 0, bits0)         # K = 0: regular, T only (C and S need a positive exponent)
    # the lattice stays the default and a lattice build with the key spelled out is the base panel
    lat = ex_build("lat-explicit", Dict{String, Any}("screen_source" => "lattice"))
    base = ex_build("base")
    @test lat.panel["orbit"] == base.panel["orbit"]
    @test lat.panel["manifest"]["cells"] == base.panel["manifest"]["cells"] && lat.panel["manifest"]["periodic"] == base.panel["manifest"]["periodic"]
    @test !isfile(lat.path * ".pools.toml")
end

# ------------------------------------------------------------------------------------------------------------------ D94e
d94_hb(M; kw...) = d94_bits(M; kw...)
const D94_HAND = [
    # (name, bits, quota, r_max, active, hand table (round, class) => (pool, drawn, status) with pool nothing = not checked, forced sets)
    (name = "pool smaller than the quota", bits = d94_hb(8; T = [2, 5, 7]), q = Dict("T" => 5), r_max = 1, active = nothing,
     tab = Dict((1, "T") => (3, 3, "short")), sets = Dict((1, "T") => [2, 5, 7])),
    (name = "pool runs out in round 2", bits = d94_hb(10; T = 1:7), q = Dict("T" => 4), r_max = 3, active = nothing,
     tab = Dict((1, "T") => (7, 4, "populated"), (2, "T") => (3, 3, "short"), (3, "T") => (0, 0, "empty")), sets = Dict()),
    (name = "a point in both the S and the C pool", bits = d94_hb(6; S = [3], C = [3, 5]), q = Dict("S" => 1, "C" => 2), r_max = 2, active = nothing,
     tab = Dict((1, "S") => (1, 1, "populated"), (1, "C") => (1, 1, "short"), (2, "S") => (0, 0, "empty"), (2, "C") => (0, 0, "empty")),
     sets = Dict((1, "S") => [3], (1, "C") => [5])),
    (name = "S and C pools identical", bits = d94_hb(9; S = [4, 6, 9], C = [4, 6, 9]), q = Dict("S" => 3, "C" => 2), r_max = 1, active = nothing,
     tab = Dict((1, "S") => (3, 3, "populated"), (1, "C") => (0, 0, "empty")), sets = Dict((1, "S") => [4, 6, 9])),
    (name = "a point with the C and T bits goes to C first", bits = d94_hb(6; C = [1, 2], T = [1, 2, 3]), q = Dict("C" => 2, "T" => 3), r_max = 1, active = nothing,
     tab = Dict((1, "C") => (2, 2, "populated"), (1, "T") => (1, 1, "short")), sets = Dict((1, "C") => [1, 2], (1, "T") => [3])),
    (name = "K = 0 gate: round 1 only", bits = d94_hb(10; T = 1:10), q = Dict("T" => 3), r_max = 3, active = [1],
     tab = Dict((1, "T") => (10, 3, "populated"), (2, "T") => (nothing, 0, "not drawn"), (3, "T") => (nothing, 0, "not drawn")), sets = Dict()),
    (name = "a size 0 class", bits = d94_hb(6; C = [1, 2, 3], T = [4, 5]), q = Dict("C" => 0, "T" => 2), r_max = 2, active = nothing,
     tab = Dict((1, "C") => (nothing, 0, "not drawn"), (2, "C") => (nothing, 0, "not drawn"), (1, "T") => (2, 2, "populated"), (2, "T") => (0, 0, "empty")),
     sets = Dict((1, "T") => [4, 5])),
    (name = "r_max = 1", bits = d94_hb(6; T = 1:5), q = Dict("T" => 2), r_max = 1, active = nothing,
     tab = Dict((1, "T") => (5, 2, "populated")), sets = Dict()),
    (name = "pool equal to the quota", bits = d94_hb(4; T = [1, 2, 3]), q = Dict("T" => 3), r_max = 2, active = nothing,
     tab = Dict((1, "T") => (3, 3, "populated"), (2, "T") => (0, 0, "empty")), sets = Dict((1, "T") => [1, 2, 3])),
    (name = "empty from the start", bits = d94_hb(5), q = Dict("T" => 2), r_max = 2, active = nothing,
     tab = Dict((1, "T") => (0, 0, "empty"), (2, "T") => (0, 0, "empty")), sets = Dict()),
    (name = "S exhausts over two rounds while C shares three points", bits = d94_hb(8; S = [2, 3, 4, 5], C = [3, 4, 5, 6, 7, 8]),
     q = Dict("S" => 2, "C" => 3), r_max = 2, active = nothing, tab = Dict((1, "S") => (4, 2, "populated")), sets = Dict()),
    (name = "seeded pools of 60 points", bits = (r = StableRNG(d94_seed('e', 40)); d94_hb(60; T = [j for j in 1:60 if rand(r) < 0.25], C = [j for j in 1:60 if rand(r) < 0.3],
                                                                                   S = [j for j in 1:60 if rand(r) < 0.15])),
     q = Dict("T" => 5, "C" => 4, "S" => 3), r_max = 3, active = nothing, tab = Dict(), sets = Dict()),
    (name = "a short cell stays short or empty", bits = d94_hb(5; T = [1, 2]), q = Dict("T" => 5), r_max = 3, active = nothing,
     tab = Dict((1, "T") => (2, 2, "short"), (2, "T") => (0, 0, "empty"), (3, "T") => (0, 0, "empty")), sets = Dict((1, "T") => [1, 2])),
    (name = "C before T on shared points", bits = d94_hb(5; C = [1, 2], T = [1, 2]), q = Dict("C" => 1, "T" => 2), r_max = 2, active = nothing,
     tab = Dict((1, "C") => (2, 1, "populated"), (1, "T") => (1, 1, "short"), (2, "C") => (0, 0, "empty"), (2, "T") => (0, 0, "empty")), sets = Dict()),
]

@testset "D94e draw on hand tables: ids, pools, statuses, used set (O23a)" begin
    @test length(D94_HAND) >= 12
    for (k, f) in enumerate(D94_HAND)
        seed = d94_seed('e', k)
        act = f.active === nothing ? (1:f.r_max) : f.active
        d94_prefetch_many([d94_round_text(seed, 1, c, r) for r in 1:f.r_max, c in ("S", "C", "T")])
        res = EX.draw_pool_rounds(seed, 1, f.bits; quota = f.q, r_max = f.r_max, active = act)
        ref = d94_ref_pool_rounds(seed, 1, f.bits, f.q, f.r_max; active = act)
        label = f.name
        # the hand table
        for ((r, c), (pool, nd, st)) in f.tab
            cell = d94_cell(res, r, c)
            @test (label, pool === nothing ? nothing : cell.pool, length(cell.drawn), cell.status) == (label, pool, nd, st)
            @test (label, cell.q) == (label, r in act ? get(f.q, c, 0) : 0)
        end
        for ((r, c), ids) in f.sets
            @test (label, sort(d94_cell(res, r, c).drawn)) == (label, ids)
        end
        # classes the table does not mention and that have no size are not drawn
        for r in 1:f.r_max, c in ("S", "C", "T")
            get(f.q, c, 0) == 0 && !haskey(f.tab, (r, c)) && @test (label, d94_cell(res, r, c).status, isempty(d94_cell(res, r, c).drawn)) == (label, "not drawn", true)
        end
        # the independent reference: ids in draw order, pools, sizes, statuses
        for rc in ref
            cell = d94_cell(res, rc.round, rc.class)
            @test (label, rc.round, rc.class, cell.q, cell.drawn, cell.status) == (label, rc.round, rc.class, rc.q, rc.drawn, rc.status)
            rc.q > 0 && @test (label, cell.pool) == (label, rc.pool)
        end
        # I29 and the draw rules: a point is drawn once whatever its class bits; it carries the bit of its class; never more than the size
        all_ids = vcat([c.drawn for c in res]...)
        @test (label, allunique(all_ids)) == (label, true)
        for c in res
            @test all(j -> f.bits[j] & D94_BIT[c.class] != 0, c.drawn)
            @test length(c.drawn) <= c.q
            c.q > 0 && @test (label, c.status) == (label, c.q == length(c.drawn) ? "populated" : isempty(c.drawn) ? "empty" : "short")
        end
        # a short or empty cell stays short or empty once its pool is spent, and is never filled from another class
        for c in ("S", "C", "T")
            sts = [d94_cell(res, r, c).status for r in 1:f.r_max if r in act && get(f.q, c, 0) > 0]
            for j in 1:length(sts), jj in (j + 1):length(sts)
                sts[j] in ("short", "empty") && @test (label, d94_cell(res, jj, c).status) in ((label, "empty"),)
            end
        end
    end
    # the draw of a round cell is a function of the stream and the candidate list only: the same call twice gives the same cells
    f = D94_HAND[12]
    @test EX.draw_pool_rounds(d94_seed('e', 12), 1, f.bits; quota = f.q, r_max = f.r_max) == EX.draw_pool_rounds(d94_seed('e', 12), 1, f.bits; quota = f.q, r_max = f.r_max)
end

@testset "D94e draw in the builder equals the independent reference, fixtures A and B" begin
    for b in (d94_build("A", D94_FA), d94_build("B", D94_FB))
        fx = b.fx
        ref = d94_check_panel(b)
        orbs = [o for o in b.panel["orbit"] if !o["control"]]
        pooled = [o for o in orbs if haskey(o, "point_id")]
        @test all(o -> o["class"] in ("T", "C", "S"), pooled) && all(o -> haskey(o, "point_id") == (o["class"] in ("T", "C", "S")), orbs)
        for i in eachindex(fx.K)
            ids = [o["point_id"] for o in pooled if o["k_index"] == i]
            @test allunique(ids)                                                           # no point id drawn twice at a K (I29)
            bits = d94_side_bits(b.pools, i)
            @test all(o -> bits[o["point_id"]] & D94_BIT[o["class"]] != 0, [o for o in pooled if o["k_index"] == i])      # the class bit of the class
            @test all(o -> 1 <= o["point_id"] <= fx.M, [o for o in pooled if o["k_index"] == i])
        end
        @test allunique([ex_start(o) for o in orbs])                                       # no start twice across rounds, classes and K (I29)
        cells = b.panel["manifest"]["round_cells"]
        for i in eachindex(fx.K), c in ("T", "C", "S")
            rows = sort([x for x in cells if x["k_index"] == i && x["class"] == c]; by = x -> x["round"])
            for j in 1:length(rows), jj in (j + 1):length(rows)
                rows[j]["status"] in ("short", "empty") && @test rows[jj]["drawn"] == 0              # spent pool: never filled again, from no other class
            end
        end
        @test length(cells) == fx.r_max * length(fx.K) * 5                                 # one record per (round, K index, class), undefined and not drawn included
        # the used set in the file: the pool of every round cell is the sidecar pool minus the ids drawn before it at that K, in the order of E34
        for i in eachindex(fx.K), r in 1:fx.r_max, c in ("S", "C", "T")
            row = only(x for x in cells if (x["round"], x["k_index"], x["class"]) == (r, i, c))
            row["status"] in ("undefined", "not drawn") && continue
            earlier = Set(o["point_id"] for o in pooled if o["k_index"] == i && (o["round"] < r || (o["round"] == r && findfirst(==(o["class"]), ("S", "C", "T")) < findfirst(==(c), ("S", "C", "T")))))
            @test row["pool"] == count(j -> !(j in earlier), d94_pool(d94_side_bits(b.pools, i), c))
        end
    end
end

# ------------------------------------------------------------------------------------------------------------------ D94f, D94g
d94_key(o) = (o["round"], o["k_index"], o["class"], o["draw_index"])
d94_sig(o) = (o["class"], o["round"], o["id"], o["orbit_seed"], o["in_D"], o["x0"], o["p0"], o["k_index"], o["draw_index"])

@testset "D94f prefix: the panel of r_max = m restricted to the rounds up to r equals the panel of r_max = r (I31)" begin
    full = d94_build("A", D94_FA)
    for r in (1, 2)
        sm = d94_build("A-r$r", merge(D94_FA, (r_max = r,)))
        small = [o for o in sm.panel["orbit"] if !o["control"]]
        big = Dict(d94_key(o) => o for o in full.panel["orbit"] if !o["control"])
        @test !isempty(small) && all(o -> o["round"] <= r, small)
        @test all(o -> haskey(big, d94_key(o)) && d94_sig(o) == d94_sig(big[d94_key(o)]), small)       # start, class, id, seed, round, D flag: bitwise
        @test length(small) == count(o -> o["round"] <= r, values(big))
        @test sm.panel["manifest"]["r_max"] == r
        @test sm.pools["k"] == full.pools["k"]                                              # the screening does not depend on r_max
    end
    # the streams are independent: a smaller size for one cell changes no start of another cell (C of K index 3 is drawn after S, before T)
    q2 = deepcopy(D94_FA.quota)
    q2[3]["C"] = 2
    other = d94_build("A-q", merge(D94_FA, (quota = q2,)))
    a = Dict(d94_key(o) => ex_start(o) for o in full.panel["orbit"] if !o["control"])
    b = Dict(d94_key(o) => ex_start(o) for o in other.panel["orbit"] if !o["control"])
    @test all(k -> a[k] == b[k], [k for k in keys(a) if !(k[2] == 3 && k[3] == "C") && haskey(b, k)])
    @test all(k -> haskey(b, k), [k for k in keys(a) if !(k[2] == 3 && k[3] == "C")])                  # the other cells keep their orbits
    @test all(k -> a[k] == b[k], [(1, 3, "C", t) for t in 1:2])                                        # round 1 of the cell: the draws t = 1, 2 are the same uniforms
    @test count(k -> k[2] == 3 && k[3] == "C" && k[1] == 1, keys(b)) == 2
end

@testset "D94g K = 0 gates and ids" begin
    for b in (d94_build("A", D94_FA), d94_build("B", D94_FB))
        fx = b.fx
        orbs = [o for o in b.panel["orbit"] if !o["control"]]
        k0 = [o for o in orbs if o["k_index"] == 1]
        @test !isempty(k0) && Set(o["round"] for o in k0) == Set(fx.gates)                     # K = 0 orbits exist in the gate rounds only
        for r in fx.gates, c in ("T", "R")
            @test count(o -> o["round"] == r && o["class"] == c, k0) == min(fx.quota[1][c], count(o -> o["round"] == r && o["class"] == c, k0)) > 0
        end
        @test all(o -> o["class"] in ("T", "R"), k0)
        cells = b.panel["manifest"]["round_cells"]
        for r in 1:fx.r_max, c in D94_CLASSES
            row = only(x for x in cells if (x["round"], x["k_index"], x["class"]) == (r, 1, c))
            r in fx.gates || @test row["q"] == 0 && row["drawn"] == 0 && row["status"] in ("not drawn", "undefined")
            c in ("C", "S", "X") && @test row["status"] == "undefined"
        end
        ids = [o["id"] for o in b.panel["orbit"]]
        @test ids == collect(1:length(ids))                                                    # 1 to N without a gap
        key(o) = (o["round"], o["k_index"], findfirst(==(o["class"]), D94_CLASSES), o["draw_index"])
        @test issorted(orbs; by = key) && allunique(key.(orbs))                                  # the order of E36: round, K index, class, draw index
    end
    # the id and the seed of every orbit of round r are equal for r_max of 2 and 3 (the build with r_max = 2 is the prefix of the build with 3)
    full, two = d94_build("A", D94_FA), d94_build("A-r2", merge(D94_FA, (r_max = 2,)))
    @test [(o["id"], o["orbit_seed"]) for o in two.panel["orbit"] if !o["control"]] == [(o["id"], o["orbit_seed"]) for o in full.panel["orbit"] if !o["control"] && o["round"] <= 2]
end

# ------------------------------------------------------------------------------------------------------------------ D94h (panel part)
@testset "D94h D in the panel: flags are a function of the seed, the cell, the number drawn and the share" begin
    base = d94_build("A", D94_FA)
    ref = d94_ref_panel(base)
    @test count(o -> o["in_D"], [o for o in base.panel["orbit"] if !o["control"]]) == count(o -> o.in_D, ref.orbits)
    for c in ref.cells                                              # the number of D orbits per round cell is floor(f n + 1/2) (O23c)
        c.status in ("undefined", "not drawn") && continue
        f = d94_rat(D94_FA.dfrac[c.k_index][c.class])
        @test count(o -> (o["round"], o["k_index"], o["class"]) == (c.round, c.k_index, c.class) && o["in_D"], base.panel["orbit"]) == d94_dsize(f, c.drawn)
    end
    flags(b) = Dict(d94_key(o) => o["in_D"] for o in b.panel["orbit"] if !o["control"])
    starts(b) = Dict(d94_key(o) => ex_start(o) for o in b.panel["orbit"] if !o["control"])
    f0, s0 = flags(base), starts(base)
    # another fraction of the same cell changes the flags of that cell only and no start; the drawn orbits are unchanged
    d2 = deepcopy(D94_FA.dfrac)
    d2[2]["C"] = "0.7"
    b2 = d94_build("A-df", merge(D94_FA, (dfrac = d2,)))
    @test starts(b2) == s0
    @test all(k -> flags(b2)[k] == f0[k], [k for k in keys(f0) if !(k[2] == 2 && k[3] == "C")])
    @test [o["id"] for o in b2.panel["orbit"]] == [o["id"] for o in base.panel["orbit"]]
    # the fraction of another cell, db_size and controls change no D flag
    d3 = deepcopy(D94_FA.dfrac)
    d3[3]["X"] = "0.9"
    b3 = d94_build("A-df3", merge(D94_FA, (dfrac = d3,)))
    @test starts(b3) == s0 && all(k -> flags(b3)[k] == f0[k], [k for k in keys(f0) if !(k[2] == 3 && k[3] == "X")])
    b4 = d94_build("A-db", merge(D94_FA, (db = 2,)))
    @test starts(b4) == s0 && flags(b4) == f0
    b5 = d94_build("A-ctrl", merge(D94_FA, (controls = "standard",)))
    @test starts(b5) == s0 && flags(b5) == f0
    @test [o for o in b5.panel["orbit"] if !o["control"]] == [o for o in base.panel["orbit"] if !o["control"]]       # the drawn orbits are bit for bit those of the base
end

# ------------------------------------------------------------------------------------------------------------------ D94j
@testset "D94j controls, counts and no screening value" begin
    forbidden = Set(["lambda", "lambda_scr", "block_exponents", "blocks", "e_med", "screening", "s_tenth", "s_n", "exponent_scr", "window"])
    a, bb, c = d94_build("A", D94_FA), d94_build("B", D94_FB), d94_build("A-ctrl", merge(D94_FA, (controls = "standard",)))
    for b in (a, bb, c)
        man = b.panel["manifest"]
        orbs = b.panel["orbit"]
        @test man["totals"]["orbits"] == length(orbs) && man["totals"]["drawn"] == count(o -> !o["control"], orbs) && man["totals"]["controls"] == count(o -> o["control"], orbs)
        @test man["controls"] == [o["id"] for o in orbs if o["control"]]
        @test isdisjoint(Set(lowercase.(ex_all_keys(b.panel))), forbidden) && isdisjoint(Set(lowercase.(ex_all_keys(b.pools))), forbidden)     # I34
        for cl in b.panel["manifest"]["round_cells"]                                            # the counts per round cell equal the statuses and the orbits
            n = count(o -> (get(o, "round", 0), o["k_index"], o["class"]) == (cl["round"], cl["k_index"], cl["class"]), orbs)
            @test cl["drawn"] == n
            @test cl["status"] == (cl["status"] in ("undefined", "not drawn") ? cl["status"] : cl["q"] == n ? "populated" : n > 0 ? "short" : "empty")
        end
    end
    # controls = "none": no control, no orbit of design class control, the check of the base holds trivially
    for b in (a, bb)
        @test isempty(b.panel["manifest"]["controls"]) && !any(o -> o["control"] || o["class"] == "control", b.panel["orbit"])
        @test b.panel["manifest"]["totals"]["controls"] == 0
    end
    # controls = "standard": the controls of the base follow the drawn orbits, the Y5 control last
    ctrl = [o for o in c.panel["orbit"] if o["control"]]
    @test !isempty(ctrl) && ctrl[end]["control_kind"] == "y5" && count(o -> o["control_kind"] == "y5", ctrl) == 1
    @test [o["id"] for o in ctrl] == collect((length(c.panel["orbit"]) - length(ctrl) + 1):length(c.panel["orbit"]))
    @test length(ctrl) <= ex_block_default()["n_ctrl_max"] + 1 && all(o -> o["in_SD"] && !o["in_SL"], ctrl)
    # the manifest of the round panel: source, points, sizes, rounds, gates, db_size, the extended stream rule, the pool counts (counts only)
    man = a.panel["manifest"]
    @test man["screen_source"] == "random" && man["screen_points"] == D94_FA.M && man["r_max"] == 3 && man["gate_rounds"] == [1] && man["db_size"] == D94_FA.db
    @test man["quota"]["2"]["C"] == 3 && man["quota"]["1"]["T"] == 3 && haskey(man, "d_fraction") && haskey(man["d_fraction"], "3")
    for t in ("chaoslab-screen-1|<seed>|<K index>", "chaoslab-panel-1|<seed>|<K index>|<class>|<round>", "chaoslab-d-1|<seed>|<round>|<K index>|<class>")
        @test occursin(t, man["stream_rule"])
    end
    for i in 1:3, cl in ("T", "C", "S")
        @test only(x for x in man["pools"] if x["k_index"] == i && x["class"] == cl)["count"] == count(x -> x & D94_BIT[cl] != 0, d94_side_bits(a.pools, i))
    end
    @test all(x -> sort(collect(keys(x))) == ["class", "count", "k_index"], man["pools"])
end

# ------------------------------------------------------------------------------------------------------------------ D94p, D94c (threads)
@testset "D94p determinism: the panel file, the sidecar and the points for 1, 4 and 23 threads and any screening chunk" begin
    base = d94_build("A", D94_FA)
    h0, hp = ex_sha256(read(base.path)), ex_sha256(read(base.path * ".pools.toml"))
    for (tag, kw) in (("A-nt1", (; ntasks = 1)), ("A-nt23", (; ntasks = 23)), ("A-ch7", (; chunk = 7)), ("A-rerun", (;)))
        b = d94_build(tag, D94_FA; kw...)
        @test ex_sha256(read(b.path)) == h0 && ex_sha256(read(b.path * ".pools.toml")) == hp
    end
    @test !occursin("threads", read(base.path, String)) && !haskey(base.panel["manifest"], "threads")
    side = TOML.parsefile(d94_build("A-nt23", D94_FA; ntasks = 23).path * ".provenance.toml")
    @test side["ntasks"] == 23 && haskey(side, "threads")                                  # the thread count is in the provenance sidecar only
end

# ------------------------------------------------------------------------------------------------------------------ D94t
@testset "D94t the suite runs with a second fixture that differs in every value" begin
    @test D94_FA.seed != D94_FB.seed && D94_FA.K != D94_FB.K && D94_FA.M != D94_FB.M && D94_FA.r_max != D94_FB.r_max && D94_FA.gates != D94_FB.gates
    @test D94_FA.quota != D94_FB.quota && D94_FA.dfrac != D94_FB.dfrac && D94_FA.db != D94_FB.db
    b = d94_build("B", D94_FB)
    ref = d94_check_panel(b)                                            # the generic checks of the helpers: hand-independent reference from the fixture
    @test length(ref.orbits) == count(o -> !o["control"], b.panel["orbit"]) > 0
    @test b.panel["manifest"]["r_max"] == 2 && b.panel["manifest"]["gate_rounds"] == [2] && b.panel["manifest"]["screen_points"] == 48
    @test all(o -> o["run_seed"] == "777", b.panel["orbit"])
    # no key of the grammar can hold a threshold or a target: checked in D94a (refused as unknown keys); the source search is a stage 04 script
end

# ------------------------------------------------------------------------------------------------------------------ guards of the references
@testset "D94 references: the independent pool draw reproduces the hand tables; the panel reference runs without the module (guards)" begin
    for (k, f) in enumerate(D94_HAND)
        seed = d94_seed('e', k)
        act = f.active === nothing ? (1:f.r_max) : f.active
        ref = d94_ref_pool_rounds(seed, 1, f.bits, f.q, f.r_max; active = act)
        for ((r, c), (pool, nd, st)) in f.tab
            cell = d94_cell(ref, r, c)
            @test (f.name, pool === nothing ? nothing : cell.pool, length(cell.drawn), cell.status) == (f.name, pool, nd, st)
        end
        for ((r, c), ids) in f.sets
            @test (f.name, sort(d94_cell(ref, r, c).drawn)) == (f.name, ids)
        end
    end
    # reference screening bits of fixture A, the sidecar form (run-length lists) and the reference of the whole panel, on the manifest of the lattice build
    # of the base fixture (the same K list: its periodic table feeds the R draw)
    fb = ex_build("base")
    bits = [d94_ref_bits((; fx = D94_FA), i) for i in 1:3]
    @test all(b -> length(b) == D94_FA.M, bits) && count(x -> x & 0x01 != 0, bits[1]) > 0 && all(x -> x & 0x06 == 0, bits[1])
    side = Dict("k" => Dict(string(i) => Dict("points" => D94_FA.M, "bits_sha256" => ex_sha256(bits[i]), "T" => d94_runs(bits[i], "T"),
                                              "C" => d94_runs(bits[i], "C"), "S" => d94_runs(bits[i], "S")) for i in 1:3))
    @test all(i -> d94_side_bits(side, i) == bits[i], 1:3)
    ref = d94_ref_panel((; fx = D94_FA, panel = Dict("manifest" => fb.panel["manifest"]), pools = side))
    @test !isempty(ref.orbits) && [o.id for o in ref.orbits] == collect(1:length(ref.orbits)) && length(ref.cells) == 3 * 3 * 5
    @test Set(o.round for o in ref.orbits if o.k_index == 1) == Set([1]) && all(o -> o.class in ("T", "R"), [o for o in ref.orbits if o.k_index == 1])
    @test all(c -> c.drawn <= c.q, ref.cells) && issorted([(o.round, o.k_index, findfirst(==(o.class), D94_CLASSES), o.t) for o in ref.orbits])
    @test allunique([(o.x0, o.p0) for o in ref.orbits if o.class in ("T", "C", "S")]) && count(o -> o.in_D, ref.orbits) > 0
    pts = d94_points(D94_FA.seed, 1, 5)
    @test all(p -> 0 <= p[1] < 2pi && 0 <= p[2] < 2pi, pts) && pts == d94_points(D94_FA.seed, 1, 5) && pts != d94_points(D94_FA.seed, 2, 5)
end
