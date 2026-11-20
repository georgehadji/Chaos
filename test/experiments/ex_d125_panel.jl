# Amendment D125, stage 02: D125d (points and starts), D125e (screening and classes), D125f (draw, D, ids), D125g (controls), D125h (determinism of the panel).
# Function-level names assumed (plan, "New API"): EX.systems_points, systems_flow_starts, screen_floor, screen_ratio, class_bits, cell_quantile, systems_draw,
# systems_dsubset, control_exponents. Builder level goes through EX.build_panel; the file layout is pinned in the plan (U1, U6, U7).

d125_ks(x::AbstractVector) = (n = length(x); s = sort(x); maximum(max(i / n - s[i], s[i] - (i - 1) / n) for i in 1:n))
d125_map_cells(fx) = [c for c in fx.cells if c["kind"] != "henon_heiles"]
d125_win(c) = (c["q2_lo"], c["q2_hi"], c["p2_lo"], c["p2_hi"])

@testset "D125d screening points of the maps (E46)" begin
    for fx in (D125_FA, D125_FB), c in d125_map_cells(fx)
        d = d125_dim(c)
        pts = EX.systems_points(fx.seed, c["id"], d, 1000)                                 # RED: UndefVarError (systems_points)
        ref = d125_points(fx.seed, c["id"], d, 1000)
        @test all(i -> Float64.(collect(pts[i])) == ref[i], 1:1000)                         # tolerance 0
        @test all(p -> all(x -> 0.0 <= x < EX_TWO_PI, p), pts)                              # all coordinates, momenta included, in [0, 2 pi)
        @test length(unique(Tuple.(pts))) == 1000
        @test Float64.(collect(EX.systems_points(fx.seed, c["id"], d, 7)[7])) == ref[7]     # prefix property in M
        @test Float64.(collect(EX.systems_points(fx.seed + 1, c["id"], d, 3)[1])) != ref[1]
        @test Float64.(collect(EX.systems_points(fx.seed, c["id"] * "x", d, 3)[1])) != ref[1]
    end
    # Kolmogorov-Smirnov statistic of every coordinate of 10^4 points of fixture A against the uniform law (the bound of the criterion)
    for c in d125_map_cells(D125_FA)
        p = EX.systems_points(D125_FA.seed, c["id"], d125_dim(c), 10_000)
        @test all(i -> d125_ks([q[i] / EX_TWO_PI for q in p]) <= 1.95 / sqrt(10_000), 1:d125_dim(c))
    end
end

@testset "D125d starts of the flow cell (E47, I49)" begin
    for fx in (D125_FA, D125_FB)
        c = fx.cells[3]
        starts, pairs = EX.systems_flow_starts(fx.seed, c["id"], c["energy"], d125_win(c), 200)    # RED: UndefVarError
        ref, rpairs = d125_flow_starts(fx.seed, c["id"], c["energy"], d125_win(c), 200)
        @test pairs == rpairs && length(starts) == 200
        @test all(i -> starts[i][1] == 0.0 && starts[i][2] == ref[i][2] && starts[i][4] == ref[i][4], 1:200)          # q1 = 0, q2 and p2 bitwise
        @test all(i -> starts[i][3] > 0 && abs(starts[i][3] - ref[i][3]) <= 2 * eps(ref[i][3]), 1:200)             # p1 within 2 ulp (reading U3)
        E = setprecision(() -> d125_big(c["energy"]), BigFloat, 256)
        defect = i -> setprecision(BigFloat, 256) do
            q2, p1, p2 = BigFloat(starts[i][2]), BigFloat(starts[i][3]), BigFloat(starts[i][4])
            abs((p1^2 + p2^2) / 2 + q2^2 / 2 - q2^3 / 3 - E)
        end
        bound = 16 * ldexp(1.0, -52) * max(1.0, Float64(E))
        @test all(i -> defect(i) <= bound, 1:200)                                                                  # I49
        @test rpairs >= 200
    end
    # a window without energy-surface points is refused after 1000 M pairs (R29 item 3 refusal class of E47), also through the builder (0 files)
    c = Dict{String, Any}(D125_FA.cells[3]); c["q2_lo"] = "0.9"; c["q2_hi"] = "0.95"; c["p2_lo"] = "0.5"; c["p2_hi"] = "0.6"; c["n_start"] = 3
    @test ex_rejects(() -> EX.systems_flow_starts(EX_SEED, "hhx", c["energy"], d125_win(c), 3), ["too many rejections"])
    blk = d125_edit(d125_block(D125_FA), ("cell", 3) => c)
    @test d125_refused(blk, ["too many rejections", "hh1"]; tag = "empty-window")
end

@testset "D125e class predicates at and one step around every rule value (E48, E49)" begin
    rules = d125_rules_nt(D125_FA.rules)
    t = 400; phi = d125_phi(D125_FA.rules["floor_abs"], D125_FA.rules["floor_log"], t)
    @test EX.screen_floor(D125_FA.rules["floor_abs"], D125_FA.rules["floor_log"], t) === phi                                # RED: UndefVarError
    # hand cases that do not use the reference (strict for regular, not strict for chaotic at phi)
    cb = EX.class_bits
    @test cb(prevfloat(phi), rules.ratio_regular_max, phi, rules, nothing) == 0x01
    @test cb(phi, rules.ratio_regular_max, phi, rules, nothing) == 0x00                                                      # lambda = phi is not regular, ratio is not chaotic
    @test cb(phi, rules.ratio_chaotic_min, phi, rules, nothing) & 0x02 == 0x02                                              # not strict for chaotic
    @test cb(prevfloat(phi), rules.ratio_chaotic_min, phi, rules, nothing) == 0x00
    @test cb(10 * phi, nextfloat(rules.ratio_regular_max), phi, rules, nothing) == 0x00                                      # ambiguous ratio gap
    @test cb(10 * phi, rules.ratio_weak_min, phi, rules, nothing) == (0x02 | 0x04)                                           # weak by the ratio
    @test cb(10 * phi, prevfloat(rules.ratio_weak_min), phi, rules, nothing) == (0x02 | 0x08)                                # strong
    @test cb(10 * phi, 7.0, phi, rules, nothing) == (0x02 | 0x08)                                                            # Lambda_q undefined: second weak test false
    Lq = 8.0 * phi; thr = setprecision(() -> Float64(BigFloat(rules.weak_fraction) * BigFloat(Lq)), BigFloat, 256)
    @test cb(thr, 7.0, phi, rules, Lq) == (0x02 | 0x04) && cb(nextfloat(thr), 7.0, phi, rules, Lq) == (0x02 | 0x08) && cb(prevfloat(thr), 7.0, phi, rules, Lq) == (0x02 | 0x04)
    # no ratio: G(N / back) <= 0 or not finite gives no ratio and bits 0 (ambiguous)
    @test EX.screen_ratio(1.0, 0.0) === nothing && EX.screen_ratio(1.0, -2.0) === nothing && EX.screen_ratio(1.0, NaN) === nothing && EX.screen_ratio(NaN, 1.0) === nothing
    @test EX.screen_ratio(6.0, 2.0) == 3.0
    @test cb(10 * phi, nothing, phi, rules, Lq) == 0x00
    # 120 synthetic records against the reference, every rule value and one step either side
    rs = vcat([[prevfloat(v), v, nextfloat(v)] for v in (rules.ratio_regular_max, rules.ratio_chaotic_min, rules.ratio_weak_min)]..., [1.0, 4.0, 7.5, 12.0, nothing])
    ls = [prevfloat(phi), phi, nextfloat(phi), prevfloat(thr), thr, nextfloat(thr), 0.5 * phi, 20 * phi]
    n = 0
    for Lqv in (nothing, Lq), r in rs, l in ls
        n += 1
        @test cb(l, r, phi, rules, Lqv) == d125_bits(l, r, phi, rules, Lqv)
    end
    @test n >= 80
    # phi equals the 1024-bit value rounded once (0 ulp) for 50 values of t
    rng = StableRNG(d125_seed('e', 2))
    for _ in 1:50
        tt = rand(rng, 2:2_000_000)
        @test EX.screen_floor(D125_FA.rules["floor_abs"], D125_FA.rules["floor_log"], tt) === d125_phi(D125_FA.rules["floor_abs"], D125_FA.rules["floor_log"], tt)
    end
end

@testset "D125e Lambda_q by the nearest rank (E49, P134)" begin
    for M in (1, 2, 9, 10, 11), q in ("9/10", "1/2", "0.25")
        lams = collect(1.0:M)                                                                 # value = rank
        qr = d125_rat(q)
        @test EX.cell_quantile(shuffle(StableRNG(d125_seed('e', 3)), lams), q) == Float64(cld(numerator(qr) * M, denominator(qr)))   # hand rank ceil(q M)
    end
    @test EX.cell_quantile(Float64[], "9/10") === nothing                                      # M_c = 0: undefined
    @test EX.cell_quantile(collect(1.0:10), "9/10") == 9.0 && EX.cell_quantile(collect(1.0:11), "9/10") == 10.0 && EX.cell_quantile(collect(1.0:9), "9/10") == 9.0
    rng = StableRNG(d125_seed('e', 4))
    for M in (1, 2, 9, 10, 11, 57), q in ("9/10", "7/10", "1/3")
        l = rand(rng, M)
        @test EX.cell_quantile(l, q) == d125_Lq(l, q)
    end
end

@testset "D125e class bits and cell values of 200 points per cell kind equal the verified Benettin (tolerance 0 on the bits)" begin
    for fx in (D125_FA, D125_FB)
        blk = d125_big_block(fx, 200)
        b = d125_build("big-" * fx.name, fx; block = blk)                                     # RED: ArgumentError (unknown key scheme)
        cells = Dict(c["id"] => d125_ref_cell(blk, c) for c in blk["cell"])
        d125_check_panel(b; cells = cells)
        @test all(c -> length(cells[c["id"]].bits) == 200, blk["cell"])
    end
end

@testset "D125f draw and subset D on hand fixtures (E51, E53, I44, I52)" begin
    seed = EX_SEED
    pools = [Int[], [5], [2, 9], collect(1:12), collect(1:50), collect(1:7), collect(10:10:100), [3, 4, 5, 6], collect(1:300), [17], collect(1:9), collect(1:11)]
    sizes = [3, 3, 3, 0, 5, 7, 4, 10, 6, 1, 9, 2]
    for (i, (pool, size)) in enumerate(zip(pools, sizes))
        got = EX.systems_draw(seed, "y3a", "reg", pool, size)                                  # RED: UndefVarError
        @test got == d125_draw(seed, "y3a", "reg", pool, size)
        @test length(got) == min(size, length(pool)) && length(unique(got)) == length(got) && all(in(pool), got)
        size >= length(pool) && @test sort(got) == pool                                         # short or complete: never filled from elsewhere
    end
    @test EX.systems_draw(seed, "y3a", "reg", collect(1:50), 5) != EX.systems_draw(seed, "y3a", "wk", collect(1:50), 5)     # the class id is in the stream text
    @test EX.systems_draw(seed, "y3a", "reg", collect(1:50), 5) != EX.systems_draw(seed, "ch2", "reg", collect(1:50), 5)
    rng = StableRNG(d125_seed('f', 1))
    for _ in 1:20                                                                                # 20 fixtures, n up to 300, d_min above n included
        n = rand(rng, 0:300); f = rand(rng, ("0", "1/3", "0.4", "1/2", "7/10", "1")); dmin = rand(rng, (0, 1, 5, 400))
        got = EX.systems_dsubset(seed, "ch2", "st", n, f, dmin)
        @test got == d125_dindices(seed, "ch2", "st", n, f, dmin)
        @test length(got) == min(n, max(dmin, d94_dsize(d125_rat(f), n))) && issorted(got) && all(in(1:n), got)
    end
    @test EX.systems_dsubset(seed, "ch2", "st", 3, "0", 10) == [1, 2, 3]                       # d_min above n: all
end

@testset "D125f built panel equals the reference draw, ids, D, seeds and statuses; streams of other cells do not interact" begin
    b = d125_build("big-A", D125_FA; block = d125_big_block(D125_FA, 200))
    orbs = d125_check_panel(b)
    @test length(unique((o.cell, o.point_id) for o in orbs)) == length(orbs)                    # no point twice (I44)
    @test length(unique(Tuple(o.start) for o in orbs)) == length(orbs)                          # no start twice
    @test all(c -> sum(1 for o in orbs if o.cell == c["id"]) <= sum(values(c["quota"])), b.block["cell"])
    # D flags unchanged when the quota and fraction of another cell change (stream independence)
    blk2 = d125_edit(b.block, ("cell", 2, "quota") => Dict("reg" => 2, "wk" => 2, "st" => 3), ("cell", 2, "d_fraction") => Dict("reg" => "1/2", "wk" => "1", "st" => "1/2"))
    b2 = d125_build("big-A2", D125_FA; block = blk2)
    sel = (b, cell) -> [(o["class"], o["point_id"], o["in_D"]) for o in b.panel["orbit"] if get(o, "cell", "") == cell]
    @test sel(b, "y3a") == sel(b2, "y3a") && sel(b, "hh1") == sel(b2, "hh1")
end

@testset "D125g controls (E54, R32)" begin
    ctl = [Dict{String, Any}("id" => "a", "kind" => "y5", "dimension" => 4, "mu" => ["3", "5/2"]),
           Dict{String, Any}("id" => "b", "kind" => "y5", "dimension" => 6, "mu" => ["2", "3/2", "4"]),
           Dict{String, Any}("id" => "c", "kind" => "chain_origin", "N" => 3, "K" => "1.5", "boundary" => "ring", "eps" => "0"),
           Dict{String, Any}("id" => "d", "kind" => "chain_origin", "N" => 2, "K" => "0.7", "boundary" => "open", "eps" => "0")]
    for c in ctl
        ex = EX.control_exponents(c)                                                           # RED: UndefVarError
        d = c["kind"] == "y5" ? c["dimension"] : 2 * c["N"]
        @test length(ex) == d && issorted(ex; rev = true) && ex == -reverse(ex)                 # d values, descending, symmetric
        ref = setprecision(BigFloat, 1024) do                                                   # independent 1024-bit values
            if c["kind"] == "y5"
                v = BigFloat[]; for m in c["mu"]; x = d125_big(m); push!(v, log(x), -log(x)); end
                sort(v; rev = true)
            else
                k = d125_big(c["K"]); x = log((2 + k + sqrt((2 + k)^2 - 4)) / 2)
                vcat(fill(x, c["N"]), fill(-x, c["N"]))
            end
        end
        @test all(i -> abs(ex[i] - ref[i]) <= 1e-16 * abs(ref[i]), 1:d)
    end
    @test count(==(EX.control_exponents(ctl[3])[1]), EX.control_exponents(ctl[3])) == 3          # multiplicity N of +x
    @test ex_rejects(() -> EX.control_exponents(Dict{String, Any}("id" => "z", "kind" => "y5", "dimension" => 3, "mu" => ["3"])), ["dimension"])
    b = d125_build("big-A", D125_FA; block = d125_big_block(D125_FA, 200))
    ctls = [o for o in b.panel["orbit"] if o["control"]]
    @test length(ctls) == 2 && ctls[1]["control_kind"] == "y5" && ctls[2]["control_kind"] == "chain_origin"
    @test Float64.(ctls[1]["exact_exponents"]) == EX.control_exponents(D125_FA.controls[1]) && Float64.(ctls[2]["exact_exponents"]) == EX.control_exponents(D125_FA.controls[2])
    @test all(o -> o["start"] == zeros(4) && o["dimension"] == 4, ctls)                          # start the origin
    @test [m["id"] for m in get(b.panel["manifest"], "control", [])] == ["c_y5", "c_ch"]
end

@testset "D125h determinism of the panel (I42) and forbidden content (I43)" begin
    fx = D125_FA
    base = d125_build("det-base", fx)
    sha = p -> ex_sha256(read(p))
    ref = (sha(base.path), sha(base.path * ".pools.toml"))
    for (i, kw) in enumerate(((; ntasks = 1), (; ntasks = 4), (; ntasks = 23), (; chunk = 64), (; chunk = 4096), (; ntasks = 4, chunk = 7)))
        out = ex_path("d125-det-$i-panel.toml")
        EX.build_panel(base.prereg; out = out, kw...)
        @test (sha(out), sha(out * ".pools.toml")) == ref
    end
    out = ex_path("d125-det-again-panel.toml"); EX.build_panel(base.prereg; out = out)
    @test sha(out) == ref[1]                                                                    # a second build
    # thread, backend and device data only in the provenance sidecar; the screening ran on the CPU with an empty reason (R30 item 1)
    side = TOML.parsefile(base.path * ".provenance.toml")
    @test side["backend"] == "cpu" && get(side, "backend_fallback_reason", "") == ""
    forbidden = ("wall", "wall_seconds", "time", "host", "hostname", "threads", "ntasks", "backend", "device", "lambda", "lambdas", "ratio", "ratios", "block_value", "created_utc")
    function keys_of(x, acc = String[])
        x isa AbstractDict && for (k, v) in x; push!(acc, lowercase(k)); keys_of(v, acc); end
        x isa AbstractVector && foreach(v -> keys_of(v, acc), x)
        return acc
    end
    @test isempty(intersect(keys_of(base.panel), forbidden)) && isempty(intersect(keys_of(base.pools), forbidden))
    @test sum(count(==("lambda_q"), keys_of(t)) for t in (base.panel, base.pools)) == length(fx.cells)   # the only screening-derived value (E50), once per cell
    @test !haskey(base.panel["manifest"], "K") && !any(o -> haskey(o, "x0"), base.panel["orbit"])
    # a second fixture that differs in every value passes the generic checks (D125t of the grammar group)
    bB = d125_build("det-B", D125_FB)
    d125_check_panel(bB)
end
