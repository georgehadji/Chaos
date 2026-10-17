# Panel-level tests on small fixture builds (3 K values, 8 x 8 grid, N_scr = 2000, q_max = 2): AC12 (cells), AC13, AC15, AC16,
# AC17, AC18, AC21, AC22, AC23, AC24, AC29 and the invariants I1 to I10. The full study panel is not built here (stage 04).
# File layout (fixed by these tests, see the plan): [manifest] with sub-tables, one [[orbit]] table per orbit.

const EX_CLS = Dict("T" => 1, "R" => 2, "C" => 3, "S" => 4, "X" => 5)
const EX_BLOCK = ex_block()
ex_man(b) = b.panel["manifest"]
ex_cells(b) = ex_man(b)["cells"]
ex_cell(b, i, k) = only(c for c in ex_cells(b) if c["k_index"] == i && c["class"] == k)
ex_of(b, i, k) = [o for o in b.panel["orbit"] if get(o, "k_index", 0) == i && o["class"] == k]
ex_ctrls(b) = [o for o in b.panel["orbit"] if o["control"]]
ex_drawn(b) = [o for o in b.panel["orbit"] if !o["control"]]
ex_sig(s) = length(lstrip(replace(split(lowercase(s), 'e')[1], r"[-+.]" => ""), '0'))
ex_bigdec(s, bits = 256) = setprecision(() -> parse(BigFloat, s), BigFloat, bits)
ex_wrapf(x) = (y = x - EX_TWO_PI * floor(x / EX_TWO_PI); y)

@testset "panel structure and invariants (I1 to I6, I10)" begin
    b = ex_build("base")
    orbs, man = b.panel["orbit"], ex_man(b)
    N = length(orbs)
    @test N >= 5                                                  # at least the controls
    @test [o["id"] for o in orbs] == collect(1:N)                 # I2: ids 1..N without gaps, in the order of E20
    key(o) = (o["k_index"], EX_CLS[o["class"]], o["draw_index"])
    dr = ex_drawn(b)
    @test issorted(dr, by = key) && allunique(key.(dr))
    @test all(o -> o["control"] == (o["class"] == "control"), orbs)
    ids_c = [o["id"] for o in ex_ctrls(b)]
    @test ids_c == collect((N - length(ids_c) + 1):N)             # controls follow the drawn orbits
    kinds = [o["control_kind"] for o in ex_ctrls(b)]
    kind_rank = Dict("fixed" => 1, "newton" => 2, "y5" => 3)
    @test issorted(kinds, by = k -> kind_rank[k])
    @test kinds[end] == "y5" && count(==("y5"), kinds) == 1
    # I1: finite lifted starts, wrapped form in [0, 2 pi)
    @test all(o -> isfinite(o["x0"]) && isfinite(o["p0"]), orbs)
    @test all(o -> 0 <= ex_wrapf(o["x0"]) < EX_TWO_PI && 0 <= ex_wrapf(o["p0"]) < EX_TWO_PI, orbs)
    # I2: membership, class equals the class of the cell it was drawn for; draw index runs 1..drawn within a cell
    for c in ex_cells(b)
        os = ex_of(b, c["k_index"], c["class"])
        @test [o["draw_index"] for o in os] == collect(1:length(os))
        @test length(os) == c["drawn"]
        @test c["drawn"] <= EX_BLOCK["n_c"]
        # I5: no cell filled from another class; the status words
        c["status"] == "populated" && @test c["drawn"] == EX_BLOCK["n_c"] && c["reason"] == ""
        c["status"] == "short" && @test 0 < c["drawn"] < EX_BLOCK["n_c"] && c["reason"] != ""
        c["status"] == "empty" && @test c["drawn"] == 0 && c["reason"] in ("no qualifying grid cell", "no elliptic periodic point")
        c["status"] == "undefined" && @test c["drawn"] == 0 && c["reason"] == "class undefined at K = 0"
    end
    # I10: at K = 0 (index 1) only T and R exist; C, S and X are recorded undefined, not empty
    @test sort([c["class"] for c in ex_cells(b) if c["k_index"] == 1]) == ["C", "R", "S", "T", "X"]
    @test all(c -> c["status"] == "undefined", (ex_cell(b, 1, k) for k in ("C", "S", "X")))
    @test all(c -> c["status"] != "undefined", (ex_cell(b, 1, k) for k in ("T", "R")))
    @test isempty([o for o in orbs if get(o, "k_index", 0) == 1 && o["class"] in ("C", "S", "X")])
    @test all(i -> all(k -> ex_cell(b, i, k)["status"] != "undefined", ("T", "R", "C", "S", "X")), 2:3)
    # I3: controls are disjoint from drawn starts after wrapping and appear once
    # keyed by K index: the same fixed point at two K values is two controls (AC18)
    w(o) = (get(o, "k_index", 0), ex_wrapf(o["x0"]), ex_wrapf(o["p0"]))
    cs = w.(ex_ctrls(b))
    @test allunique(cs) && isdisjoint(cs, w.(dr))
    # I6 and the number of controls
    @test length(ids_c) <= EX_BLOCK["n_ctrl_max"] + 1
    sd, sl = [o for o in orbs if o["in_SD"]], [o for o in orbs if o["in_SL"]]
    @test all(o -> o["in_SD"], ex_ctrls(b)) && !any(o -> o["in_SL"], ex_ctrls(b))
    @test all(o -> o["in_SD"], sl) && length(sd) <= EX_BLOCK["cap_D"] && length(sl) <= EX_BLOCK["cap_L"]
end

@testset "AC16 draw in the build: starts are screening points" begin
    b = ex_build("base")
    gc = Set((c[1], c[2]) for c in EX.grid_centres(EX_BLOCK["grid_side"]))
    for i in 1:3
        starts = Dict(k => [ex_start(o) for o in ex_of(b, i, k)] for k in ("T", "C", "S"))
        for k in ("T", "C", "S")
            @test all(in(gc), starts[k])
            @test allunique(starts[k])                              # without replacement
        end
        @test isdisjoint(starts["S"], starts["C"])                  # a grid cell drawn for S is removed from the C candidates
    end
    # a cell whose class qualifies nowhere is empty with its reason; the other classes of that K are unchanged (K = 10, q_max = 1)
    e = ex_build("noell", Dict{String, Any}("K" => ["10"], "q_max" => 1))
    @test ex_cell(e, 1, "R")["status"] == "empty" && ex_cell(e, 1, "R")["reason"] == "no elliptic periodic point"
    @test isempty(ex_of(e, 1, "R"))
    @test all(k -> ex_cell(e, 1, k)["status"] != "undefined", ("T", "C", "S", "X"))
    @test length(ex_of(e, 1, "X")) == EX_BLOCK["n_c"]
end

@testset "AC12 K = 0 cells" begin
    # closed-form screening exponent at K = 0 (AC12 above): about 3.6e-3 here; thresholds on each side make the T cell full or empty
    full = ex_build("k0-full", Dict{String, Any}("K" => ["0"]))
    @test ex_cell(full, 1, "T")["status"] == "populated" && ex_cell(full, 1, "T")["drawn"] == EX_BLOCK["n_c"]
    @test ex_cell(full, 1, "T")["qualifying"] == EX_BLOCK["grid_side"]^2
    empty = ex_build("k0-empty", Dict{String, Any}("K" => ["0"], "theta_T" => 0.001))
    @test ex_cell(empty, 1, "T")["status"] == "empty" && ex_cell(empty, 1, "T")["reason"] == "no qualifying grid cell"
    @test ex_cell(empty, 1, "T")["qualifying"] == 0 && isempty(ex_of(empty, 1, "T"))
    @test ex_cell(empty, 1, "R")["status"] == "populated"            # the other class of K = 0 is unchanged
end

@testset "AC13 K = 0 resonant class" begin
    b = ex_build("k0-full", Dict{String, Any}("K" => ["0"]))
    rs = ex_of(b, 1, "R")
    @test length(rs) == EX_BLOCK["n_c"]
    for o in rs
        a, bb = o["a"], o["b"]
        @test gcd(a, bb) == 1 && 1 <= bb <= EX_BLOCK["q_max"] && 0 <= a < bb
        p_ref = setprecision(() -> Float64(2 * BigFloat(pi) * a / bb), BigFloat, EX_BITS)
        @test o["p0"] == p_ref                                       # nearest binary64 to 2 pi a / b
        x = o["x0"]
        for _ in 1:bb
            x = x + o["p0"]                                          # K = 0: x' = x + p, p' = p
        end
        @test abs(x - (o["x0"] + EX_TWO_PI * a)) <= bb * eps(max(EX_TWO_PI * a, x))
    end
    @test sort(unique(c["class"] for c in ex_cells(b))) == ["C", "R", "S", "T", "X"]
    @test all(k -> ex_cell(b, 1, k)["status"] == "undefined", ("C", "S", "X"))
end

@testset "AC15 R starts at K > 0" begin
    b = ex_build("base")
    per = Dict(p["id"] => p for p in ex_man(b)["periodic"])
    dR = ex_man(b)["inputs"]["delta_R"]["value"]
    for i in 2:3
        rs = ex_of(b, i, "R")
        for o in rs
            p = per[o["source_orbit"]]
            @test p["class"] == "elliptic" && p["q"] <= EX_BLOCK["q_max"] && p["k_index"] == i
            pt = p["points"][o["source_point"]]
            setprecision(BigFloat, 256) do
                zx, zp = ex_bigdec(pt[1]), ex_bigdec(pt[2])
                phi = o["phi"]
                x = zx + BigFloat(dR) * cos(BigFloat(phi))
                pp = zp + BigFloat(dR) * sin(BigFloat(phi))
                @test o["x0"] == Float64(x) && o["p0"] == Float64(pp)            # rebuilt bit for bit
                d = ex_tdist((o["x0"], o["p0"]), (Float64(zx), Float64(zp)))
                @test abs(d - dR) <= 1e-12 * dR + 16 * EX_EPS
            end
        end
    end
    # X starts: torus distance r from the hyperbolic fixed point (0, 0), log10 r inside the input range
    for i in 2:3, o in ex_of(b, i, "X")
        @test EX_BLOCK["log10_r"][1] - 1e-12 <= log10(o["r"]) <= EX_BLOCK["log10_r"][2] + 1e-12
        @test abs(ex_tdist((o["x0"], o["p0"]), (0.0, 0.0)) - o["r"]) <= 1e-12 * o["r"]
    end
end

@testset "AC17 seeds in the file" begin
    b = ex_build("base")
    seed = EX_SEED
    for o in ex_drawn(b)
        @test parse(UInt64, o["run_seed"]) == seed
        @test parse(UInt64, o["orbit_seed"]) == EX.orbit_seed(seed, o["id"])
        @test parse(UInt64, o["cell_seed"]) == ex_cell_seed(seed, o["k_index"], o["class"])
        @test o["draw_index"] >= 1
    end
    for o in ex_ctrls(b)
        @test parse(UInt64, o["run_seed"]) == seed && parse(UInt64, o["orbit_seed"]) == EX.orbit_seed(seed, o["id"])
        @test !haskey(o, "cell_seed") && !haskey(o, "draw_index")      # controls are not drawn: no cell, no draw index
    end
    @test ex_man(b)["stream_rule"] isa String && occursin("chaoslab-panel-1", ex_man(b)["stream_rule"])
end

@testset "AC18 controls" begin
    b = ex_build("base")
    cs = ex_ctrls(b)
    fixed = [o for o in cs if o["control_kind"] == "fixed"]
    @test [(o["k_index"], o["x0"], o["p0"]) for o in fixed] == [(2, 0.0, 0.0), (2, Float64(pi), 0.0), (3, 0.0, 0.0), (3, Float64(pi), 0.0)]
    setprecision(BigFloat, EX_BITS) do
        for o in fixed
            K = BigFloat(o["K"])
            t = o["x0"] == 0 ? 2 + K : 2 - K
            want = abs(t) > 2 + 1e-12 ? Float64(log((abs(t) + sqrt(t^2 - 4)) / 2)) : 0.0
            @test o["period"] == 1 && (o["wind_m"], o["wind_j"]) == (0, 0)
            o["exponent"] == 0 ? (@test want == 0) : (@test abs(o["exponent"] - want) <= 1e-16 * want)
            @test o["stability"] == (abs(t) > 2 + 1e-12 ? "hyperbolic" : (abs(t) < 2 - 1e-12 ? "elliptic" : "parabolic"))
        end
    end
    # Y5 control
    y = only(o for o in cs if o["control_kind"] == "y5")
    @test y["system"] == "linear" && y["x0"] == 0.0 && y["p0"] == 0.0
    mu = setprecision(() -> ex_eval(ex_man(b)["inputs"]["mu"]["value"]), BigFloat, EX_BITS)
    @test abs(y["exponent"] - Float64(log(mu))) <= 16 * EX_EPS * Float64(log(mu))
    @test abs(y["mu"] - Float64(mu)) <= 16 * EX_EPS * Float64(mu)
    # the list is cut at n_ctrl_max, fixed points first (n_ctrl_max = 2 keeps the two fixed points of the first K > 0); Y5 is extra
    c2 = ex_build("ctrl2", Dict{String, Any}("n_ctrl_max" => 2))
    @test [(o["control_kind"], o["k_index"]) for o in ex_ctrls(c2)][1:2] == [("fixed", 2), ("fixed", 2)]
    @test length(ex_ctrls(c2)) == 3 && ex_ctrls(c2)[end]["control_kind"] == "y5"
    # beyond the fixed points the Newton orbits follow, ordered by (period, K index, representative), each with its analytic exponent
    c7 = ex_build("ctrl7", Dict{String, Any}("n_ctrl_max" => 7))
    nw = [o for o in ex_ctrls(c7) if o["control_kind"] == "newton"]
    @test length(nw) == 3 && length(ex_ctrls(c7)) == 4 + 3 + 1
    @test issorted([(o["period"], o["k_index"], o["x0"], o["p0"]) for o in nw])
    # the first Newton orbits of the table (fixed points excluded): the controls are a prefix of that order
    rep(p) = (Float64(ex_bigdec(p["points"][1][1])), Float64(ex_bigdec(p["points"][1][2])))
    isfixed(p) = p["q"] == 1 && any(z -> ex_tdist(rep(p), z) < 1e-9, ((0.0, 0.0), (Float64(pi), 0.0)))
    tab = [p for p in ex_man(c7)["periodic"] if p["k_index"] >= 2 && !isfixed(p)]
    order = sort(tab, by = p -> (p["q"], p["k_index"], rep(p)...))
    @test [(p["q"], p["k_index"], rep(p)...) for p in order[1:3]] == [(o["period"], o["k_index"], o["x0"], o["p0"]) for o in nw]
    setprecision(BigFloat, EX_BITS) do
        for o in nw
            K = BigFloat(o["K"])
            z = (ex_bigdec(o["x0_dec"], EX_BITS), ex_bigdec(o["p0_dec"], EX_BITS))
            @test ex_residual(K, z, o["period"]) <= EX_RHO * (1 + ldexp(1.0, -200))
            t = ex_orbit_trace(K, z, o["period"])
            want = abs(t) > 2 + 1e-12 ? Float64(log((abs(t) + sqrt(t^2 - 4)) / 2) / o["period"]) : 0.0
            o["exponent"] == 0 ? (@test want == 0) : (@test abs(o["exponent"] - want) <= 1e-15 * want)
            @test (o["exponent"] > 0) == (o["stability"] == "hyperbolic")
        end
    end
end

@testset "AC21 subsets in the build" begin
    base = ex_build("base")
    sd = [o for o in base.panel["orbit"] if o["in_SD"]]
    @test ex_man(base)["subsets"]["capped_D"] == false && isempty(ex_man(base)["subsets"]["removed_D"])
    # before capping, every populated cell contributes min(s_D, drawn) members of rank 1..s_D (s_L: rank 1..s_L)
    for c in ex_cells(base)
        os = ex_of(base, c["k_index"], c["class"])
        @test count(o -> o["in_SD"], os) == min(EX_BLOCK["s_D"], length(os))
        @test count(o -> o["in_SL"], os) == min(EX_BLOCK["s_L"], length(os))
        @test all(o -> o["in_SD"] == (o["draw_index"] <= EX_BLOCK["s_D"]), os)
    end
    cap = ex_build("capped", Dict{String, Any}("cap_D" => 12, "cap_L" => 3))
    orbs = cap.panel["orbit"]
    sub = ex_man(cap)["subsets"]
    SD, SL = [o["id"] for o in orbs if o["in_SD"]], [o["id"] for o in orbs if o["in_SL"]]
    @test length(SD) <= 12 && length(SL) <= 3 && issubset(SL, SD)
    @test sub["capped_D"] == true && !isempty(sub["removed_D"])
    @test all(o -> o["in_SD"], ex_ctrls(cap))                         # controls are never removed
    # the manifest's removed ids equal the recomputed ones (the cap procedure of E17 applied to the file's own members)
    members = [(id = o["id"], cell = 5 * (o["k_index"] - 1) + EX_CLS[o["class"]], rank = o["draw_index"], control = false) for o in ex_drawn(cap)]
    append!(members, [(id = o["id"], cell = 0, rank = 0, control = true) for o in ex_ctrls(cap)])
    ref = EX.select_subsets(members; s_D = 2, cap_D = 12, s_L = 1, cap_L = 3)
    @test sort(ref.SD) == sort(SD) && sort(ref.SL) == sort(SL)
    @test sort(ref.removed_D) == sort(sub["removed_D"]) && sort(ref.removed_L) == sort(sub["removed_L"])
end

@testset "AC22 counts and file" begin
    b = ex_build("base")
    man, orbs = ex_man(b), b.panel["orbit"]
    tot = man["totals"]
    @test tot["orbits"] == length(orbs)
    @test tot["orbits"] == sum(c["drawn"] for c in man["cells"]) + length(ex_ctrls(b)) == tot["drawn"] + tot["controls"]
    @test all(c -> c["drawn"] == length(ex_of(b, c["k_index"], c["class"])), man["cells"])
    # every required field present for every orbit
    req_all = ("id", "class", "control", "control_kind", "system", "x0", "p0", "x0_dec", "p0_dec", "run_seed", "orbit_seed", "in_SD", "in_SL")
    for o in orbs
        @test all(k -> haskey(o, k), req_all)
        if o["control_kind"] != "y5"
            @test all(k -> haskey(o, k), ("k_index", "K_text", "K"))
        end
        if !o["control"]
            @test all(k -> haskey(o, k), ("cell_seed", "draw_index"))
            o["class"] == "X" && @test haskey(o, "r") && haskey(o, "phi")
            o["class"] == "R" && o["k_index"] == 1 && @test haskey(o, "a") && haskey(o, "b")
            o["class"] == "R" && o["k_index"] > 1 && @test haskey(o, "source_orbit") && haskey(o, "source_point") && haskey(o, "phi")
        elseif o["control_kind"] != "y5"
            @test all(k -> haskey(o, k), ("period", "wind_m", "wind_j", "stability", "exponent"))
        else
            @test haskey(o, "exponent") && haskey(o, "mu")
        end
    end
    # read back: binary64 fields are Float64, decimal text rebuilds the 256-bit value it was written from
    gc = Set((c[1], c[2]) for c in EX.grid_centres(EX_BLOCK["grid_side"]))
    for o in orbs
        @test o["x0"] isa Float64 && o["p0"] isa Float64
        for (v, d) in ((o["x0"], o["x0_dec"]), (o["p0"], o["p0_dec"]))
            @test Float64(ex_bigdec(d)) == v
            @test ex_sig(d) <= 79
        end
        if !o["control"] && o["class"] in ("T", "C", "S")             # formed in binary64: the decimal is the exact value
            @test ex_bigdec(o["x0_dec"]) == BigFloat(o["x0"]) && ex_bigdec(o["p0_dec"]) == BigFloat(o["p0"])
        end
    end
    pi_ctrl = only(o for o in ex_ctrls(b) if o["control_kind"] == "fixed" && o["k_index"] == 2 && o["x0"] != 0)
    @test ex_sig(pi_ctrl["x0_dec"]) == 79 && abs(ex_bigdec(pi_ctrl["x0_dec"]) - BigFloat(pi)) <= 4 * ldexp(BigFloat(1), -254)
    # the manifest holds every input with its source, the values read back equal the preregistration block
    inp = man["inputs"]
    for (k, v) in EX_BLOCK
        @test inp[k]["value"] == v && inp[k]["source"] == "preregistration"
    end
    @test inp["seed"]["value"] == EX_SEED && inp["seed"]["source"] == "default"
    @test inp["bigfloat_bits"]["value"] == 256 && inp["tau"]["value"] == 1
    @test man["K_binary64"] == [ex_f64(s) for s in EX_BLOCK["K"]]
    @test man["prereg"]["header_sha256"] == man["prereg"]["computed_sha256"] == ex_header_hash(read(b.prereg, String))
    @test man["prereg"]["amendments"] == 2
    @test all(k -> haskey(man, k), ("format_version", "package_version", "git_sha", "git_dirty", "stream_rule", "draw_order", "cells",
                                    "periodic", "newton", "min_separation", "subsets", "totals", "controls"))
    @test man["min_separation"] >= 10 * EX_DUP
end

@testset "AC23 seed fields and recorded reasons" begin
    b = ex_build("base")
    @test all(o -> haskey(o, "run_seed") && haskey(o, "orbit_seed"), b.panel["orbit"])
    for c in ex_cells(b)
        @test c["drawn"] > 0 || c["reason"] != ""                      # every class of the preregistration has an orbit or a reason
    end
end

@testset "AC24 portrait content" begin
    b = ex_build("base")
    P = EX.portrait_content(b.panel; n_plot = 20)
    @test length(P) == 3
    @test P == EX.portrait_content(b.path; n_plot = 20)                # from the dictionary and from the file: one description
    man, orbs = ex_man(b), b.panel["orbit"]
    sym = Dict("hyperbolic" => :cross, "elliptic" => :circle, "parabolic" => :diamond)
    for i in 1:3
        p = P[i]
        @test p.k_index == i && p.K == man["K_binary64"][i]
        here = [o for o in orbs if get(o, "k_index", 0) == i]
        @test sort([s.id for s in p.starts]) == sort([o["id"] for o in here])          # one start marker per orbit of that K
        for s in p.starts
            o = orbs[s.id]
            @test s.class == o["class"] && ex_tdist((s.x, s.p), (o["x0"], o["p0"])) <= 1e-12
            @test 0 <= s.x < EX_TWO_PI && 0 <= s.p < EX_TWO_PI
        end
        @test sort([t.id for t in p.iterates]) == sort([o["id"] for o in here])
        for t in p.iterates
            o = orbs[t.id]
            @test length(t.xy) == 20 && t.class == o["class"]
            @test all(z -> 0 <= z[1] < EX_TWO_PI && 0 <= z[2] < EX_TWO_PI, t.xy)
            @test ex_tdist(t.xy[1], (o["x0"], o["p0"])) <= 1e-12                          # the first plotted point is the start
            @test ex_tdist(t.xy[2], ex_T(o["K"], (o["x0"], o["p0"]))) <= 1e-9
        end
        per = [q for q in man["periodic"] if q["k_index"] == i]
        @test length(p.periodic) == sum(q["q"] for q in per)                              # one marker per periodic orbit point
        @test sort([m.symbol for m in p.periodic]) == sort(vcat([fill(sym[q["class"]], q["q"]) for q in per]...))
        @test all(m -> 0 <= m.x < EX_TWO_PI && 0 <= m.p < EX_TWO_PI, p.periodic)
        @test all(c -> any(n -> occursin(c["reason"], n), p.notes),                       # empty and short cells are noted with their reason
                  (c for c in ex_cells(b) if c["k_index"] == i && c["status"] in ("empty", "short")))
        # no data derived from screening
        names_ = Set(String.(vcat(collect(propertynames(p)), [collect(propertynames(s)) for s in p.starts]...,
                                  [collect(propertynames(t)) for t in p.iterates]...)))
        @test isdisjoint(names_, Set(["lambda", "lambda_scr", "blocks", "block_exponents", "qualifying"]))
    end
end

@testset "AC29 provenance" begin
    mktempdir() do d
        git(args...) = read(`git -C $d -c user.name=ex -c user.email=ex@example.org -c commit.gpgsign=false $args`, String)
        git("init", "-q")
        write(joinpath(d, "a.txt"), "x")
        git("add", "a.txt")
        git("commit", "-q", "-m", "m")
        head = strip(git("rev-parse", "HEAD"))
        s = EX.git_state(d)
        @test s.sha == head && s.dirty == false
        write(joinpath(d, "a.txt"), "changed")
        @test EX.git_state(d).dirty == true && EX.git_state(d).sha == head
    end
    b = ex_build("base")
    root = normpath(joinpath(@__DIR__, "..", ".."))
    @test ex_man(b)["git_sha"] == strip(read(`git -C $root rev-parse HEAD`, String))
    tracked_dirty = !isempty(strip(read(`git -C $root status --porcelain --untracked-files=no`, String)))
    @test ex_man(b)["git_dirty"] == tracked_dirty
    proj = TOML.parsefile(joinpath(root, "Project.toml"))["version"]
    @test ex_man(b)["package_version"] == proj == string(EX.package_version())
end

@testset "D51 build at K = 2 pi (binary64)" begin
    # binary64 2 pi lies 2.4e-16 below the accelerator-mode saddle-node: residual-accepted ghosts near (pi/2, 0) exist
    b = ex_build("sn2pi", Dict{String, Any}("K" => ["2pi"], "q_max" => 1, "grid_side" => 16))   # D51c proxy: completes
    man = ex_man(b)
    ill = [r for r in man["ill_conditioned"] if r["k_index"] == 1 && r["q"] == 1]
    @test !isempty(ill) && ill[1]["count"] > 0 && 0 <= ill[1]["min_abs_2_minus_t"] < 1e-5
    @test man["newton"]["ill_conditioned"] >= ill[1]["count"]
    # no accepted period-1 point within 1e-5 of (pi/2, 0) or its mirror (3 pi/2, 0)
    for q in man["periodic"], pt in q["points"]
        z = setprecision(() -> (Float64(ex_bigdec(pt[1])), Float64(ex_bigdec(pt[2]))), BigFloat, 256)
        q["q"] == 1 && @test min(ex_tdist(z, (pi / 2, 0.0)), ex_tdist(z, (3pi / 2, 0.0))) > 1e-5
    end
    @test man["min_separation"] >= 10 * EX_DUP
    # D51d: the D50 outputs are gone
    @test !haskey(man, "merged") && !haskey(man["inputs"], "sn_trace")
    @test all(q -> q["class"] in ("elliptic", "parabolic", "hyperbolic"), man["periodic"])
    @test all(c -> c in ("elliptic", "parabolic", "hyperbolic"), (o["stability"] for o in ex_ctrls(b) if haskey(o, "stability")))
    @test !isdefined(Chaoslab.Experiments, :SN_CLASS) && !isdefined(Chaoslab.Experiments, :MergedPair)
    base = ex_man(ex_build("base"))
    @test haskey(base, "ill_conditioned") && all(r -> r["count"] > 0, base["ill_conditioned"])
end
