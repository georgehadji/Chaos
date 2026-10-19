# D53b plan, D53c hashes, D53d recording grid, D53g chunks and memory (planning part), D53h backend. Spec: Amendment D53, R1 to R4.
# The plan is a pure function of the block, the panel bytes, the defaults and the code version; nothing is computed here.

ex_label(p) = replace(p, r"\(\d+\)" => "")
ex_item_key(it) = (it.orbit, it.method, ex_label(it.precision), it.horizon, it.tau, it.eps === nothing ? nothing : parse(Float64, it.eps),
                   it.direction, it.coordinate_mode)
const EX_PREC_RANK = Dict("Float64" => 1, "Double64" => 2, "BigFloat" => 3)
function ex_cell_map(orbits, run)
    m = Dict{Int, String}(o["id"] => ex_ref_cell(o) for o in orbits)
    nP = maximum(keys(m))
    for i in 1:length(get(run, "h4", [])); m[nP + i] = "h4"; end
    return m
end
ex_runtype(it) = (it.method, ex_label(it.precision), it.horizon, it.tau, it.eps, it.direction, it.coordinate_mode)

# synthetic panel at the sizes of the study caps (1100 orbits with 13 controls, S_D 200, S_L 100) for the plan-only proxy of D53b [full]
function ex_panel_synth(; nP = 1100, nSD = 200, nSL = 100, nctrl = 13)
    rng = StableRNG(ex_s53('b', 7))
    orbs = Dict{String, Any}[]
    classes = ["T", "R", "C", "S", "X"]
    nd = nP - nctrl
    for id in 1:nP
        ctrl = id > nd
        d = Dict{String, Any}("id" => id, "class" => ctrl ? "control" : classes[mod1(id, 5)], "control" => ctrl,
            "control_kind" => ctrl ? "fixed" : "none", "k_index" => 1 + (id - 1) % 6, "K_text" => "1.5", "K" => 1.5, "system" => "standard",
            "x0" => 6 * rand(rng), "p0" => 6 * rand(rng), "run_seed" => string(EX_RUN_SEED), "orbit_seed" => string(EX.orbit_seed(EX_RUN_SEED, id)),
            "in_SD" => ctrl || id <= nSD - nctrl, "in_SL" => !ctrl && id <= nSL, "draw_index" => 1)
        d["x0_dec"] = string(d["x0"]); d["p0_dec"] = string(d["p0"])
        ctrl && merge!(d, Dict{String, Any}("stability" => "hyperbolic", "exponent" => 0.5, "period" => 1))
        push!(orbs, d)
    end
    return orbs
end

@testset "D53b plan on fixture $kind" for kind in (:small, :alt, :full)
    fx = ex_fix(kind)
    plan = ex_plan(fx)
    ref = ex_ref_items(fx.run, fx.orbits)
    M = length(plan.items)
    @test M == plan.counts.items == ref.nunique
    kind === :small && @test M == 52                                  # hand count: 15 x 2 + 9 x 2 + 2 + 2
    kind === :alt && @test M == 66                                    # hand count: 15 x 3 + 1 + 2 + 9 x 2
    @test [it.id for it in plan.items] == collect(1:M)                # ids 1..M, no gap
    # the set of items and their arms equal the Cartesian products of R1 item 2, merged by content
    @test Dict(ex_item_key(it) => it.arms for it in plan.items) == ref.items
    # order of R1 item 4: precision first, then the first arm; orbit order inside a run type of one first arm
    armidx = Dict(a["name"] => i for (i, a) in enumerate(fx.run["arm"]))
    keyf(it) = (EX_PREC_RANK[ex_label(it.precision)], armidx[first(it.arms)])
    @test issorted(plan.items; by = keyf)
    groups = Dict{Any, Vector{Int}}()
    for it in plan.items
        push!(get!(groups, (keyf(it), ex_runtype(it)), Int[]), it.orbit)
    end
    @test all(issorted, values(groups))
    # planned steps per arm before the merge equal the sum of products; the merge reports its saving
    for (name, b) in ref.before
        s = plan.steps.arms[name]
        @test s.items_before == b.items && s.steps_before == b.steps
        @test s.steps_after == sum(k[4] for (k, arms) in ref.items if first(arms) == name; init = 0)
    end
    @test plan.steps.total_before == ref.steps_before && plan.steps.total_after == ref.steps_after
    @test plan.steps.saving == ref.steps_before - ref.steps_after
    for label in ("Float64", "Double64", "BigFloat")
        want_b = sum(b.steps for (a, b) in ref.before if label in only(arm for arm in fx.run["arm"] if arm["name"] == a)["precisions"]; init = 0)
        haskey(plan.steps.precisions, label) || (@test want_b == 0; continue)
        @test plan.steps.precisions[label].steps_before == want_b
    end
    # chunks: one per (run type, cell); ordinals 1..C in plan order; every item in exactly one chunk
    cell = ex_cell_map(fx.orbits, fx.run)
    @test plan.counts.chunks == length(plan.chunks) == length(unique((ex_runtype(it), cell[it.orbit]) for it in plan.items))
    @test [c.ordinal for c in plan.chunks] == collect(1:length(plan.chunks))
    @test sort(vcat((c.items for c in plan.chunks)...)) == collect(1:M)
    @test issorted([minimum(c.items) for c in plan.chunks])
    for c in plan.chunks
        its = plan.items[c.items]
        @test length(unique(ex_runtype.(its))) == 1 && all(it -> cell[it.orbit] == c.cell, its)
        @test sort(c.arms) == sort(unique(vcat((it.arms for it in its)...)))
    end
    # determinism: bytes identical on repeated planning, for 1, 4 and all threads, and for another textual order of the tables
    txt = EX.plan_text(plan)
    @test plan.hash == ex_sha256(txt)
    @test !occursin(EX_DIR, txt)                                      # no path, host or time in the plan
    for nt in (1, 4, Threads.nthreads())
        p2 = ex_plan(fx; ntasks = nt)
        @test EX.plan_text(p2) == txt && p2.hash == plan.hash
    end
    rev = ex_fix_with(fx.run, "b-rev-$kind"; rev = true)
    @test read(rev.prereg) != read(fx.prereg)
    @test EX.plan_text(ex_plan(rev)) == txt
end

@testset "D53b merged items list every arm" begin
    plan = ex_plan(ex_fix(:full))
    want(arms) = count(it -> it.arms == arms, plan.items)
    @test want(["main", "dup", "dirs", "tauw"]) == 9                  # L5 at the primary values on S_D: four arms, one item
    @test want(["main", "dirs"]) == 9                                 # L1 at direction pi/4 on S_D
    @test want(["main", "epsw"]) == 9                                 # L2 at eps = the primary value
    @test want(["main", "tauw"]) == 9                                 # L4 at tau = 1
    @test plan.steps.saving > 0 && plan.steps.total_before == plan.steps.total_after + plan.steps.saving
    # arm names and arm order do not enter chunk keys or hashes; only the item lists change
    fx = ex_fix(:small)
    p1 = ex_plan(fx)
    p2 = ex_plan(ex_fix_with(ex_rb_with(fx.run, ("arm", 1, "name") => "renamed", ("arm", 3, "name") => "other"), "b-renamed"))
    @test [it.hash for it in p1.items] == [it.hash for it in p2.items]
    @test [(c.key, c.hash) for c in p1.chunks] == [(c.key, c.hash) for c in p2.chunks]
    @test p1.hash != p2.hash || EX.plan_text(p1) != EX.plan_text(p2)  # the plan text lists the arm names
end

@testset "D53b study block on a synthetic panel of the study caps (plan only)" begin
    real = ex_real_run()
    orbs = ex_panel_synth()
    pp = ex_write_panel(ex_path("b-synth-panel.toml"); orbits = orbs, accepted = real["panel_prereg_sha256_accepted"][1])
    plan = EX.plan_run(EX_REAL_PREREG, pp; code_version = EX_CODE, ntasks = 4)
    ref = ex_ref_items(real, orbs)
    @test plan.counts.items == ref.nunique
    @test plan.steps.arms["main"].steps_before + plan.steps.arms["main_l2"].steps_before == 1100 * 9 * 10^7     # P x 9 methods x N_1: the first row of the Budget table (L2 in its own arm since A5)
    @test plan.steps.arms["long"].steps_before + plan.steps.arms["long_l2"].steps_before == 200 * 9 * 10^8
    for (name, b) in ref.before
        @test plan.steps.arms[name].items_before == b.items && plan.steps.arms[name].steps_before == b.steps
    end
    @test plan.steps.total_before == ref.steps_before && plan.steps.total_after == ref.steps_after
    @test sum(p.steps_before for p in values(plan.steps.precisions)) == ref.steps_before
end

@testset "D53c hashes" begin
    plan = ex_plan(ex_fix(:alt))                                      # has L2 items (eps) and extended precision
    it5 = first(it for it in plan.items if it.method == "L3")
    it2 = first(it for it in plan.items if it.method == "L2")
    @test it2.eps !== nothing && it5.eps === nothing
    @test occursin("eps=", EX.item_text(it2)) && !occursin("eps=", EX.item_text(it5))     # absent, not defaulted
    @test all(it -> it.hash == EX.item_hash(it), plan.items)
    sample = plan.items[round.(Int, range(1, length(plan.items); length = 25))]
    @test [it.hash for it in sample] == ex_sha256_many(EX.item_text.(sample))              # O8, one platform process
    # the text carries every hashed field as name=value on its own line
    for f in (:orbit, :system, :param, :x0, :p0, :run_seed, :orbit_seed, :method, :precision, :horizon, :tau, :direction,
              :coordinate_mode, :vectors, :grid_sha256, :code_version, :hash_version)
        @test occursin("$f=$(getfield(it5, f))\n", EX.item_text(it5))
    end
    # one case per field: changing it changes the hash (18 fields)
    h0, h2 = EX.item_hash(it5), EX.item_hash(it2)
    chg = [(:orbit, it5.orbit + 1), (:system, "linear"), (:param, "9.75"), (:x0, "0.5"), (:p0, "0.25"), (:run_seed, "1"),
           (:orbit_seed, "1"), (:method, "L1"), (:precision, "Double64"), (:horizon, it5.horizon + it5.tau), (:tau, it5.tau + 1),
           (:direction, "0.4"), (:coordinate_mode, "wrapped"), (:vectors, 3), (:grid_sha256, "0"^64), (:code_version, "other"),
           (:hash_version, it5.hash_version + 1)]
    @test length(chg) + 1 >= 14
    for (f, v) in chg
        @test EX.item_hash(merge(it5, NamedTuple{(f,)}((v,)))) != h0
    end
    @test EX.item_hash(merge(it2, (eps = "1.0e-9",))) != h2
    # excluded from the hash: backend request, thread count, arm names (item level), time, launch capacity
    base = ex_fix(:alt)
    other = ex_plan(ex_fix_with(ex_rb_with(base.run, ("backend", "Float64") => "auto", ("k19", "threads") => 8,
                                            ("arm", 1, "name") => "renamed"), "c-excl"))
    @test [it.hash for it in other.items] == [it.hash for it in plan.items]
    @test [it.hash for it in ex_plan(base; ntasks = 1).items] == [it.hash for it in plan.items]
    # the code version enters every item; the chunk hash changes with the values that decide its status
    p_code = EX.plan_run(base.prereg, base.panel; code_version = "other-code")
    @test all(i -> p_code.items[i].hash != plan.items[i].hash, eachindex(plan.items))
    for (path, v) in ((("health", "fail_fraction"), 0.2), (("health", "control_bound"), 0.3), (("health", "drift_bound"), 1e4),
                      (("health", "setup_check_tol"), 1e-8))
        p = ex_plan(ex_fix_with(ex_rb_with(base.run, path => v), "c-chunk-$(path[2])"))
        @test [it.hash for it in p.items] == [it.hash for it in plan.items]
        @test all(i -> p.chunks[i].hash != plan.chunks[i].hash && p.chunks[i].key == plan.chunks[i].key, eachindex(plan.chunks))
    end
    # chunk hash: independent implementation of the text (O8); the text names the key and lists the sorted item hashes
    for c in plan.chunks
        t = EX.chunk_text(plan, c.ordinal)
        @test occursin(c.key, t)
        listed = [m.match for m in eachmatch(r"\b[0-9a-f]{64}\b", t) if m.match in Set(plan.items[i].hash for i in c.items)]
        @test listed == sort([plan.items[i].hash for i in c.items])
    end
    @test [c.hash for c in plan.chunks] == ex_sha256_many([EX.chunk_text(plan, c.ordinal) for c in plan.chunks])
    @test EX.plan_text(plan) isa String && plan.hash == ex_sha256(EX.plan_text(plan))
end

@testset "D53c independent hash on 1000 seeded items" begin
    rng = StableRNG(ex_s53('c', 1))
    function rand_item()
        method = rand(rng, EX_ALL9)
        return (; orbit = rand(rng, 1:1100), system = rand(rng, ["standard", "linear"]), param = string(round(10 * rand(rng); digits = 3)),
                x0 = string(6 * rand(rng)), p0 = string(6 * rand(rng)), run_seed = "20261012", orbit_seed = string(rand(rng, UInt64)),
                method = method, precision = rand(rng, ["Float64", "Double64", "BigFloat(256)"]),
                horizon = rand(rng, [10^3, 10^5, 10^7, 10^8]), tau = rand(rng, [1, 2, 4, 8, 16]),
                eps = method == "L2" ? string(rand(rng, [1e-6, 1e-8, 1e-10])) : nothing,
                direction = rand(rng, ["pi/4", "0.3", "0.3 + pi/5"]), coordinate_mode = rand(rng, ["lifted", "wrapped"]),
                vectors = ex_vectors(method), grid_sha256 = bytes2hex(rand(rng, UInt8, 32)), code_version = EX_CODE, hash_version = 1)
    end
    items = [rand_item() for _ in 1:1000]
    texts = [EX.item_text(it) for it in items]
    @test length(unique(texts)) >= 990
    @test [EX.item_hash(it) for it in items] == ex_sha256_many(texts)
end

@testset "D53c code version" begin
    root = mktempdir()
    mkpath(joinpath(root, "src", "Sub")); mkpath(joinpath(root, "docs")); mkpath(joinpath(root, "test"))
    write(joinpath(root, "src", "A.jl"), "a() = 1\n"); write(joinpath(root, "src", "Sub", "B.jl"), "b() = 2\n")
    write(joinpath(root, "Manifest.toml"), "# lock 1\n"); write(joinpath(root, "Project.toml"), "name = \"X\"\n")
    v0 = EX.code_version(root)
    @test v0 isa String && length(v0) >= 8 && EX.code_version(root) == v0
    write(joinpath(root, "README.md"), "docs\n"); write(joinpath(root, "docs", "x.md"), "docs\n"); write(joinpath(root, "test", "t.jl"), "@test true\n")
    @test EX.code_version(root) == v0                                  # a file outside the library source and the lock: unchanged
    write(joinpath(root, "src", "Sub", "B.jl"), "b() = 3\n")
    v1 = EX.code_version(root)
    @test v1 != v0                                                     # a source change
    write(joinpath(root, "src", "Sub", "B.jl"), "b() = 2\n")
    @test EX.code_version(root) == v0
    write(joinpath(root, "Manifest.toml"), "# lock 2\n")
    @test EX.code_version(root) != v0                                  # the dependency lock
end

@testset "D53d recording grid" begin
    kw = (points_per_decade = 20, first_steps = 100, extra_fractions = [0.1, 0.5])
    diffs = 0
    for N in (10^3, 10^4, 10^5, 10^6, 10^7), tau in (1, 2, 4, 8, 16)
        N % tau == 0 || continue
        g = EX.record_grid(N, tau; kw...)
        ref = ex_ref_grid(N, tau; ppd = 20, first_steps = 100, fracs = [0.1, 0.5])
        diffs += g == ref ? 0 : 1
        @test issorted(g; lt = <=) && all(n -> n % tau == 0, g)         # strictly increasing, multiples of tau
        @test all(s -> tau * cld(s, tau) in g, 1:100)
        @test N in g && tau * cld(N ÷ 10, tau) in g && tau * cld(N ÷ 2, tau) in g
    end
    @test diffs == 0
    @test length(EX.record_grid(10^7, 1; kw...)) == 201                # derived in R2 item 1
    @test length(EX.record_grid(10^8, 1; kw...)) == 221
    hs = [EX.grid_hash(EX.record_grid(10^5, t; kw...)) for t in (1, 2, 4, 8, 16)]
    @test length(unique(hs)) == 5 && all(h -> occursin(r"^[0-9a-f]{64}$", h), hs)
    @test EX.grid_hash(EX.record_grid(10^5, 4; kw...)) == hs[3]
end

@testset "D53g chunk ids, keys and launch capacity" begin
    # capacity = floor(f free / (M + reducer + series bytes)), series bytes = 8 columns grid steps
    cap(free; f = 0.5) = EX.launch_capacity(free; f = f, M = 240, reducer = 96, columns = 12, grid_steps = 57)
    per = 240 + 96 + 8 * 12 * 57
    for free in (2^30, 3 * 2^28 + 1234, 10^9 + 7)
        @test cap(free) == floor(Int, Rational{BigInt}(1 // 2) * free / per)
    end
    @test cap(2^31) == cap(2^30) * 2 || cap(2^31) == 2 * cap(2^30) + 1   # doubling the memory at most doubles (floor)
    for free in (2^30, 2^30 + 2 * 12345)
        @test cap(free ÷ 2) == cap(free) ÷ 2                            # halves (floor) when free memory halves
    end
    @test cap(2^30; f = 0.25) == floor(Int, Rational{BigInt}(1 // 4) * 2^30 / per)
    # groups: consecutive chunks of one run type up to the capacity; ids and keys do not depend on it
    plan = ex_plan(ex_fix(:small))
    keys0 = [(c.ordinal, c.key, c.hash) for c in plan.chunks]
    for free in (2^20, 2^24, 2^30)
        c = cap(free)
        groups = EX.launch_groups(plan; capacity = max(c, 1))
        @test vcat((g.chunks for g in groups)...) == [c.ordinal for c in plan.chunks]
        for g in groups
            @test length(unique(plan.chunks[o].run_type for o in g.chunks)) == 1
            @test g.orbits <= max(c, 1) || (length(g.chunks) == 1 && g.launches == cld(g.orbits, max(c, 1)))
        end
    end
    for capacity in (1, 7, 10^6)
        groups = EX.launch_groups(plan; capacity = capacity)
        @test all(g -> g.launches >= 1, groups) && sum(g.orbits for g in groups) == length(plan.items)
    end
    @test keys0 == [(c.ordinal, c.key, c.hash) for c in ex_plan(ex_fix(:small)).chunks]
end

@testset "D53h backend resolution (R4)" begin
    devs = Dict("no device" => ExDev(false, false, true), "device not functional" => ExDev(true, false, true),
                "kernel does not compile" => ExDev(true, true, false), "device" => ExDev(true, true, true))
    recs = Dict("no record" => nothing, "record, threads equal, faster" => (threads = 4, gpu_faster = true),
                "record, threads equal, not faster" => (threads = 4, gpu_faster = false),
                "record at another thread count" => (threads = 8, gpu_faster = true))
    bad = String[]
    reasons = Set{String}()
    for T in (Float64, Double64, BigFloat), req in (:cpu, :gpu, :auto), (dn, dev) in devs, (rn, rec) in recs, method in ("L5", "L6")
        method == "L6" && T !== Float64 && continue
        got = ex_with_device(() -> EX.resolve_run_backend(req, T, method; threads = 4, record = rec), dev)
        want = ex_expect_backend(req, T, method == "L6", dev, rec, 4)
        push!(reasons, got.reason)
        (got.backend, got.reason) == want || push!(bad, "$T $req $dn $rn $method: got $((got.backend, got.reason)), want $want")
    end
    @test isempty(bad)
    @test reasons <= Set(vcat([""], Chaoslab.FALLBACK_REASONS))        # no tenth reason
    @test EX.cpu_only_method("L6") && !EX.cpu_only_method("L5") && !EX.cpu_only_method("L1")
    # a record taken at another thread count resolves to CPU with the existing reason and is noted
    r = ex_with_device(() -> EX.resolve_run_backend(:auto, Float64, "L5"; threads = 4, record = (threads = 8, gpu_faster = true)),
                       ExDev(true, true, true))
    @test r.backend === :cpu && r.reason == Chaoslab.FALLBACK_REASONS[7] && occursin("thread", r.note)
    # the table of a plan: Float64 auto may use the device, Double64 and BigFloat resolve to CPU under their block requests
    plan = ex_plan(ex_fix(:full))
    tbl = ex_with_device(() -> EX.backend_table(plan; threads = 4, record = (threads = 4, gpu_faster = true)), ExDev(true, true, true))
    @test tbl[("Float64", "L5")].backend === :gpu && tbl[("Float64", "L5")].request === :auto && tbl[("Float64", "L5")].reason == ""
    @test tbl[("Float64", "L6")].backend === :cpu && tbl[("Float64", "L6")].reason == "reducer is CPU-only"
    @test tbl[("Double64", "L5")].backend === :cpu && tbl[("Double64", "L5")].reason == ""
    @test tbl[("BigFloat", "L5")].backend === :cpu && tbl[("BigFloat", "L5")].reason == ""
    @test Set(keys(tbl)) == Set((ex_label(it.precision), it.method) for it in plan.items)
    tbl0 = ex_with_device(() -> EX.backend_table(plan; threads = 4, record = nothing), ExDev(true, true, true))
    @test tbl0[("Float64", "L5")].backend === :cpu && tbl0[("Float64", "L5")].reason == "GPU not faster for the element type"
end
