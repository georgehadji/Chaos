# Amendment D86, criteria D86i (stretching), D86j (frequency feed), D86k (frequency values and flags), D86l (phase advance), D86m (weighted Birkhoff).
# Spec: experiments-spec.md, "Amendment D86": E21, E25, E26, E28, E29, R14, I22, I23, O16 to O18. API assumed: experiments-test-plan.md, "API assumed by the D86 tests".

# ---------------------------------------------------------------- blocks
function d86_win_block()
    specs = [d86_spec("frequency"; id = "fq", windows = [16, 24, 40]), d86_spec("frequency"; id = "fqp", observable = "q_plus_i_p", windows = [16, 24, 40]),
             d86_spec("phase_advance"; id = "pa", windows = [16, 24, 32]), d86_spec("phase_advance"; id = "pal", observable = "lifted", windows = [16, 24, 32]),
             d86_spec("birkhoff"; id = "bk", windows = [16, 24, 40])]
    arms = [d86_arm("w1", "S_D", ["fq", "fqp", "bk"], ["Float64"], [80]; coordinate_mode = ["wrapped"]),
            d86_arm("w1l", "S_D", ["fqp"], ["Float64"], [80]; coordinate_mode = ["lifted"]),
            d86_arm("w2", "S_D", ["pa", "pal"], ["Float64"], [65]; coordinate_mode = ["wrapped"]),
            d86_arm("w2l", "S_D", ["pal"], ["Float64"], [65]; coordinate_mode = ["lifted"]),
            d86_arm("w3", "h4", ["fq", "fqp", "bk"], ["Double64"], [80]; coordinate_mode = ["wrapped"]),
            d86_arm("w3p", "h4", ["pa"], ["Double64"], [65]; coordinate_mode = ["wrapped"]),
            d86_arm("l1a", "S_D", String[], ["Float64"], [80]; methods = ["L1"], coordinate_mode = ["wrapped", "lifted"]),   # a catalog method: N / 10 = 8 must be a multiple of tau
            d86_arm("l1c", "h4", String[], ["Double64"], [80]; methods = ["L1"], coordinate_mode = ["wrapped"])]
    for a in arms
        isempty(a["indicators"]) && delete!(a, "indicators")
    end
    return d86_block(arms, specs)
end
d86_win_run() = d86_run("win", d86_win_block())

# stretching: all orbits of S_D with the counts of S_L, tau 1 and 4, plus the same spec on orbit 1 alone
d86_st_spec() = d86_spec("stretching"; id = "s1")
function d86_st_block(; orbits = "S_D", bins = 20, range = "derived", direction = nothing)
    sp = d86_spec("stretching"; id = "s1", bins = bins, range = range)
    kw = direction === nothing ? (;) : (; direction = direction)
    return d86_block([d86_arm("st", orbits, ["s1"], ["Float64"], [100]; tau = [1, 4], kw...), d86_arm("stw", orbits, ["s1"], ["Float64"], [100]; coordinate_mode = ["wrapped"], kw...)], [sp])
end

# ================================================================ D86i
@testset "D86i stretching: series, scalars, counts, pooled counts, ranges, closed forms" begin
    run = d86_run("st", d86_st_block())
    plan = run.plan
    @test run.res.status == "complete"
    spec = d86_st_spec()
    rng = StableRNG(d86_seed('i', 1))
    bad = Dict{String, Int}()
    miss!(tag, ok) = ok || (bad[tag] = get(bad, tag, 0) + 1)
    pooled = Dict{Int, Matrix{Int}}()
    for file in ex_chunk_files(run.store)
        c = plan.chunks[ex_ordinal(file)]
        raw = d86_raw(file)
        want = zeros(Int, spec["bins"] + 2, 2)
        for i in c.items
            it = plan.items[i]; o = it.orbit
            t = raw["orbits/$o/t"]
            ser = d86_ds(raw, o, "s1", "series")
            miss!("columns", d86_at(raw, o, "s1", "series", "columns") == ["stretching.D"] && size(ser) == (length(t), 1))
            for k in sort(unique([rand(rng, 1:length(t), 5); length(t)]))
                ref = d86_stretch_final(it, spec, t[k])
                miss!("D", reinterpret(UInt64, ser[k, 1]) == reinterpret(UInt64, Float64(ref.fin.D)))
            end
            last = d86_stretch_final(it, spec, it.horizon)
            sc = d86_ds(raw, o, "s1", "scalars")
            f = last.fin
            # E25: "the first grid step at which the handoff count reaches min_per_cell x bins"; when the run never reaches it (tau = 4: 25 handoffs against 40) the
            # spec names no value, so that entry is not compared
            ifull = findfirst(x -> x ÷ it.tau >= spec["min_per_cell"] * spec["bins"], t)
            expected = Float64[f.handoffs; f.sum; f.comp; f.sum2; f.comp2; f.mean; f.std; f.below; f.above]
            ifull === nothing || push!(expected, t[ifull])
            miss!("scalars length", length(sc) == 18)
            miss!("scalars values", all(v -> any(s -> reinterpret(UInt64, s) == reinterpret(UInt64, v), sc), expected))
            rg = d86_at(raw, o, "s1", "scalars", "range")                                   # two binary64, lo then hi, on `scalars` (correction 5)
            miss!("range", rg isa Vector{Float64} && rg == [last.range[1], last.range[2]])
            # counts exist only for the orbits of the named subset (S_L = orbit 1), sum to the handoff count per column, equal the stand-alone counts
            has = haskey(raw, "orbits/$o/indicators/s1/counts")
            miss!("counts subset", has == plan.panel.orbits[o]["in_SL"])
            if has
                cnt = d86_ds(raw, o, "s1", "counts")
                miss!("counts", cnt isa Matrix{UInt32} && size(cnt) == (spec["bins"] + 2, 2) && all(==(it.horizon ÷ it.tau), vec(sum(Int.(cnt); dims = 1))) && cnt == f.counts)
            end
            want .+= Int.(last.fin.counts)
        end
        # the pooled counts of a chunk: the exact integer sum over the chunk's orbits of the counts at the last grid step
        pc = raw["pooled/s1/counts"]
        miss!("pooled", pc isa Matrix{Int64} && pc == want && raw["pooled/s1@orbits"] == length(c.items))
        pooled[c.ordinal] = pc
    end
    @test isempty(bad)
    # unchanged when the order of the chunks is reversed
    rev = d86_run("st-rev", d86_st_block(); order = :reverse)
    revp = Dict(ex_ordinal(f) => d86_raw(f)["pooled/s1/counts"] for f in ex_chunk_files(rev.store))
    @test revp == pooled
    # an item recorded alone and in a group with other orbits has identical stored data (the workspace is per orbit)
    alone = d86_run("st-alone", d86_st_block(; orbits = "S_L"))
    ra = Dict{Tuple, Any}(); rb = Dict{Tuple, Any}()
    d86_each(alone) do it, o, raw, c, file
        ra[(it.tau, it.coordinate_mode)] = Dict(n => raw["orbits/$o/indicators/s1/$n"] for n in ("series", "scalars", "counts"))
    end
    d86_each(run) do it, o, raw, c, file
        o == 1 && (rb[(it.tau, it.coordinate_mode)] = Dict(n => raw["orbits/$o/indicators/s1/$n"] for n in ("series", "scalars", "counts")))
    end
    @test !isempty(ra) && all(k -> all(n -> d86_same(ra[k][n], rb[k][n]), ("series", "scalars", "counts")), keys(ra))
end

@testset "D86i stretching closed forms through the recorder: Y5 hyperbolic and the K = 0 shear (O18)" begin
    # Y5 hyperbolic (stretching SN8b): for a unit column v of the frame the stretching number of handoff k is z_k = ln(norm(A^k v) / norm(A^(k-1) v)) (tau = 1), computed here in
    # BigFloat(256) from the stored matrix; the counts are those of the cells of z_k "tolerance 0 under the edge condition" and D = ||N1 - N2|| / (k sqrt(h)) (SN5) within 4 eps.
    # An exact eigenframe (direction pi/4) cannot be used: the column along the contracting direction is rounded off it, and the rounding grows like mu^k, so after about 20
    # handoffs the recorded column is the expanding direction (the stand-alone reducer does the same); the closed form is therefore taken for the direction 0.3, where it
    # holds for all k. The derived range is the exact bound padded by 8 eps, so z of Y5 sits at its edge by construction; a given range [-2, 2] puts the z away from the edges.
    run = d86_run("st-cf", d86_st_block(; bins = 21, range = [-2.0, 2.0], direction = ["0.3"]))
    plan = run.plan
    bad = String[]
    nchecked = Ref(0)
    d86_each(run) do it, o, raw, c, file
        (o == 15 && it.tau == 1 && it.coordinate_mode == "lifted") || return
        nchecked[] += 1
        lo, hi = d86_at(raw, o, "s1", "scalars", "range")
        (lo == -2.0 && hi == 2.0) || push!(bad, "range")
        D = d86_ds(raw, o, "s1", "series")[:, 1]
        n = it.horizon
        W = EX.deviation_frame(Float64, it.direction, 2)
        A = Chaoslab.y5_control(parse(Float64, it.param)).sys.A
        counts, margin = setprecision(256) do
            B = BigFloat(21); h = (BigFloat(hi) - BigFloat(lo)) / B
            cnt = zeros(Int, 23, 2); mar = BigFloat(1)
            for j in 1:2
                v = BigFloat.(collect(W[:, j])); Ab = BigFloat.(Matrix(A))
                prev = log(sqrt(sum(abs2, v)))
                for k in 1:n
                    v = Ab * v
                    cur = log(sqrt(sum(abs2, v)))
                    z = cur - prev; prev = cur
                    u = (z - BigFloat(lo)) / h
                    cell = z < lo ? 0 : z > hi ? 22 : min(Int(floor(u)) + 1, 21)
                    cnt[cell + 1, j] += 1
                    mar = min(mar, abs(u - round(u)))
                end
            end
            cnt, mar
        end
        margin > 1e-9 || push!(bad, "edge condition")                                       # a z within 1e-9 of a cell edge would make the counts depend on rounding
        pc = d86_stretch_final(it, d86_spec_of(plan, it), n).fin.counts
        Int.(pc) == counts || push!(bad, "counts")
        dcf = setprecision(() -> Float64(sqrt(sum((counts[:, 1] .- counts[:, 2]) .^ 2)) / (n * sqrt((BigFloat(4) / 21)))), BigFloat, 256)
        abs(D[end] - dcf) <= 4 * EX_EPS * dcf || push!(bad, "D")
    end
    @test nchecked[] == 1
    @test isempty(bad)
    # K = 0 shear, identity frame: column 1 = (1, 0) has b = 0, so every stretching number is exactly 0 and all its counts sit in the central cell (bins odd)
    orbs = ex_panel_orbits(:std)
    for o in orbs
        o["control"] || (o["K_text"] = "0"; o["K"] = 0.0)
    end
    blk = d86_one(d86_spec("stretching"; id = "s1", bins = 21); orbits = "S_L", direction = ["0"])
    pp = ex_write_panel(ex_path("d86-i-shear-panel.toml"); orbits = orbs)
    pr = ex_write_prereg_run(ex_path("d86-i-shear-prereg.md"), blk)
    store = ex_store("d86-i-shear")
    res = ex_campaign((; run = blk, panel = pp, prereg = pr, tag = "i-shear", orbits = orbs); store = store)
    @test res.status == "complete"
    raw = d86_raw(first(ex_chunk_files(store)))
    cnt = d86_ds(raw, 1, "s1", "counts"); k = raw["orbits/1/t"][end]
    @test count(!iszero, cnt[:, 1]) == 1 && cnt[12, 1] == k                                   # cell 11 sits at row 12 (row = cell + 1)
    lo, hi = d86_at(raw, 1, "s1", "scalars", "range"); h = (hi - lo) / 21
    Dstored = d86_ds(raw, 1, "s1", "series")[end, 1]
    Dalg = sqrt(sum((Int.(cnt[:, 1]) .- Int.(cnt[:, 2])) .^ 2)) / (k * sqrt(h))               # SN5 from the stored counts
    @test abs(Dstored - Dalg) <= 4 * EX_EPS * Dalg && 0 <= Dstored <= sqrt(2 / h) * (1 + 4 * EX_EPS)
end

# ================================================================ D86j
@testset "D86j frequency feed: capacity, stored bytes, resident memory, series identity, no series stored" begin
    run = d86_win_run()
    plan = run.plan
    @test run.res.status == "complete"
    fqs = [it for it in plan.items if it.method == "fq"]
    # (a) capacity = floor(f free / bytes(orbit)) exactly, for 3 free-memory values and halving; B = 32 W_max (Float64), 64 W_max (Double64)
    for it in (first(i for i in fqs if i.precision == "Float64"), first(i for i in fqs if i.precision == "Double64"))
        b = EX.orbit_bytes(plan, it)
        wm = 40
        @test b.B == (it.precision == "Float64" ? 32 : 64) * wm
        @test b.total == b.M + b.R + b.S + b.B
        for free in (4 * 10^9, 8 * 10^9, 16 * 10^9 + 122)
            cap = EX.group_capacity(plan, it; free_bytes = free, f = 0.5)
            @test cap == Int(fld(big(free), 2 * big(b.total)))                          # floor(f free / bytes(orbit)) with f = 1/2, exact integers
            @test EX.group_capacity(plan, it; free_bytes = free ÷ 2, f = 0.5) == Int(fld(big(free ÷ 2), 2 * big(b.total)))
        end
    end
    # (b) 96 stored bytes per orbit and window: eleven binary64 and one flag word
    bad = String[]
    d86_each(run; methods = ["fq", "fqp"]) do it, o, raw, c, file
        fin, flg, bits = d86_final(raw, o, it.method)
        (fin isa Matrix{Float64} && size(fin) == (3, 11) && flg isa Vector{Int64} && length(flg) == 3 && 8 * size(fin, 2) + 8 == 96) || push!(bad, "layout $(it.method) $o")
        d86_at(raw, o, it.method, "final", "columns") == D86_FREQ_COLS || push!(bad, "columns")
    end
    @test isempty(bad)
    # (c) the buffer is released: resident bytes after a chunk exceed those before by at most the stored results plus 64 kB (full GC, allocation counters)
    # 40 orbits (the panel :big): the buffers of a chunk would be 40 x 32 x 250 = 320 kB if they were kept, five times the 64 kB of the criterion
    bigblk = d86_one(d86_spec("frequency"; id = "fq", windows = [250]); orbits = "P", mode = "wrapped", name = "big")
    fb = d86_fx("j-big", bigblk; panel = :big)
    plan2 = ex_plan(fb)
    ctx = EX.make_ctx(plan2)
    items = NamedTuple[it for it in plan2.items if it.method == "fq"]
    # warm-up outside the measurement: one full pass over the same items (the constant series of the control orbits fill a one-time cache of the transform; the criterion
    # is that the buffer is not kept, so what is measured is the growth of a second pass)
    EX.compute_records(ctx, items, 1)
    GC.gc(true); before = Base.gc_live_bytes()
    recs = EX.compute_records(ctx, items, 1)
    GC.gc(true); after = Base.gc_live_bytes()
    @test after - before <= Base.summarysize(recs) + 64 * 1024 && length(items) * 32 * 250 > 64 * 1024
    # (d) stored frequencies, delta and D equal the pure function on the independently built series (O17), both observables, three windows, both precisions
    bad = String[]
    seen = Set{Tuple}()
    d86_each(run; methods = ["fq", "fqp"]) do it, o, raw, c, file
        spec = d86_spec_of(plan, it)
        fin, flg, bits = d86_final(raw, o, it.method)
        ref = d86_ref_freq(it, spec)
        push!(seen, (it.method, it.precision))
        for (r, row) in enumerate(ref)
            d86_isequal(fin[r, :], row.vals) || push!(bad, "values $(it.method) $o w$r")
            names = d86_flagset(flg[r], bits)
            all(f -> (f in row.flags) == (f in names), D86_FLAGS_CHECKED) || push!(bad, "flags $(it.method) $o w$r")
        end
    end
    @test isempty(bad)
    @test seen == Set([("fq", "Float64"), ("fq", "Double64"), ("fqp", "Float64"), ("fqp", "Double64")])
    # (e) the orbit points of the item are bitwise those of the L1 item (the feed's orbit is the stream's)
    l1 = Dict{Tuple, Any}()
    d86_each(run; methods = ["L1"]) do it, o, raw, c, file
        it.horizon == 80 && (l1[(o, it.precision, it.coordinate_mode)] = raw)
    end
    bad = String[]
    d86_each(run; methods = ["fq", "fqp"]) do it, o, raw, c, file
        r1 = get(l1, (o, it.precision, it.coordinate_mode), nothing)
        r1 === nothing && return push!(bad, "no twin")
        for n in ("t", "state", "diagnostics/det", "diagnostics/symp")
            d86_same(raw["orbits/$o/$n"], r1["orbits/$o/$n"]) || push!(bad, "$n")
        end
    end
    @test isempty(bad)
    # (f) no dataset of the chunk file holds the series: only the grid rows, the table, the flags and the scalars
    allowed = r"^(orbits/\d+/(t|state|winding|diagnostics/\w+|cost/\w+|indicators/(fq|fqp)/(final|flags)))$"
    bad = String[]
    d86_each(run; methods = ["fq", "fqp"]) do it, o, raw, c, file
        for (k, v) in raw
            (v isa AbstractArray && startswith(k, "orbits/$o/") && !occursin("@", k)) || continue          # datasets only (attributes such as `columns` are string vectors)
            occursin(allowed, k) || push!(bad, k)
            v isa AbstractArray && length(v) >= 2 * 40 && !(k in ("orbits/$o/state", "orbits/$o/winding") || occursin(r"/(t|diagnostics/\w+)$", k)) && push!(bad, "long $k")
        end
    end
    @test isempty(bad)
end

# ================================================================ D86k
# The Y5 rotation block of naff-fma criterion 7(a) cannot be written in the panel format (correction 6d): the criterion keeps the exact fixed points and the K = 0 map;
# the Y5 origin start covers the `undefined window` (D86l) and `degenerate` bits.
function d86_k0_orbits()
    orbs = ex_panel_orbits(:std)
    rng = StableRNG(d86_seed('k', 1))
    p0s = [0.1, 1.0, 2.0, 3.3, 5.0, 6.1, Float64(pi), 1.0, 2.0, 3.3, 5.0]
    for o in orbs
        o["control"] && continue
        id = o["id"]
        o["K_text"] = "0"; o["K"] = 0.0
        q = 2pi * rand(rng)
        p = id <= 11 ? p0s[id] : 0.0
        if id == 12                                                         # the island orbit: close to the elliptic fixed point (pi, 0) at K = 0.1
            o["K_text"] = "0.1"; o["K"] = 0.1; q = Float64(pi) + 0.03
        end
        o["x0"] = q; o["p0"] = p; o["x0_dec"] = string(q); o["p0_dec"] = string(p)
    end
    return orbs
end
function d86_k0_block()
    sp = [d86_spec("frequency"; id = "f0", windows = [32, 64], mean_removal = false), d86_spec("frequency"; id = "fm", windows = [32, 64], mean_removal = true),
          d86_spec("frequency"; id = "fp", observable = "q_plus_i_p", windows = [32, 64], mean_removal = false)]
    arms = [d86_arm("k0", "P", ["f0", "fm", "fp"], ["Float64"], [128]; coordinate_mode = ["wrapped"])]
    return d86_block(arms, sp)
end
function d86_k0_run(tag, blk; kw...)
    orbs = d86_k0_orbits()
    pp = ex_write_panel(ex_path("d86-$tag-panel.toml"); orbits = orbs)
    pr = ex_write_prereg_run(ex_path("d86-$tag-prereg.md"), blk)
    fx = (; run = blk, panel = pp, prereg = pr, tag, orbits = orbs)
    store = ex_store("d86-$tag")
    res = ex_campaign(fx; store = store, kw...)
    return (; fx, store, res, plan = ex_plan(fx))
end

@testset "D86k frequency values: the K = 0 map, flags, fixed-point controls, island orbit, injected sample" begin
    run = d86_k0_run("k0", d86_k0_block())
    plan = run.plan
    @test run.res.status == "complete"
    tp = setprecision(() -> 2 * BigFloat(pi), BigFloat, 256)
    bad = String[]
    d86_each(run) do it, o, raw, c, file
        fin, flg, bits = d86_final(raw, o, it.method)
        names = d86_flagset.(flg, Ref(bits))
        if it.method == "f0" && o <= 11
            p0 = plan.panel.orbits[o]["p0"]
            nuref = setprecision(() -> (x = BigFloat(p0) / tp; x - round(x)), BigFloat, 256)
            for r in 1:2, col in (1, 2)
                setprecision(() -> (d = abs(BigFloat(fin[r, col]) - nuref); min(d, 1 - d)), BigFloat, 256) <= 16 * EX_EPS || push!(bad, "nu orbit $o window $r")
            end
            "below_floor" in names[1] || push!(bad, "below_floor orbit $o")
        end
        # The hyperbolic fixed point (0, 0) is exact for the binary64 map: exp(i q) is the constant 1. naff-fma 11(b) (replacement of item 8): a window of equal samples is
        # degenerate "when either mean removal is on or the common value is zero"; with mean removal off and a nonzero common value it "is a line at frequency 0: status ok"
        # with the near-zero flag. So `fm` (mean removal on) gives the `degenerate` bit, and `f0` (mean removal off) a zero frequency with `near_zero1`.
        if it.method == "fm" && o == 14
            "degenerate" in names[1] || push!(bad, "degenerate orbit $o")
        end
        if it.method == "f0" && o == 14
            (fin[1, 1] == 0.0 && "near_zero1" in names[1] && !("degenerate" in names[1])) || push!(bad, "zero frequency orbit $o")
        end
        # whatever the orbit: the table and the flags are those of the pure function
        spec = d86_spec_of(plan, it)
        ref = d86_ref_freq(it, spec)
        for (r, row) in enumerate(ref)
            d86_isequal(fin[r, :], row.vals) || push!(bad, "pure function $(it.method) orbit $o")
            all(f -> (f in row.flags) == (f in names[r]), D86_FLAGS_CHECKED) || push!(bad, "flags $(it.method) orbit $o")
        end
        # the island orbit (12) with the mean removed: the pure function sets mirror and near_zero, and so do the stored bits
        if it.method == "fm" && o == 12
            any(r -> "mirror1" in r.flags || "mirror2" in r.flags, ref) || push!(bad, "fixture: island orbit shows no mirror")
            any(r -> "near_zero1" in r.flags, ref) || push!(bad, "fixture: island orbit shows no near_zero")
        end
    end
    @test isempty(bad)
    # the fixed-point controls count in the manifest, the chunk is ok
    man = EX.read_json(run.res.manifest)
    ch = only(c for c in man["chunks"] if c["status"] == "ok" && get(c, "key", "") != "" && occursin("method=fm", c["key"]) && occursin("cell=controls", c["key"]))
    @test ch["indicator_flags"]["degenerate"] >= 1
    # an injected non-finite sample fails the chunk and names the orbit id and the sample index (naff-fma criterion 11(a))
    blk = d86_one(d86_spec("frequency"; id = "fq"); orbits = "S_D", mode = "wrapped")
    fx = d86_fx("k-nan", blk)
    store = ex_store("d86-k-nan")
    inject(info, orbit, buffer) = orbit == 2 ? (buffer[7] = complex(NaN, 0.0)) : nothing
    ex_campaign(fx; store = store, hooks = EX.RunHooks(clock = ex_clock(info -> 0.0), series = inject))
    chn = EX.read_chunk(first(f for f in ex_chunk_files(store) if EX.read_chunk(f).attrs["status"] == "failed"))
    @test any(r -> occursin("orbit 2", r) && occursin("sample 7", r), chn.attrs["failure_reasons"])
end

# ================================================================ D86l
@testset "D86l phase advance: the table, the answers on the K = 0 map, flags, plan-time refusals" begin
    run = d86_win_run()
    plan = run.plan
    bad = String[]
    d86_each(run; methods = ["pa", "pal"]) do it, o, raw, c, file
        spec = d86_spec_of(plan, it)
        fin, flg, bits = d86_final(raw, o, it.method)
        names = d86_flagset.(flg, Ref(bits))
        Wm = maximum(spec["windows"])
        # the stand-alone reducer on the same orbit (O16), chain of one: planar, and lifted from the state coordinate and its winding count (wrapped) or the lifted coordinate
        # (the linear control Y5 has no torus: its state has two components in both modes and no winding count, E26)
        obs = spec["observable"] == "planar" ? Chaoslab.PlanarObservable(1, 2) :
              (it.coordinate_mode == "wrapped" && it.system != "linear") ? Chaoslab.LiftedAngleObservable(1, 3) : Chaoslab.LiftedAngleObservable(1)
        fr = ex_with_prec(it.precision) do T
            sys, x = ex_sys_x(it, T)
            red = Chaoslab.WeightedPhaseAdvance(window = Tuple(spec["windows"]), observable = obs, birkhoff_order = spec["birkhoff_order"])
            st = ex_Tg.propagate_orbit(sys, x, SMatrix{2, 1, T}(one(T), zero(T)), it.tau * (2Wm + 1), it.tau, (red,))[3][1]
            ex_Tg.reducer_final(red, st)
        end
        for r in eachindex(spec["windows"])
            d86_isequal(Float64[fr.nu1[r], fr.nu2[r], fr.delta[r], fr.D[r]], fin[r, :]) || push!(bad, "reducer $(it.method) orbit $o w$r")
            (("below_floor" in names[r]) == fr.below_floor[r]) || push!(bad, "below_floor $(it.method) orbit $o")
        end
        if it.method == "pa"                                                    # planar: also the pure function on the independently built series (O17)
            for (r, row) in enumerate(d86_ref_pa(it, spec))
                d86_isequal(fin[r, :], row.vals) || push!(bad, "pure function orbit $o w$r")
            end
        end
        sc = d86_ds(raw, o, it.method, "scalars"); cn = d86_at(raw, o, it.method, "scalars", "columns")
        (cn == ["handoffs", "near_pi", "complete"] && sc[1] == 2Wm + 1 && sc[2] == fr.near_pi && sc[3] == 1.0) || push!(bad, "scalars $o")
        d86_pinfo(plan, o) !== nothing && d86_pinfo(plan, o)["control_kind"] == "y5" && it.method == "pa" && !("undefined window" in names[1]) && push!(bad, "undefined window Y5")     # the Y5 start is the origin: a planar sample at the origin
    end
    @test isempty(bad)
    # The lifted observable on the K = 0 map gives p_0 / (2 pi). naff-fma 7(b) gates 16 eps for "wrapped state with winding counts ... lifted post-kick momentum rebuilt from the
    # winding counts": the wrapped mode. In the lifted mode the coordinate itself is rounded at each step (at most |q| eps / 2 per increment; 7(c): "the series built from the
    # lifted coordinate is reported with its error (about |q| eps), not gated"), so the lifted mode is gated by 16 eps plus |q|max eps, the scale 7(c) names.
    blk = d86_block([d86_arm("k0", "P", ["pl"], ["Float64"], [65]; coordinate_mode = ["wrapped", "lifted"])], [d86_spec("phase_advance"; id = "pl", observable = "lifted", windows = [16, 32])])
    k0 = d86_k0_run("l-k0", blk)
    @test k0.res.status == "complete"
    tp = setprecision(() -> 2 * BigFloat(pi), BigFloat, 256)
    bad = String[]
    d86_each(k0) do it, o, raw, c, file
        (o <= 11 && !k0.plan.panel.orbits[o]["control"]) || return
        fin, _, _ = d86_final(raw, o, "pl")
        p0 = k0.plan.panel.orbits[o]["p0"]
        qmax = maximum(abs, raw["orbits/$o/state"][:, 1])
        bound = it.coordinate_mode == "wrapped" ? 16 * EX_EPS : (16 + qmax) * EX_EPS
        for r in 1:2, col in 1:2
            abs(BigFloat(fin[r, col]) - BigFloat(p0) / tp) <= bound || push!(bad, "nu orbit $o $(it.coordinate_mode)")
        end
    end
    @test isempty(bad)
    # an observable index beyond the dimension and a run of another length are refused before anything is written (see D86f(c))
end

# ================================================================ D86m
@testset "D86m weighted Birkhoff: 31 and 32 observables, the table, the digit bounds, constants, zero, the K = 0 map" begin
    run = d86_win_run()
    plan = run.plan
    spec = d86_spec("birkhoff"; windows = [16, 24, 40])
    bad = String[]
    d86_each(run; methods = ["bk"]) do it, o, raw, c, file
        s = d86_spec_of(plan, it)
        fin, flg, bits = d86_final(raw, o, "bk")
        ref = d86_ref_bk(it, s)
        K = length(s["observables"])
        size(fin) == (3K, 8) || push!(bad, "shape")
        for (r, row) in enumerate(ref.rows)
            d86_isequal(fin[r, :], row) || push!(bad, "values orbit $o row $r")
            names = d86_flagset(flg[r], bits)
            (("below_floor" in names) == ref.below[r]) && (("zero" in names) == ref.zero[r]) || push!(bad, "flags orbit $o row $r")
            (("reldig_undefined" in names) == (isnan(row[7]) && !ref.zero[r])) || push!(bad, "reldig_undefined orbit $o row $r")
        end
        # absdig never exceeds -log10(phi), reldig never -log10(phi / m), phi = 32 eps S
        for r in eachindex(ref.rows)
            avg1, avg2, S, _, _, absdig, reldig, _ = fin[r, :]
            phi = 32 * (it.precision == "Float64" ? EX_EPS : 2.0^-100) * S
            isnan(absdig) || absdig <= -log10(phi) * (1 + 1e-12) || push!(bad, "absdig bound orbit $o")
            m = (abs(avg1) + abs(avg2)) / 2
            isnan(reldig) || reldig <= -log10(phi / m) * (1 + 1e-12) || push!(bad, "reldig bound orbit $o")
        end
        if it.precision == "Float64" && it.coordinate_mode == "wrapped"
            raw["@observable_names"] == [Chaoslab.observable_name(ob) for ob in d86_observables(s, "wrapped")] || push!(bad, "names")
        end
        if o == 14                                                                                  # the hyperbolic fixed point (0, 0): cos q = 1 and sin q = 0 for every sample
            n1, n2 = d86_flagset(flg[1], bits), d86_flagset(flg[2], bits)
            (reinterpret(UInt64, fin[1, 5]) == 0x0000000000000000 && "below_floor" in n1 && !("zero" in n1)) || push!(bad, "constant observable")
            ("zero" in n2 && "below_floor" in n2) || push!(bad, "zero observable")
        end
    end
    @test isempty(bad)
    # K = 31, H = 1 at 20 seeded starts; 32 observables are refused at validation
    obs31 = vcat([Dict{String, Any}("kind" => k, "coordinate" => c, "frequency" => f) for f in 1:8 for c in 1:2 for k in ("cos", "sin")])[1:31]
    s31 = d86_spec("birkhoff"; id = "b31", observables = obs31, windows = [16])
    @test_throws ArgumentError EX.run_block(d86_one(d86_spec("birkhoff"; id = "b32", observables = vcat(obs31, [Dict{String, Any}("kind" => "lifted", "coordinate" => 1)]), windows = [16]); mode = "wrapped"))
    blk = d86_one(s31; orbits = "P", mode = "wrapped")
    r31 = d86_run("m31", blk; panel = :big)
    plan31 = r31.plan
    bad = String[]
    rng = StableRNG(d86_seed('m', 1))
    starts = sort(unique(rand(rng, 1:40, 20)))
    n31 = Ref(0)
    d86_each(r31) do it, o, raw, c, file
        o in starts || return
        n31[] += 1
        fin, flg, bits = d86_final(raw, o, "b31")
        ref = d86_ref_bk(it, s31)
        size(fin) == (31, 8) && all(r -> d86_isequal(fin[r, :], ref.rows[r]), 1:31) || push!(bad, "K=31 orbit $o")
    end
    @test isempty(bad) && n31[] == length(starts) >= 15
    # cos and sin in lifted mode are refused
    @test d86_rej(() -> EX.run_block(d86_one(d86_spec("birkhoff"); mode = "lifted")), "observables", "lifted")
    # the K = 0 map (weighted-birkhoff criterion 3, N = 1000, run length 2000; the one check that keeps the window 1000, correction 2 of the D86 stage 02 checkpoint): cos q and sin q averages within 4 eps, |Delta| <= 8 eps, lifted p averages p_0 within 2 eps p_0
    kb = d86_block([d86_arm("k0", "P", ["bk0"], ["Float64"], [2000]; coordinate_mode = ["wrapped"])],
                   [d86_spec("birkhoff"; id = "bk0", windows = [1000], birkhoff_order = 1,
                             observables = [Dict{String, Any}("kind" => "cos", "coordinate" => 1, "frequency" => 1), Dict{String, Any}("kind" => "sin", "coordinate" => 1, "frequency" => 1),
                                            Dict{String, Any}("kind" => "lifted", "coordinate" => 2)])])
    k0 = d86_k0_run("m-k0", kb)
    bad = String[]
    d86_each(k0) do it, o, raw, c, file
        p0 = k0.plan.panel.orbits[o]["p0"]
        (o in (2, 3, 4, 5) || (o in 8:11)) || return                                   # p_0 in {1.0, 2.0, 3.3, 5.0}
        fin, _, _ = d86_final(raw, o, "bk0")
        for r in 1:2
            (abs(fin[r, 1]) <= 4 * EX_EPS && abs(fin[r, 5]) <= 8 * EX_EPS) || push!(bad, "cos/sin orbit $o row $r")
        end
        abs(fin[3, 1] - p0) <= 2 * EX_EPS * p0 || push!(bad, "lifted p orbit $o")
    end
    @test isempty(bad)
end
