# Reference panel of Amendment D125 (prefix d125_ref_): the screening of one point by the verified Benettin on the identical orbit (E48), the class bits,
# the cell quantile, the draw, D, ids and seeds of a whole `systems` build (E49 to E53), written from the spec text and the lower layers only.
# Definitions only. Pins of the panel file and the pool sidecar are in the test plan, "Amendment D125 (stage 02)", readings U1, U6, U7.

d125_dim(c) = c["kind"] == "chain" ? 2 * c["N"] : 4
d125_cell_params(c) = c["kind"] == "froeschle" ? [c["nu"], c["kappa"], c["mu"]] : c["kind"] == "chain" ? [string(c["N"]), c["K"], c["eps"], c["boundary"]] : String[]
d125_flow_h(c) = parse(Float64, c["integrator"]["step"])                   # fixtures use decimal steps that are exact in binary64

# lambda and r of one screening point (maps: wrapped mode, flow: the cell's integrator map), by Benettin k = 1 with the seeded deviation of E48
function d125_ref_screen(c::AbstractDict, seed, j::Integer, x0::AbstractVector)
    d = d125_dim(c); N = c["scr_steps"]; tau = c["scr_tau"]; Nb = N ÷ c["scr_back"]
    sc = d94_sseed(d125_dev_text(seed, c["id"]))
    W = Chaoslab.initial_deviations(sc, j, Val(d), Val(1), Float64)
    sys, x = if c["kind"] == "henon_heiles"
        d125_flow_map(c["integrator"], Float64), Chaoslab.SVector{4, Float64}(x0...)
    else
        s = d125_system(c["kind"], d125_cell_params(c), Float64; wrapped = true)
        s, Chaoslab.torus_state(s, Chaoslab.SVector{d, Float64}(x0...))
    end
    function G(n)
        r = Chaoslab.Benettin()
        st = ex_Tg.propagate_orbit(sys, x, W, n, tau, (r,))[3]
        return Float64(only(ex_Tg.reducer_final(r, st[1]).lambda)) * n
    end
    GN, Gb = G(N), G(Nb)
    h = c["kind"] == "henon_heiles" ? d125_flow_h(c) : 1.0
    r = (isfinite(GN) && isfinite(Gb) && Gb > 0) ? GN / Gb : nothing
    return (lambda = GN / (N * h), r = r, t = N * h)
end

d125_rules_nt(r) = (ratio_regular_max = r["ratio_regular_max"], ratio_chaotic_min = r["ratio_chaotic_min"], ratio_weak_min = r["ratio_weak_min"],
                    weak_fraction = r["weak_fraction"])

# the whole screening of a cell: points, records, bits, Lambda_q (E46 to E50) and the counts of R33
function d125_ref_cell(fx_block::AbstractDict, c::AbstractDict)
    seed = fx_block["seed"]; M = c["n_start"]; rules = fx_block["rules"]
    pts, pairs = c["kind"] == "henon_heiles" ?
        d125_flow_starts(seed, c["id"], c["energy"], (c["q2_lo"], c["q2_hi"], c["p2_lo"], c["p2_hi"]), M) :
        (d125_points(seed, c["id"], d125_dim(c), M), 0)
    recs = [d125_ref_screen(c, seed, j, pts[j]) for j in 1:M]
    t = recs[1].t
    phi = d125_phi(rules["floor_abs"], rules["floor_log"], round(Int, t) == t ? Int(t) : error("t not an integer in the fixture"))
    nt = d125_rules_nt(rules)
    chaotic = [rc.lambda for rc in recs if rc.r !== nothing && rc.r >= nt.ratio_chaotic_min && rc.lambda >= phi]
    Lq = d125_Lq(chaotic, rules["weak_quantile"])
    bits = UInt8[d125_bits(rc.lambda, rc.r, phi, nt, Lq) for rc in recs]
    cnt(b) = count(x -> x & b != 0, bits)
    counts = (regular = cnt(0x01), chaotic = cnt(0x02), weak = cnt(0x04), strong = cnt(0x08), ambiguous = count(==(0x00), bits),
              no_ratio = count(rc -> rc.r === nothing, recs))
    return (pts = pts, pairs = pairs, recs = recs, bits = bits, Lq = Lq, counts = counts, phi = phi)
end

# draw, ids, D of a whole build from the reference cells: one entry per orbit in the order of E52, then the controls
function d125_ref_panel(blk::AbstractDict, cells::AbstractDict)
    seed = blk["seed"]; cls = (blk["classes"]["regular"], blk["classes"]["weak"], blk["classes"]["strong"])
    bitof = (0x01, 0x04, 0x08)
    orbs = NamedTuple[]; status = Dict{Tuple{String, String}, NamedTuple}()
    for c in blk["cell"]
        rc = cells[c["id"]]
        for (cid, bit) in zip(cls, bitof)
            pool = [j for j in eachindex(rc.bits) if rc.bits[j] & bit != 0]
            size = c["quota"][cid]
            drawn = d125_draw(seed, c["id"], cid, pool, size)
            din = d125_dindices(seed, c["id"], cid, length(drawn), c["d_fraction"][cid], c["d_min"])
            status[(c["id"], cid)] = (size = size, pool = length(pool), drawn = length(drawn), status = d125_status(size, length(drawn)))
            for (t, j) in enumerate(drawn)
                id = length(orbs) + 1
                push!(orbs, (id = id, cell = c["id"], class = cid, point_id = j, start = rc.pts[j], draw_index = t, in_D = t in din,
                             orbit_seed = string(EX.orbit_seed(seed, id)), run_seed = string(seed)))
            end
        end
    end
    return orbs, status
end

# the generic checks of a built panel against the reference (used for both fixtures, D125d to D125h)
function d125_check_panel(b; cells = nothing)
    blk = b.block
    cells === nothing && (cells = Dict(c["id"] => d125_ref_cell(blk, c) for c in blk["cell"]))
    orbs, status = d125_ref_panel(blk, cells)
    got = [o for o in b.panel["orbit"] if !get(o, "control", false)]
    @test length(got) == length(orbs)
    @test all(i -> got[i]["id"] == orbs[i].id && got[i]["cell"] == orbs[i].cell && got[i]["class"] == orbs[i].class &&
                   got[i]["point_id"] == orbs[i].point_id && got[i]["draw_index"] == orbs[i].draw_index && got[i]["in_D"] == orbs[i].in_D &&
                   got[i]["orbit_seed"] == orbs[i].orbit_seed && got[i]["run_seed"] == orbs[i].run_seed, eachindex(orbs))
    @test all(i -> Float64.(got[i]["start"]) == orbs[i].start, eachindex(orbs))                  # bitwise
    # controls after the drawn orbits, file order, exact exponents (R32)
    ctl = [o for o in b.panel["orbit"] if get(o, "control", false)]
    @test [o["cell"] for o in ctl] == fill("controls", length(blk["control"]))
    @test all(o -> o["class"] == "control" && !haskey(o, "draw_index") && !(o["in_D"]), ctl)
    @test [o["id"] for o in b.panel["orbit"]] == collect(1:(length(orbs) + length(ctl)))
    # manifest: per cell counts, Lambda_q, per class size, pool, drawn, status (R33)
    mc = b.panel["manifest"]["cell"]
    for c in blk["cell"]
        m = mc[c["id"]]; rc = cells[c["id"]]
        @test m["kind"] == c["kind"] && m["n_start"] == c["n_start"]
        @test Int(m["regular"]) == rc.counts.regular && Int(m["chaotic"]) == rc.counts.chaotic && Int(m["weak"]) == rc.counts.weak &&
              Int(m["strong"]) == rc.counts.strong && Int(m["ambiguous"]) == rc.counts.ambiguous && Int(m["no_ratio"]) == rc.counts.no_ratio
        rc.Lq === nothing ? @test(!haskey(m, "Lambda_q")) : @test(parse(Float64, m["Lambda_q"]) == rc.Lq)
        for cid in (blk["classes"]["regular"], blk["classes"]["weak"], blk["classes"]["strong"])
            s = status[(c["id"], cid)]; mm = m["class"][cid]
            @test (Int(mm["size"]), Int(mm["pool"]), Int(mm["drawn"]), mm["status"]) == (s.size, s.pool, s.drawn, s.status)
        end
        c["kind"] == "henon_heiles" && @test Int(m["pairs_consumed"]) == rc.pairs
    end
    # pool sidecar: bit vector hash and run-length lists (reading U6)
    for c in blk["cell"]
        sc = b.pools["cell"][c["id"]]; bits = cells[c["id"]].bits
        @test Int(sc["points"]) == length(bits) && sc["bits_sha256"] == ex_sha256(bits)
        for (name, bit) in (("regular", 0x01), ("chaotic", 0x02), ("weak", 0x04), ("strong", 0x08))
            @test d94_runs_ids(get(sc, name, [])) == [j for j in eachindex(bits) if bits[j] & bit != 0]
        end
    end
    return orbs
end

# a fixture block with n_start changed in every cell (the class bits of 200 points per kind, D125e)
d125_big_block(fx, n) = d125_edit(d125_block(fx), [("cell", i, "n_start") => n for i in eachindex(fx.cells)]...)
