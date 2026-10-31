# Helpers for the Amendment D94 suite (prefix d94_). Definitions only: nothing here calls the D94 API at load time, so a missing function
# errors the testsets that need it and nothing else. Reuses ex_helpers.jl and ex_run_helpers.jl without changing them.
# Spec: workspaces/module-build/stages/01-spec/output/experiments-spec.md, "Amendment D94". API assumed: the test plan, "Amendment D94 (stage 02)".
# Every reference below is written from the text of the spec (E32 to E43, R22 to R27), not from the source of the module.

d94_seed(letter::Char, sub::Integer = 0) = EX_SEED + 31000 + 50 * (letter - 'a' + 1) + sub     # testset D94<letter> uses 31000 + 50 n + sub, sub < 50

const D94_CLASSES = ("T", "R", "C", "S", "X")                  # numbering order (E36)
const D94_BIT = Dict("T" => 0x01, "C" => 0x02, "S" => 0x04)

# ---------------------------------------------------------------- streams by the independent SHA-256 (E32, E33, E37, E38)
const D94_SEEDS = Dict{String, UInt64}()
d94_u64(hex::AbstractString) = parse(UInt64, hex[1:16]; base = 16)                  # first 8 bytes, big-endian
function d94_prefetch_many(texts)
    todo = unique([t for t in texts if !haskey(D94_SEEDS, t)])
    isempty(todo) && return
    if Sys.which("sha256sum") === nothing
        foreach(t -> D94_SEEDS[t] = d94_u64(ex_sha256(t)), todo)
    else
        for lo in 1:100:length(todo)                          # the command line of one sha256sum call holds at most a few hundred names
            part = todo[lo:min(end, lo + 99)]
            for (t, h) in zip(part, ex_sha256_many(part))
                D94_SEEDS[t] = d94_u64(h)
            end
        end
    end
end
d94_sseed(text::AbstractString) = get!(() -> d94_u64(ex_sha256(text)), D94_SEEDS, String(text))
d94_stream(text) = StableRNG(d94_sseed(text))

d94_screen_text(seed, i) = "chaoslab-screen-1|$(seed)|$(i)"
d94_round_text(seed, i, c, r) = "chaoslab-panel-1|$(seed)|$(i)|$(c)|$(r)"
d94_d_text(seed, r, i, c) = "chaoslab-d-1|$(seed)|$(r)|$(i)|$(c)"
d94_db_text(seed, R, i) = "chaoslab-db-1|$(seed)|$(R)|$(i)"

# E32: points of K index i; u_1, u_2, ... in order, point j = (RN(2 pi u_(2j-1)), RN(2 pi u_(2j))), RN of the 256-bit value
function d94_points(seed, i, M)
    rng = d94_stream(d94_screen_text(seed, i))
    return setprecision(BigFloat, 256) do
        tp = 2 * BigFloat(pi)
        [(Float64(tp * BigFloat(rand(rng))), Float64(tp * BigFloat(rand(rng)))) for _ in 1:M]
    end
end

# E14 on a candidate list with a given stream (the draw of the base is the same code of the spec): partial Fisher-Yates, one uniform per draw
d94_fy(rng, cands, q) = ex_ref_draw(rng, cands, q)

# ---------------------------------------------------------------- exact arithmetic (I36)
d94_rat(x::Integer) = Rational{BigInt}(x)
function d94_rat(s::AbstractString)                  # the rational of the decimal text, or of a/b
    if occursin('/', s)
        a, b = split(s, '/')
        return Rational{BigInt}(parse(BigInt, strip(a)), parse(BigInt, strip(b)))
    end
    m = match(r"^(\d+)(?:\.(\d+))?$", strip(s))
    frac = m.captures[2] === nothing ? "" : m.captures[2]
    return Rational{BigInt}(parse(BigInt, m.captures[1] * frac), BigInt(10)^length(frac))
end
d94_rat(x::AbstractFloat) = d94_rat(string(x))
# m = floor(f n + 1/2) in integers: floor((2 a n + b) / (2 b)) for f = a / b
d94_dsize(f::Rational, n::Integer) = Int(div(2 * numerator(f) * n + denominator(f), 2 * denominator(f)))

# E38 largest remainder allocation; exact rationals for the remainders
function d94_alloc(n::AbstractVector{<:Integer}, b::Integer)
    T = sum(n)
    T <= b && return collect(Int, n)
    q = [floor(BigInt, Rational{BigInt}(b * ni, T)) for ni in n]
    rem = [Rational{BigInt}(b * ni, T) - qi for (ni, qi) in zip(n, q)]
    extra = b - sum(q)
    order = sort(collect(eachindex(n)); by = i -> (-rem[i], i))            # ties to the smaller K index
    out = collect(Int, q)
    for i in order[1:extra]
        out[i] += 1
    end
    return out
end

# E37: the draw indices of D of one round cell (n orbits, share f)
d94_dindices(seed, r, i, c, n, f) = sort(d94_fy(d94_stream(d94_d_text(seed, r, i, c)), 1:n, d94_dsize(f, n)))

# E38: D_B from D_R = [(id, k_index)], b and R
function d94_db(seed, R, DR, b)
    nK = isempty(DR) ? 0 : maximum(x.k_index for x in DR)
    n = [count(x -> x.k_index == i, DR) for i in 1:nK]
    al = d94_alloc(n, b)
    out = Int[]
    for i in 1:nK
        ids = sort([x.id for x in DR if x.k_index == i])
        al[i] == 0 && continue
        append!(out, d94_fy(d94_stream(d94_db_text(seed, R, i)), ids, al[i]))
    end
    return sort(out)
end

# ---------------------------------------------------------------- pools: class bit vectors and the draw of the pool classes (E34, E35)
function d94_bits(M::Integer; T = Int[], C = Int[], S = Int[])
    b = zeros(UInt8, M)
    for (c, ids) in (("T", T), ("C", C), ("S", S)), j in ids
        b[j] |= D94_BIT[c]
    end
    return b
end
d94_pool(bits, c) = [j for j in eachindex(bits) if bits[j] & D94_BIT[c] != 0]

# rounds of one K: the order inside a round is S, C, T; the used set holds every id drawn at that K; the recorded pool is the candidate count
# at the moment of the draw (reading pinned in the plan: after removing the used set, earlier classes of the same round included)
function d94_ref_pool_rounds(seed, i, bits, quota, r_max; active = 1:r_max)
    used = Set{Int}()
    cells = NamedTuple[]
    for r in 1:r_max, c in ("S", "C", "T")
        q = r in active ? get(quota, c, 0) : 0
        cand = [j for j in d94_pool(bits, c) if !(j in used)]
        drawn = q > 0 ? d94_fy(d94_stream(d94_round_text(seed, i, c, r)), cand, q) : Int[]
        union!(used, drawn)
        status = q == 0 ? "not drawn" : length(drawn) == q ? "populated" : isempty(drawn) ? "empty" : "short"
        push!(cells, (round = r, class = c, pool = length(cand), q = q, drawn = drawn, status = status))
    end
    return cells
end
d94_cell(cells, r, c) = only(x for x in cells if x.round == r && x.class == c)

# ---------------------------------------------------------------- the pool sidecar (R23 item 3), read by an independent reader
# layout pinned in the plan: [k.<i>] points, bits_sha256, T, C, S; each class a list of [first, last] runs of consecutive point ids
d94_runs_ids(runs) = isempty(runs) ? Int[] : sort(vcat([collect(Int(r[1]):Int(r[2])) for r in runs]...))
# run-length encoding of one class of a bit vector, the form the sidecar is read in (maximal runs of consecutive ids)
function d94_runs(bits, c)
    ids = d94_pool(bits, c)
    runs = Vector{Int}[]
    for j in ids
        !isempty(runs) && runs[end][2] == j - 1 ? (runs[end][2] = j) : push!(runs, [j, j])
    end
    return runs
end
function d94_side_bits(side, i)
    k = side["k"][string(i)]
    bits = zeros(UInt8, Int(k["points"]))
    for c in ("T", "C", "S"), j in d94_runs_ids(k[c])
        bits[j] |= D94_BIT[c]
    end
    return bits
end

# ---------------------------------------------------------------- fixtures of the builder
# A fixture is a NamedTuple; fixture A and fixture B differ in every value (D94t).
d94_qtab(q) = Dict{String, Any}(string(i) => Dict{String, Any}(c => n for (c, n) in row) for (i, row) in q)
d94_ftab(f) = Dict{String, Any}(string(i) => Dict{String, Any}(c => v for (c, v) in row) for (i, row) in f)

const D94_FA = (name = "A", seed = EX_SEED, K = ["0", "1.5", "5"], M = 64, r_max = 3, gates = [1], controls = "none", db = 6,
                quota = Dict(1 => Dict("T" => 3, "R" => 2), 2 => Dict("T" => 3, "R" => 2, "C" => 3, "S" => 2, "X" => 2),
                             3 => Dict("T" => 2, "R" => 2, "C" => 3, "S" => 2, "X" => 3)),
                dfrac = Dict(1 => Dict("T" => "0.4", "R" => "0.5"), 2 => Dict("T" => "0.5", "R" => "1", "C" => "0.35", "S" => "0.5", "X" => "0.5"),
                             3 => Dict("T" => "0.75", "R" => "0.5", "C" => "0.25", "S" => "0.5", "X" => "1/3")))
const D94_FB = (name = "B", seed = 777, K = ["0", "0.9", "3"], M = 48, r_max = 2, gates = [2], controls = "none", db = 3,
                quota = Dict(1 => Dict("T" => 2, "R" => 1), 2 => Dict("T" => 2, "R" => 1, "C" => 2, "S" => 1, "X" => 1),
                             3 => Dict("T" => 3, "R" => 1, "C" => 2, "S" => 2, "X" => 2)),
                dfrac = Dict(1 => Dict("T" => "0.5", "R" => "0.5"), 2 => Dict("T" => "0.2", "R" => "0.5", "C" => "0.5", "S" => "1", "X" => "0.5"),
                             3 => Dict("T" => "0.4", "R" => "0.5", "C" => "0.6", "S" => "0.5", "X" => "0.5")))

function d94_over(fx)
    return Dict{String, Any}("screen_source" => "random", "screen_points" => fx.M, "r_max" => fx.r_max, "gate_rounds" => fx.gates,
                             "quota" => d94_qtab(fx.quota), "d_fraction" => d94_ftab(fx.dfrac), "db_size" => fx.db, "controls" => fx.controls,
                             "seed" => fx.seed, "K" => fx.K)
end

const D94_BUILDS = Dict{String, Any}()
function d94_build(tag::AbstractString, fx; kw...)
    get!(D94_BUILDS, tag) do
        block = ex_block(d94_over(fx); drop = ["n_c"])
        pre = ex_path("d94-prereg-" * ex_sha256(ex_prereg_text(block))[1:12] * ".md")
        isfile(pre) || ex_write_prereg(pre, block)
        out = ex_path("d94-$tag-panel.toml")
        EX.build_panel(pre; out = out, kw...)
        (panel = TOML.parsefile(out), path = out, prereg = pre, pools = TOML.parsefile(out * ".pools.toml"), fx = fx)
    end
end

# lattice builds of the base fixtures (D94b): the same tags as ex_panel.jl, so the builds are shared
d94_lattice_builds() = Dict("base" => ex_build("base"), "k0-full" => ex_build("k0-full", Dict{String, Any}("K" => ["0"])),
                            "noell" => ex_build("noell", Dict{String, Any}("K" => ["10"], "q_max" => 1)),
                            "k0-empty" => ex_build("k0-empty", Dict{String, Any}("K" => ["0"], "theta_T" => 0.001)))
const D94_SKIP = ("git_sha", "git_dirty", "package_version", "julia_version", "stablerngs_version")
function d94_canon_sha(path)
    d = TOML.parsefile(path)
    foreach(k -> delete!(d["manifest"], k), D94_SKIP)
    io = IOBuffer()
    TOML.print(io, d; sorted = true)
    return ex_sha256(take!(io))
end

# ---------------------------------------------------------------- reference screening bits (the verified screening on the identical start, E11)
function d94_ref_bits(b, i)
    fx = b.fx
    blk = ex_block_default()
    K = Float64(ex_f64(fx.K[i]))
    ps = ex_pset(K; q_max = blk["q_max"], N_g = blk["grid_side"])
    alpha = ex_f64(blk["alpha"])
    pts = d94_points(fx.seed, i, fx.M)
    bits = zeros(UInt8, fx.M)
    for (j, (x, p)) in enumerate(pts)
        r = EX.screen_point(K, x, p; N_scr = blk["N_scr"], w = blk["w"], alpha = alpha, tau = 1)
        f = 0x00
        if r.lambda <= blk["theta_T"] && r.S_tenth > 0 &&
           EX.pred_T(r.lambda, r.S_tenth, r.S_N, EX.periodic_distance((x, p), ps); theta_T = blk["theta_T"], g_T = blk["g_T"], r_T = blk["r_T"])
            f |= 0x01
        end
        if r.lambda >= blk["theta_C"] && EX.pred_C(r.lambda, EX.periodic_distance((x, p), ps; elliptic_only = true); theta_C = blk["theta_C"], r_C = blk["r_C"])
            f |= 0x02
        end
        EX.pred_S(r.lambda, r.blocks; theta_S = blk["theta_S"], w = blk["w"], f_S = blk["f_S"], L_S = blk["L_S"]) && (f |= 0x04)
        bits[j] = f
    end
    return bits
end

# ---------------------------------------------------------------- reference of a whole built panel (E32 to E37) from its sidecar and manifest
# returns (orbits, cells): orbits in the order of E36 (without controls), cells = every (round, K index, class)
function d94_ref_panel(b)
    fx, man = b.fx, b.panel["manifest"]
    seed, nK = fx.seed, length(fx.K)
    blk = ex_block_default()
    pts = [d94_points(seed, i, fx.M) for i in 1:nK]
    bits = [d94_side_bits(b.pools, i) for i in 1:nK]
    texts = vec([d94_round_text(seed, i, c, r) for r in 1:fx.r_max, i in 1:nK, c in D94_CLASSES])
    append!(texts, vec([d94_d_text(seed, r, i, c) for r in 1:fx.r_max, i in 1:nK, c in D94_CLASSES]))
    d94_prefetch_many(texts)
    orbits, cells = NamedTuple[], NamedTuple[]
    pooled = [d94_ref_pool_rounds(seed, i, bits[i], Dict(c => get(get(fx.quota, i, Dict()), c, 0) for c in ("T", "C", "S")), fx.r_max;
                                  active = (fx.K[i] == "0" ? fx.gates : 1:fx.r_max)) for i in 1:nK]
    ell = [[(row["id"], k, row["points"][k]) for row in man["periodic"] if row["k_index"] == i && row["class"] == "elliptic" for k in eachindex(row["points"])] for i in 1:nK]
    frs = [(a, bb) for bb in 1:blk["q_max"] for a in 0:(bb - 1) if gcd(a, bb) == 1]
    lo, hi = blk["log10_r"]
    for r in 1:fx.r_max, i in 1:nK, c in D94_CLASSES
        zerok = fx.K[i] == "0"
        undefined = zerok && c in ("C", "S", "X")
        q = undefined ? 0 : (zerok && !(r in fx.gates)) ? 0 : get(get(fx.quota, i, Dict()), c, 0)
        items = NamedTuple[]
        pool = 0
        if c in ("T", "C", "S")
            pc = d94_cell(pooled[i], r, c)
            pool = pc.pool
            for (t, j) in enumerate(pc.drawn)
                push!(items, (t = t, point_id = j, x0 = pts[i][j][1], p0 = pts[i][j][2], source = nothing))
            end
        elseif c == "R" && q > 0
            rng = d94_stream(d94_round_text(seed, i, "R", r))
            for t in 1:(zerok || !isempty(ell[i]) ? q : 0)
                u1, u2 = rand(rng), rand(rng)
                if zerok
                    a, bb = frs[1 + floor(Int, u1 * length(frs))]
                    push!(items, (t = t, point_id = 0, x0 = 2pi * u2, p0 = setprecision(() -> Float64(2 * BigFloat(pi) * a / bb), BigFloat, 1024), source = (a, bb)))
                else
                    k = 1 + floor(Int, u1 * length(ell[i]))
                    push!(items, (t = t, point_id = 0, x0 = NaN, p0 = NaN, source = (ell[i][k][1], ell[i][k][2], 2pi * u2)))
                end
            end
        elseif c == "X" && q > 0
            rng = d94_stream(d94_round_text(seed, i, "X", r))
            for t in 1:q
                u1, u2 = rand(rng), rand(rng)
                rr = 10.0^(lo + (hi - lo) * u1)
                phi = 2pi * u2
                push!(items, (t = t, point_id = 0, x0 = rr * cos(phi), p0 = rr * sin(phi), source = (rr, phi)))
            end
        end
        status = undefined ? "undefined" : q == 0 ? "not drawn" : length(items) == q ? "populated" : isempty(items) ? "empty" : "short"
        push!(cells, (round = r, k_index = i, class = c, q = q, drawn = length(items), status = status, pool = pool, items = items))
        f = (haskey(fx.dfrac, i) && haskey(fx.dfrac[i], c)) ? d94_rat(fx.dfrac[i][c]) : d94_rat(0)
        sel = Set(d94_dindices(seed, r, i, c, length(items), f))
        for it in items
            push!(orbits, (; it..., round = r, k_index = i, class = c, in_D = it.t in sel))
        end
    end
    return (orbits = [merge(o, (id = n,)) for (n, o) in enumerate(orbits)], cells = cells)
end

# ---------------------------------------------------------------- hand panel for the plan and the runner (no builder): round, in_D, no controls
# rounds: (round, K index, class, count, in_D draw indices); ids run in the order of E36
const D94_CELLS = [(1, 1, "T", 2, [1]), (1, 2, "T", 2, [1]), (1, 2, "C", 2, [1]),
                   (2, 2, "T", 2, [1]), (2, 2, "C", 1, [1]), (2, 3, "S", 2, [1]),
                   (3, 2, "T", 2, [1]), (3, 2, "C", 2, [1]), (3, 3, "S", 1, [1])]
const D94_KTEXT = ["0.5", "1.5", "2.5"]

function d94_panel_orbits(cells = D94_CELLS; with_D::Bool = true)
    rng = StableRNG(d94_seed('l'))
    orbs = Dict{String, Any}[]
    for (r, i, c, n, din) in cells, t in 1:n
        id = length(orbs) + 1
        x, p = 2pi * rand(rng), 2pi * rand(rng)
        o = Dict{String, Any}("id" => id, "class" => c, "control" => false, "control_kind" => "none", "k_index" => i, "K_text" => D94_KTEXT[i],
                              "K" => parse(Float64, D94_KTEXT[i]), "system" => "standard", "x0" => x, "p0" => p, "x0_dec" => string(x), "p0_dec" => string(p),
                              "draw_index" => t, "run_seed" => string(EX_RUN_SEED), "orbit_seed" => string(EX.orbit_seed(EX_RUN_SEED, id)),
                              "in_SD" => false, "in_SL" => false, "round" => r)
        with_D && (o["in_D"] = t in din)
        push!(orbs, o)
    end
    return orbs
end

function d94_write_panel(path::AbstractString; orbits = d94_panel_orbits(), r_max = 3, db_size = 5, accepted = EX_PANEL_HASH)
    man = Dict{String, Any}("format_version" => 1, "controls" => Int[], "r_max" => r_max, "gate_rounds" => [1], "db_size" => db_size,
                            "prereg" => Dict{String, Any}("header_sha256" => accepted, "computed_sha256" => accepted))
    r_max === nothing && delete!(man, "r_max")
    open(path, "w") do io
        TOML.print(io, Dict{String, Any}("manifest" => man, "orbit" => orbits))
    end
    return path
end

# run blocks of the D94 fixture: each_round arm "ref" on P, closure arms "dd" on D and "db" on D_B (different methods, so no item merges)
d94_arms(; rel = true) = [ex_arm("ref", "P", ["L1"], ["Float64"], [20]; (rel ? (; release = "each_round") : (;))...),
                          ex_arm("dd", "D", ["L5"], ["Float64"], [20]; (rel ? (; release = "closure") : (;))...),
                          ex_arm("db", "D_B", ["L4"], ["Float64"], [20]; (rel ? (; release = "closure") : (;))...)]
function d94_block(arms = d94_arms())
    d = ex_rb_common()
    d["pilot"]["horizon"] = 20
    d["arm"] = arms
    return d
end

const D94_FIXES = Dict{String, Any}()
function d94_fix(tag::AbstractString, run = d94_block(); rev::Bool = false, orbits = d94_panel_orbits(), r_max = 3, db_size = 5)
    return get!(D94_FIXES, tag) do
        pp = d94_write_panel(ex_path("d94-$tag-panel.toml"); orbits = orbits, r_max = r_max, db_size = db_size)
        pr = ex_write_prereg_run(ex_path("d94-$tag-prereg.md"), run; rev = rev)
        (; run, panel = pp, prereg = pr, tag, orbits)
    end
end
d94_std() = d94_fix("std")

# ---------------------------------------------------------------- the round log (E40): text format pinned in the plan
const D94_EMPTY_SHA = "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855"      # SHA-256 of the empty text
function d94_entry(; round, decision, panel, prereg, digest, prev, counts = Dict("n_ts" => 11, "n_w2" => 4))
    io = IOBuffer()
    println(io, "round=", round)
    println(io, "decision=", decision)
    for k in sort(collect(keys(counts)))
        println(io, "count.", k, "=", counts[k])
    end
    println(io, "panel_sha256=", panel)
    println(io, "prereg_sha256=", prereg)
    println(io, "ref_digest=", digest)
    println(io, "prev_sha256=", prev)
    return String(take!(io))
end

# digest of the ref arm: SHA-256 of the sorted chunk hashes (each followed by a newline)
d94_digest(hashes) = ex_sha256(join(sort(collect(hashes)) .* "\n"))

# a log from specs (round, decision, digest-or-nothing); `digest(r)` gives the right digest of the rounds 1 to r; mutations by keyword
function d94_log_text(fx, specs, digest; panel = nothing, prereg = nothing)
    pnl = panel === nothing ? ex_sha256(read(fx.panel)) : panel
    pre = prereg === nothing ? ex_header_hash(read(fx.prereg, String)) : prereg
    ents = String[]
    prev = D94_EMPTY_SHA
    for s in specs
        e = d94_entry(; round = s.round, decision = s.decision, panel = get(s, :panel, pnl), prereg = get(s, :prereg, pre),
                      digest = get(s, :digest, digest(s.round)), prev = get(s, :prev, prev),
                      counts = get(s, :counts, Dict("n_ts" => 11, "n_w2" => 4)))
        push!(ents, e)
        prev = ex_sha256(e)
    end
    return join(ents), ents
end
function d94_write_log(path, fx, specs, digest; kw...)
    txt, ents = d94_log_text(fx, specs, digest; kw...)
    write(path, txt)
    return path
end

# raw chunk file root attribute (independent of the package reader)
function d94_attr(file, name)
    H = EX.HDF5
    return H.h5open(f -> H.read_attribute(f, name), String(file), "r")
end
d94_ordinal(f) = parse(Int, match(r"chunk-(\d{5})\.h5", basename(f)).captures[1])
d94_cell_round(cell) = parse(Int, match(r":r(\d+)$", cell).captures[1])

# ref arm chunk hashes of the rounds up to r from the files of a store (root attributes chunk_hash, chunk_key, arms)
function d94_ref_hashes(store, r)
    out = String[]
    for f in ex_chunk_files(store)
        arms = d94_attr(f, "arms")
        "ref" in arms || continue
        m = match(r"cell=\S+:r(\d+)$", String(d94_attr(f, "chunk_key")))
        m !== nothing && parse(Int, m.captures[1]) <= r && push!(out, String(d94_attr(f, "chunk_hash")))
    end
    return out
end

d94_ref_chunks(plan) = [c for c in plan.chunks if "ref" in c.arms]
d94_clo_chunks(plan) = [c for c in plan.chunks if !("ref" in c.arms)]

# ---------------------------------------------------------------- the plan and the campaign of the fixture
d94_plan(fx; log = nothing, kw...) = EX.plan_run(fx.prereg, fx.panel; code_version = EX_CODE, round_log = log, kw...)

# the shared campaign: round by round through the log, written by the test (stage = 0..3 chunks of the stage are run and checked)
const D94_CAMP = Ref{Any}(nothing)
function d94_campaign()
    D94_CAMP[] === nothing || return D94_CAMP[]
    fx = d94_std()
    store = ex_store("d94-camp")
    logp = ex_path("d94-camp.log")
    write(logp, "")
    reads = NamedTuple[]
    hook = EX.RunHooks(clock = ex_clock(info -> 0.0), read = (kind, path) -> push!(reads, (kind = kind, path = String(path))))
    specs = NamedTuple[]
    digest = r -> d94_digest(d94_ref_hashes(store, r))
    stages = NamedTuple[]
    for (stage, dec) in ((1, "continue"), (2, "continue"), (3, "close"))
        n0 = length(reads)
        res = ex_campaign(fx; store = store, round_log = logp, hooks = hook)
        files = ex_chunk_files(store)
        push!(stages, (; stage, res, files = copy(files), reads = reads[(n0 + 1):end]))
        push!(specs, (round = stage, decision = dec))
        d94_write_log(logp, fx, specs, digest)
    end
    n0 = length(reads)
    res = ex_campaign(fx; store = store, round_log = logp, hooks = hook)
    push!(stages, (; stage = 4, res, files = copy(ex_chunk_files(store)), reads = reads[(n0 + 1):end]))
    D94_CAMP[] = (; fx, store, logp, stages, specs, digest)
    return D94_CAMP[]
end

# reference counts of items and chunks of a hand-panel plan from the panel itself (no plan code): per phase and round
function d94_ref_counts(R_plan::Integer, closed::Bool, dbids)
    ref_items = Dict(r => 0 for r in 1:R_plan)
    ref_chunks = Dict(r => 0 for r in 1:R_plan)
    for (r, i, c, n, din) in D94_CELLS
        r <= R_plan || continue
        ref_items[r] += n
        ref_chunks[r] += 1
    end
    return ref_items, ref_chunks
end

# D and D_B ids of the hand panel (E36 order; the D flags of D94_CELLS)
function d94_hand_D(R)
    ids = Int[]
    k = Int[]
    n = 0
    for (r, i, c, cnt, din) in D94_CELLS, t in 1:cnt
        n += 1
        if t in din && r <= R
            push!(ids, n)
            push!(k, i)
        end
    end
    return [(id = a, k_index = b) for (a, b) in zip(ids, k)]
end
d94_hand_cell(id) = (n = 0; for (r, i, c, cnt, din) in D94_CELLS, t in 1:cnt; (n += 1) == id && return (r, i, c); end; nothing)

# an ArgumentError whose message holds all the texts
function d94_rej(f, texts::AbstractString...)
    try
        f()
    catch e
        return e isa ArgumentError && all(t -> occursin(t, e.msg), texts)
    end
    return false
end
d94_rej_any(f, texts) = ex_rejects(f, texts)

# the generic checks of a built panel against the independent reference (used for both fixtures, D94t)
function d94_check_panel(b)
    fx = b.fx
    ref = d94_ref_panel(b)
    orbs = [o for o in b.panel["orbit"] if !o["control"]]
    man = b.panel["manifest"]
    @test length(orbs) == length(ref.orbits)
    @test [o["id"] for o in b.panel["orbit"]] == collect(1:length(b.panel["orbit"]))                 # no gap, E36 order
    for (o, r) in zip(orbs, ref.orbits)
        @test (o["round"], o["k_index"], o["class"], o["draw_index"], o["id"]) == (r.round, r.k_index, r.class, r.t, r.id)
        @test o["in_D"] == r.in_D
        @test o["orbit_seed"] == string(EX.orbit_seed(fx.seed, o["id"]))
        if r.class in ("T", "C", "S")
            @test o["point_id"] == r.point_id && ex_start(o) == (r.x0, r.p0)                           # bitwise
        elseif r.class == "X"
            @test !haskey(o, "point_id") && abs(o["x0"] - r.x0) <= 1e-12 * abs(r.source[1]) && abs(o["p0"] - r.p0) <= 1e-12 * abs(r.source[1])
            @test abs(o["r"] - r.source[1]) <= 1e-12 * r.source[1] && abs(o["phi"] - r.source[2]) <= 1e-12
        elseif r.class == "R" && fx.K[r.k_index] == "0"
            @test (o["a"], o["b"]) == r.source && o["p0"] == r.p0 && abs(o["x0"] - r.x0) <= 1e-12
        elseif r.class == "R"
            @test (o["source_orbit"], o["source_point"]) == (r.source[1], r.source[2]) && abs(o["phi"] - r.source[3]) <= 1e-12
        end
    end
    # round cells of the manifest: pool (pool classes), q, drawn, status, and the recount of the orbits
    cells = man["round_cells"]
    @test length(cells) == length(ref.cells)
    for c in ref.cells
        row = only(x for x in cells if (x["round"], x["k_index"], x["class"]) == (c.round, c.k_index, c.class))
        @test (row["q"], row["drawn"], row["status"]) == (c.q, c.drawn, c.status)
        c.class in ("T", "C", "S") && c.status != "undefined" && c.status != "not drawn" && @test row["pool"] == c.pool
        @test row["drawn"] == count(o -> (o["round"], o["k_index"], o["class"]) == (c.round, c.k_index, c.class), orbs)
        row["status"] in ("short", "empty") && @test !isempty(row["reason"])
        row["status"] == "populated" && @test row["reason"] == ""
    end
    return ref
end
