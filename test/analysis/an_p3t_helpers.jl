# Amendment D95, stage 02: helpers of the test side of the p3-abstention suite (prefix p3t_). Definitions only: nothing here calls Chaoslab.Analysis at load time (except the
# existing D93 JSON writer `AN.p3_json`, used by the record builder), so a missing function errors the testsets that need it and nothing else.
# Spec: workspaces/module-build/stages/01-spec/output/analysis-spec.md, "Amendment D95". API assumed: analysis-test-plan.md, "Amendment D95 (stage 02)", "Contract proposed by the tests".
#
# The toy campaign (scaled, 31 orbits, 3 rounds): the orbits are copies of designed orbits of the D93 fixture (templates: grid orbits 21 to 43, whose labels rho_new and verdict
# histories are known), written by the real chunk writer in the SHAPE OF THE RUNNER (Amendment D97, E82; before D97 the toy had the shape of the module, D95). One chunk per cell
# `k<K index>:<class>:r<round>` and item, in six arms: ref (L1 and L5, Float64, N_1, every orbit), ind (tg, fq1, bk1, Float64, N_c, every orbit), ind_pa1 (pa1, Float64, N_c + 1),
# d64 (tg, fq1, bk1, Double64, N_c, the orbits of D), d64_pa1 (pa1, Double64, N_c + 1), chk (L5 in Double64 and BigFloat(256), at its own horizon N_B, the orbits of D_B). The keys of
# the chunks and of the manifest entries are those of the runner, the manifest status is `complete with failed chunks` when a chunk failed, the specs sit on the indicator arms only.
# The shape is guarded by an_p3t_shape.jl against a real campaign of the runner. The panel file carries the rounds (D94 R23). Uses an_helpers.jl, the p2 helpers and the D93 helpers
# (an_p3_*.jl) without changing them.

p3t_seed(letter::Char, sub::Integer = 0) = AN_SEED + 36000 + 50 * (letter - 'a' + 1) + sub      # D95<letter> draws use 36000 + 50 n + sub, sub < 50
const P3T_KTEXT = ["0", "0.5", "1.5", "2.5", "3.5"]                  # K text of the K index 1 to 5 (index 1 is the K = 0 gate cell)
const P3T_A1 = "Amendment A1 (freeze of the calibration, 2026-11-01)"
const P3T_KEYS = [:pilot_cost_l1, :pilot_cost_tg, :pilot_cost_pa1, :pilot_cost_bk1, :pilot_cost_fq1, :calibration_record_sha256, :p3_table_sha256, :target_t_s, :target_w2, :target_f,
                  :theta_voters, :theta_q_c, :theta_q_r, :theta_w, :theta_v, :power_hb_planned, :power_hc_planned, :cap_flag]
const P3T_PILOT = [:pilot_cost_l1, :pilot_cost_tg, :pilot_cost_pa1, :pilot_cost_bk1, :pilot_cost_fq1]
const P3T_THETA = (V = [:SALI, :FLI, :MEGNO, :E2], qC = 2, qR = 3, w = 1, v = 0)            # theta*: four voters, depth 1 (b_P is the budget index 6, so b_P / 2 is in B)
const P3T_COUNT_KEYS = ["drift", "excluded", "orbits", "t_s", "target_t_s", "target_w2", "w2"]   # ascending byte order (E62 a)

# the 71-key table of the toy campaign: scaled floors and yields so that the campaign meets them (n_t_min 3, n_w_min 5, yields 3 and 5 per round), 100 resamples, r_max 3
p3t_q(; kw...) = p3_params(; r_max = 3.0, n_t_min = 3, n_w_min = 5, y_t = 3, y_w = 5, bootstrap_resamples = 100, pp_runs = 40, kw...)

# the values of the test-side table (A1); the targets decide the rounds of the log (see p3t_counts_by_round)
p3t_values(; kw...) = merge((; pilot_cost_l1 = 210.3, pilot_cost_tg = 2701.0, pilot_cost_pa1 = 1790.2, pilot_cost_bk1 = 3050.8, pilot_cost_fq1 = 16990.4,
    calibration_record_sha256 = "0"^64, p3_table_sha256 = "0"^64, target_t_s = 5, target_w2 = 11, target_f = 6, theta_voters = String.(P3T_THETA.V), theta_q_c = 2, theta_q_r = 3,
    theta_w = 1, theta_v = 0, power_hb_planned = 0.31, power_hc_planned = 0.82, cap_flag = false), NamedTuple(kw))

# ---------------------------------------------------------------- the designed orbits
# (id, template, round, K index, tag): tag in :n (normal), :gate, :drift, :U (L5 item flagged "diagnostics exceeded"), :excl (ref chunk failed twice), :nov (ind chunk failed twice)
const P3T_LIST = [
    (1, 27, 1, 1, :gate), (2, 28, 1, 1, :gate), (3, 23, 1, 2, :n), (4, 24, 1, 2, :n), (5, 21, 1, 2, :n), (6, 33, 1, 3, :n), (7, 34, 1, 3, :n), (8, 27, 1, 3, :n), (9, 29, 1, 2, :n),
    (10, 40, 1, 4, :n), (11, 36, 1, 4, :n), (12, 41, 1, 4, :n), (13, 22, 1, 2, :n), (14, 35, 1, 3, :drift),
    (15, 25, 2, 2, :n), (16, 39, 2, 2, :n), (17, 26, 2, 3, :n), (18, 37, 2, 3, :n), (19, 38, 2, 4, :n), (20, 30, 2, 4, :n), (21, 31, 2, 2, :n), (22, 42, 2, 3, :n), (23, 24, 2, 4, :U),
    (24, 21, 2, 5, :excl),
    (25, 23, 3, 2, :n), (26, 21, 3, 5, :nov), (27, 34, 3, 4, :n), (28, 26, 3, 2, :n), (29, 43, 3, 3, :n), (30, 22, 3, 4, :n), (31, 33, 3, 2, :n)]
const P3T_D = [3, 4, 5, 6, 7, 8, 9, 12, 13, 14, 15, 16, 18, 20, 21, 25]       # in_D (14 is a drift orbit, 8 and 9 and 20 and 21 are regular, 3 4 15 16 25 are sticky T_s orbits)
const P3T_DB = [3, 6, 8, 12]                                                  # D_B: the orbits of the arm chk (all in D)
const P3T_ROUNDS = 3
const P3T_RCLASS = [9, 13, 28]                                                 # panel class R (all others C): the cells of a round and a K index split by class (E74)
const P3T_ARMS = ["ref", "ind", "ind_pa1", "d64", "d64_pa1", "chk"]            # the six arms of the test run (E81)
const P3T_NB = 1000                                                            # the horizon of the arm chk = N_B of the P1 table of the fixture (an_fx_params)
const P3T_GC = sort(unique(vcat(filter(<=(P3T_NB), P2_GA), [P3T_NB])))         # recorded steps of the chk items
# the L5 values of the chk items are the reference value plus an offset (Double64, BigFloat): orbit 6 flips the ladder (0.02 > 5e-3), orbit 8 differs below the bound, orbits 3 and 12 equal
const P3T_CHK_OFF = Dict{Int, Tuple{Float64, Float64}}(3 => (0.0, 0.0), 6 => (0.0, 0.02), 8 => (0.004, 0.004), 12 => (0.0, 0.0))
# the Float64 votes of orbit 12 (a regular template) are changed to three chaotic ones at b_P (a failure of HC), and the Double64 verdicts of some voters differ at b_P (s_P > 0)
const P3T_OV = Dict{Tuple{Int, Symbol}, String}((12, :SALI) => "rrrrrrcccc", (12, :FLI) => "rrrrrrcccc", (12, :MEGNO) => "rrrrrrcccc")
const P3T_OVD = Dict{Tuple{Int, Symbol}, Tuple{Int, Char}}((3, :SALI) => (7, 'c'), (4, :SALI) => (7, 'c'), (15, :FLI) => (7, 'r'), (15, :MEGNO) => (7, 'c'), (25, :E2) => (7, 'r'),
    (6, :SALI) => (7, 'r'), (8, :MEGNO) => (7, 'c'))                        # (orbit, voter) => (position, verdict) of the Double64 sequence

# the design of the campaign: the D93 design remapped to the new ids, the Float64 design and the Double64 design (verdict sequences with the flips)
function p3t_design()
    base = p3_design()
    rho = Dict{Int, Symbol}(); lref = Dict{Int, Float64}(); lns = Dict{Int, Float64}(); V = Dict{Int, Dict{Symbol, Vector{Symbol}}}(); Vd = Dict{Int, Dict{Symbol, Vector{Symbol}}}()
    drift = Int[]
    for (id, t, r, k, tag) in P3T_LIST
        rho[id] = tag === :U ? :U : P3_RHO_NEW[t]
        lref[id] = base.lref[t]; lns[id] = base.lns[t]
        V[id] = Dict{Symbol, Vector{Symbol}}(X => copy(base.V[t][X]) for X in keys(base.V[t]))
        for ((i, X), s) in P3T_OV
            i == id && (V[id][X] = p2_sym(s))
        end
        Vd[id] = Dict{Symbol, Vector{Symbol}}(X => copy(v) for (X, v) in V[id])
        for ((i, X), (pos, c)) in P3T_OVD
            i == id && (Vd[id][X][pos] = c == 'c' ? :chaotic : :regular)
        end
        tag === :drift && push!(drift, id)
    end
    mk(Vx) = (variant = :p3t, prof = Dict{Int, Symbol}(), rho = rho, lref = lref, lns = lns, V = Vx, d64 = Dict{Tuple{Int, Symbol}, Symbol}(),
              dir = Dict{Tuple{Int, Symbol, Int}, Symbol}(), drift = drift, nan_nu = Dict{Int, Vector{Int}}())
    return (f = mk(V), d = mk(Vd))
end

# the panel orbits of the toy campaign (start points from a recorded seed)
function p3t_orbits()
    rng = StableRNG(p3t_seed('a'))
    out = NamedTuple[]
    for (id, t, r, k, tag) in P3T_LIST
        x0, p0 = 2pi * rand(rng), 2pi * rand(rng)
        push!(out, (id = id, tmpl = t, round = r, kidx = k, tag = tag, x0 = x0, p0 = p0, ktext = P3T_KTEXT[k], ic = "$(x0),$(p0)", in_D = id in P3T_D, in_DB = id in P3T_DB,
                    cls = id in P3T_RCLASS ? "R" : "C"))
    end
    return out
end
p3t_info(o) = (x0 = o.x0, p0 = o.p0, ic = o.ic, params = o.ktext, cls = o.cls, system = "standard", kind = "none", K = parse(Float64, o.ktext), id = o.id)
p3t_cell(o) = "k$(o.kidx):$(o.cls):r$(o.round)"                       # the cell of an orbit (E74)

function p3t_write_panel(path, orbits, h; r_max = P3T_ROUNDS, drop = Symbol[], extra = Dict{String, Any}())
    man = Dict{String, Any}("format_version" => 1, "controls" => Int[], "r_max" => r_max, "gate_rounds" => [1], "db_size" => length(P3T_DB),
                            "prereg" => Dict{String, Any}("header_sha256" => h, "computed_sha256" => h))
    tab = Dict{String, Any}[]
    for o in orbits
        d = Dict{String, Any}("id" => o.id, "class" => o.cls, "control" => false, "control_kind" => "none", "k_index" => o.kidx, "K_text" => o.ktext, "K" => parse(Float64, o.ktext),
            "system" => "standard", "x0" => o.x0, "p0" => o.p0, "x0_dec" => string(o.x0), "p0_dec" => string(o.p0), "draw_index" => o.id, "run_seed" => "20261030",
            "orbit_seed" => string(2000 + o.id), "in_SD" => false, "in_SL" => false, "round" => o.round, "in_D" => o.in_D)
        for k in drop; delete!(d, string(k)); end
        merge!(d, extra)
        push!(tab, d)
    end
    open(io -> TOML.print(io, Dict{String, Any}("manifest" => man, "orbit" => tab)), path, "w")
    return path
end

# ---------------------------------------------------------------- prereg text with the pilot sentence and the Amendment A1 section; the test-side table with the quotes
p3t_texts(t) = Dict{Symbol, String}(
    :pilot_cost_l1 => string(t.pilot_cost_l1), :pilot_cost_tg => string(t.pilot_cost_tg), :pilot_cost_pa1 => string(t.pilot_cost_pa1), :pilot_cost_bk1 => string(t.pilot_cost_bk1),
    :pilot_cost_fq1 => string(t.pilot_cost_fq1), :calibration_record_sha256 => t.calibration_record_sha256, :p3_table_sha256 => t.p3_table_sha256, :target_t_s => string(t.target_t_s),
    :target_w2 => string(t.target_w2), :target_f => string(t.target_f), :theta_voters => join(t.theta_voters, ", "), :theta_q_c => string(t.theta_q_c), :theta_q_r => string(t.theta_q_r),
    :theta_w => string(t.theta_w), :theta_v => string(t.theta_v), :power_hb_planned => string(t.power_hb_planned), :power_hc_planned => string(t.power_hc_planned),
    :cap_flag => string(t.cap_flag))
p3t_section(k) = k in P3T_PILOT ? P3_SEC.COST : P3T_A1
function p3t_quotes(t; texts = Dict{Symbol, String}())
    T = merge(p3t_texts(t), texts)
    return Dict{Symbol, Tuple{String, String, String}}(k => (p3t_section(k), k in P3T_PILOT ? "pilot cost of the p2 pilot, $(k) = $(T[k]) ns per orbit step" : "$(k) = $(T[k])", T[k]) for k in P3T_KEYS)
end

# the preregistration of the campaign: the sections of the D93 fixture, the pilot sentence in the cost section, the section Amendment A1 (SELF hash rule as p3_prereg_text)
function p3t_prereg_text(q, t; slug = P3_SLUG, hashrow = true, a1 = true, texts = Dict{Symbol, String}())
    qq = p3_quotes(q); tq = p3t_quotes(t; texts)
    body = IOBuffer()
    println(body, "# Preregistration: ", slug, "\n\n## Header\n")
    println(body, "| Field | Value |\n|-------|-------|\n| Study slug | ", slug, " |")
    hashrow && println(body, "| SHA-256 of this file (current version) | ", "0"^64, " (rule as p2: fixture) |")
    println(body, "| Amendment log | A1 |\n\nFixture. Not a study.\n")
    for sk in P3_SEC_ORDER
        println(body, "## ", getfield(P3_SEC, sk), "\n")
        for k in P3_PARAM_KEYS
            qq[k][1] == getfield(P3_SEC, sk) && println(body, "- ", qq[k][2])
        end
        if getfield(P3_SEC, sk) == P3_SEC.COST
            for k in P3T_PILOT; println(body, "- ", tq[k][2]); end
        end
        println(body)
    end
    if a1
        println(body, "## ", P3T_A1, "\n")
        for k in P3T_KEYS
            k in P3T_PILOT || println(body, "- ", tq[k][2])
        end
        println(body)
    end
    text = String(take!(body))
    prep(s) = replace(replace(s, "\r\n" => "\n"), r"(?m)^(\| SHA-256 of this file \(current version\) \|)[^\n|]*\|" => s"\1 SELF |")
    h = p3ref_sha256(prep(text))
    hashrow && (text = replace(text, "0"^64 => h; count = 1))
    return text, h
end

# the test-side table (E60): value, section, quote, text per key; edit[k] = Dict(field => value or nothing to remove)
function p3t_write_values(path, t, q; drop = Symbol[], extra = Dict{String, Any}(), edit = Dict{Symbol, Any}(), texts = Dict{Symbol, String}())
    tq = p3t_quotes(t; texts)
    d = Dict{String, Any}()
    for k in P3T_KEYS
        k in drop && continue
        v = getfield(t, k)
        e = Dict{String, Any}("value" => v, "section" => tq[k][1], "quote" => tq[k][2], "text" => tq[k][3])
        haskey(edit, k) && for (kk, vv) in edit[k]; vv === nothing ? delete!(e, kk) : (e[kk] = vv); end
        d[string(k)] = e
    end
    merge!(d, extra)
    open(io -> TOML.print(io, d), path, "w")
    return path
end

# ---------------------------------------------------------------- the store, the manifest and the panel
# the items of an arm (E75, E81): `kind` is the arm name; ref and chk have no indicator specs, the indicator arms hold the specs of the runner (ids tg, fq1, bk1, pa1)
function p3t_specs(kind::Symbol, p)
    S = NamedTuple[]
    add(role, id, prec, hor, grid, arm; kd = "") = push!(S, (; role, id, precision = prec, horizon = hor, grid, arms = [arm], dindex = 1, direction = P2_DIRS[1], kind = kd))
    if kind === :ref
        add(:L1, "L1", "Float64", P2_N1, P2_GA, "ref"); add(:L5, "L5", "Float64", P2_N1, P2_GA, "ref")
    elseif kind in (:ind, :d64)
        pr = kind === :ind ? "Float64" : "Double64"; a = string(kind)
        add(:tg, "tg", pr, P2_NC, P2_GB, a; kd = "tangent"); add(:fq1, "fq1", pr, P2_NC, P2_GB, a; kd = "frequency"); add(:bk1, "bk1", pr, P2_NC, P2_GB, a; kd = "birkhoff")
    elseif kind in (:ind_pa1, :d64_pa1)
        add(:pa1, "pa1", kind === :ind_pa1 ? "Float64" : "Double64", P2_NC + 1, P2_GB1, string(kind); kd = "phase_advance")
    else
        add(:L5, "L5", "Double64", P3T_NB, P3T_GC, "chk"); add(:L5, "L5", "BigFloat(256)", P3T_NB, P3T_GC, "chk")
    end
    return S
end
const P3T_KINDS = [:ref, :ind, :ind_pa1, :d64, :d64_pa1, :chk]
p3t_allspecs(p) = vcat((p3t_specs(k, p) for k in P3T_KINDS)...)

# the key of a chunk as the runner writes it: `method=<id>[ kind=<kind> spec=<16 hex>] precision=<p> horizon=<h> tau=1 eps=- direction=<d or -> mode=wrapped cell=<cell>`
function p3t_key(sp, cell)
    spec = ""
    if sp.kind != ""
        kd, _ = P2_SPEC_INFO[sp.id]
        spec = " kind=$(kd) spec=" * p3ref_sha256("indicator.id = \"$(sp.id)\"\nindicator.kind = \"$kd\"\n")[1:16]
    end
    return "method=$(sp.id)$(spec) precision=$(sp.precision) horizon=$(sp.horizon) tau=1 eps=- direction=$(sp.role in (:L1, :L5, :tg) ? sp.direction : "-") mode=wrapped cell=$(cell)"
end

# the root attributes of a chunk (as the p2 fixture), with the cell in the chunk key; a failed chunk has status "failed" and attempt 2
function p3t_root(sp, ord, p, ph, h, cell; status = "ok", attempt = 1, qh = nothing)
    key = p3t_key(sp, cell)
    ch = p3ref_sha256("p3t-chunk-$ord-$key")
    root = Dict{String, Any}("schema_version" => 2, "study" => P3_SLUG, "chunk" => ord, "chunk_key" => key, "chunk_hash" => ch, "plan_sha256" => p.plan_sha256,
        "panel_sha256" => ph, "prereg_sha256" => h, "git_sha" => "a"^40, "code_version" => "an-code-1", "backend_request" => "cpu", "backend" => "cpu",
        "device" => "cpu", "backend_fallback_reason" => "", "precision" => sp.precision, "method" => sp.id, "horizon" => sp.horizon, "tau" => 1,
        "vectors" => sp.role in (:tg, :L5) ? 2 : 1, "coordinate_mode" => "wrapped", "grid_sha256" => "g"^64,
        "status" => status, "attempt" => attempt, "setup_check" => "not applicable", "first_attempt_reasons" => String[], "failure_reasons" => status == "ok" ? String[] : ["failed twice (fixture)"],
        "pilot" => false, "arms" => sp.arms, "session" => 1, "control_stop" => false, "dirty" => false, "created_utc" => "2026-11-01T00:00:00Z")
    sp.role in (:L1, :L5, :tg) && (root["direction"] = sp.direction)
    if sp.kind != ""
        kd, mem = P2_SPEC_INFO[sp.id]
        st_text = "indicator.id = \"$(sp.id)\"\nindicator.kind = \"$kd\"\n"
        merge!(root, Dict{String, Any}("kind" => kd, "spec" => st_text, "spec_sha256" => p3ref_sha256(st_text), "chain" => [kd], "members" => mem))
        kd in ("frequency", "phase_advance", "birkhoff") && merge!(root, Dict{String, Any}("windows" => p.windows, "run_length" => sp.horizon))
        kd == "birkhoff" && (root["observable_names"] = P2_BKROWS)
    end
    return root, ch
end

p3t_flags(o, sp) = (f = String[]; o.tag === :drift && push!(f, "drift"); (o.tag === :U && sp.arms == ["ref"] && sp.role === :L5) && push!(f, "diagnostics exceeded"); f)

# ---------------------------------------------------------------- the records of the items at the horizons of the runner
# the first recorded step not below x and the recorded step nearest to x (a tie goes to the larger step), written here from E26
p3t_first(t, x) = t[searchsortedfirst(t, x)]
function p3t_near(t, x)
    i = searchsortedfirst(t, x)
    i > length(t) && return t[end]
    i == 1 && return t[1]
    return x - t[i - 1] < t[i] - x ? t[i - 1] : t[i]
end
# the tangent items on a grid of any horizon: the encoding of p2_tg (an_p2_study.jl) with n(b_k) the first recorded step not below b_k and n' the recorded step nearest to n / 10; on
# the grid of N_1 this is p2_tg, on the grid of N_c (n = ..., 5632, 10240, 20480) it is the item of the arm ind (E75)
function p3t_tg(id, t, V)
    M = repeat([5.0 5.0 2.0 2.0 2.0 2.0 2.0 2.0 0.1 -1.0], length(t))
    for k in 0:9
        n = p3t_first(t, P2_BUDGETS[k + 1]); np = p3t_near(t, n / 10)
        i, ip = p2_ix(t, n), p2_ix(t, np)
        vs, vf, vm = V[:SALI][k + 1], V[:FLI][k + 1], V[:MEGNO][k + 1]
        M[ip, 10], M[i, 10] = vs === :undecided ? (0.5, -3.0) : (-1.0, -p2_rv(vs, id, k))
        M[i, 9] = vs === :chaotic ? 1e-10 : (vs === :regular && id % 3 == 0) ? 1e-9 : 0.1
        M[ip, 1], M[i, 1] = vf === :undecided ? (0.0, 3.0) : vf === :chaotic ? (5.0, 5.0 * p2_rv(vf, id, k)) : (2.0, 2.0 * p2_rv(vf, id, k))
        M[ip, 5], M[i, 5] = vm === :undecided ? (0.0, 2.0) : vm === :chaotic ? (5.0, 5.0 * p2_rv(vm, id, k)) : (2.0 / p2_rv(vm, id, k), 2.0)
    end
    return M
end

# an item record: the p2 record, except the tangent item (any grid) and the L5 item of the arm chk (own horizon N_B, the values of P3T_CHK_OFF)
function p3t_record(sp, info, design, p; flags = String[])
    if sp.role === :tg
        id, t = info.id, collect(sp.grid)
        Vd = Dict{Symbol, Vector{Symbol}}(X => copy(design.V[id][X]) for X in keys(design.V[id]))
        state, wind = p2_points(info, t, design)
        at = Dict{String, Any}("system" => info.system, "params" => info.params, "ic" => info.ic, "class_label" => info.cls, "seed" => string(1000 + id), "item_hash" => "i"^64,
            "design_class" => info.cls, "flags" => String.(flags), "steps_recorded" => length(t))
        n = length(t)
        diag = Dict{String, Vector{Float64}}(nm => fill(1e-16, n) for nm in ("det", "symp")); thr = Dict{String, Float64}(nm => 1e-12 for nm in ("det", "symp")); fe = Dict{String, Int}(nm => 0 for nm in ("det", "symp"))
        return AX.OrbitRecord(id, at, Int64.(t), state, wind, Dict{String, Matrix{Float64}}(), diag, thr, fe, Dict(sp.id => Float64[sp.horizon, 1, 0.01, 0.1, 1e5]), "",
                              Dict{String, NamedTuple}("tg" => AX._ind_nt(series = p3t_tg(id, t, Vd), columns = P2_TGCOLS)))
    elseif sp.arms == ["chk"]
        d2 = merge(design, (; rho = Dict{Int, Symbol}(k => :C for k in keys(design.rho))))          # no R shape of the L5 item of the p2 fixture (the grid has no step N_1 / 10)
        rec = p2_orbit_record(merge(sp, (; arms = ["x"])), info, d2, p; flags)
        off = get(P3T_CHK_OFF, info.id, (0.0, 0.0))[sp.precision == "Double64" ? 1 : 2]
        rec.lyap["L5"][:, 1] .+= off; rec.lyap["L5"][:, 2] .-= off
        return rec
    end
    return p2_orbit_record(sp, info, design, p; flags)
end

# the cells of a set of orbits: (cell, ids) in the order of round, K index, class
function p3t_cells_of(orbits, ids)
    cells = Dict{String, Vector{Int}}()
    for i in ids; push!(get!(cells, p3t_cell(orbits[i]), Int[]), i); end
    key(c) = (m = match(r"^k(\d+):(\w+):r(\d+)$", c); (parse(Int, m.captures[3]), parse(Int, m.captures[1]), m.captures[2]))
    return [(c, cells[c]) for c in sort(collect(keys(cells)); by = key)]
end

const P3T_NS_PER_STEP = 250.0                       # the thread-ns per orbit-step of the ok chunks of the toy: the cost of L1 of E78 is 250.0

# the whole campaign: dir, panel, store (one chunk per cell and item), manifest; returns the table of chunks (arm, round, hash, cell, status, ids)
function p3t_build(tag::AbstractString; q = p3t_q(), h = nothing)
    dir = an_path("p3t-$tag")
    rm(dir; recursive = true, force = true); mkpath(dir)
    orbits = p3t_orbits()
    des = p3t_design()
    panel = p3t_write_panel(joinpath(dir, "$P3_SLUG-test-panel.toml"), orbits, h === nothing ? "0"^64 : h)
    ph = p3ref_sha256(read(panel, String))
    p = p2_params(; panel_sha256 = ph)
    store = joinpath(dir, "store"); mkpath(store)
    chunks = Dict{String, Any}[]; table = NamedTuple[]
    ord = 0
    h0 = h === nothing ? "0"^64 : h
    function emit(sp, ids, cell, design; status = "ok", rnd = 0, flagfn = (id, sp) -> String[])
        ord += 1
        recs = AX.OrbitRecord[]
        status == "ok" && for id in ids
            o = orbits[id]
            push!(recs, p3t_record(sp, p3t_info(o), design, p; flags = flagfn(id, sp)))
        end
        root, ch = p3t_root(sp, ord, p, ph, h0, cell; status, attempt = status == "ok" ? 1 : 2)
        AX.write_chunk(joinpath(store, "chunk-" * lpad(string(ord), 5, '0') * ".h5"), root, recs)
        steps = length(recs) * sp.horizon
        # a failed chunk entry carries the steps and the seconds of its attempts (large on purpose: a cost that includes it is wrong, E78)
        ent = Dict{String, Any}("ordinal" => ord, "key" => root["chunk_key"], "hash" => ch, "arms" => sp.arms, "status" => status, "attempt" => status == "ok" ? 1 : 2, "orbits" => length(recs),
            "backend" => "cpu", "steps" => status == "ok" ? steps : 3 * sp.horizon, "wall_seconds" => 1.0, "thread_seconds" => status == "ok" ? steps * P3T_NS_PER_STEP * 1e-9 : 7.0,
            "fallback_reason" => "", "first_attempt_reasons" => String[], "setup_check" => "not applicable", "control_stop" => false, "control_info" => Any[],
            "flags" => Dict{String, Any}("diagnostics exceeded" => Int[], "drift" => [i for i in ids if orbits[i].tag === :drift && status == "ok"], "excluded" => status == "ok" ? Int[] : collect(ids), "breakdown_tau" => Int[]))
        sp.kind != "" && (ent["indicator_flags"] = Dict{String, Any}(k => v for (k, v) in p2_flag_counts(recs)))
        push!(chunks, ent)
        push!(table, (; ordinal = ord, arm = sp.arms[1], role = sp.role, spec = sp.id, precision = sp.precision, round = rnd, cell, hash = ch, status, ids = collect(ids), horizon = sp.horizon))
    end
    flags(id, sp) = p3t_flags(orbits[id], sp)
    allids = [o.id for o in orbits]
    # per arm: the orbits it holds, the failed cells (a cell of ref fails with its excluded orbit, a cell of ind and ind_pa1 with its orbit without verdict data), the design of its items
    for (kind, ids, design) in ((:ref, allids, des.f), (:ind, allids, des.f), (:ind_pa1, allids, des.f), (:d64, P3T_D, des.d), (:d64_pa1, P3T_D, des.d), (:chk, P3T_DB, des.f))
        for (cell, cids) in p3t_cells_of(orbits, ids)
            failed = (kind === :ref && any(i -> orbits[i].tag === :excl, cids)) || (kind in (:ind, :ind_pa1) && any(i -> orbits[i].tag === :nov, cids))
            rnd = orbits[first(cids)].round
            for sp in p3t_specs(kind, p)
                emit(sp, cids, cell, design; status = failed ? "failed" : "ok", rnd, flagfn = flags)
            end
        end
    end
    manifest = joinpath(dir, "$P3_SLUG-test-manifest.json")
    p3t_write_manifest(manifest, p, chunks, h0, ph, p3t_arms(chunks, p))
    return (; dir, store, panel, manifest, panel_hash = ph, orbits, design = des, chunks = table, q, p, files = length(table))
end

# the arm blocks of the manifest (the keys of the runner: budget_ns, cost_ns_per_orbit_step, status, steps, thread_seconds, wall_seconds; `specs` on the indicator arms only, E78)
function p3t_arms(chunks, p)
    arms = Dict{String, Any}()
    for a in P3T_ARMS
        es = [c for c in chunks if a in c["arms"]]
        steps = sum(c["steps"] for c in es); thr = sum(c["thread_seconds"] for c in es)
        arms[a] = Dict{String, Any}("status" => "complete", "steps" => steps, "cost_ns_per_orbit_step" => 1e9 * thr / max(steps, 1), "budget_ns" => 415.2, "thread_seconds" => thr, "wall_seconds" => thr / 4)
    end
    for sp in p3t_allspecs(p)
        sp.kind == "" && continue
        # the specs of the Double64 arms differ from those of the Float64 arms (a reader of the wrong arm is caught, P119)
        mult = startswith(sp.arms[1], "d64") ? 10.0 : 1.0
        haskey(P2_COST, sp.id) && (get!(arms[sp.arms[1]], "specs", Dict{String, Any}())[sp.id] = Dict{String, Any}("budget_line_ns" => P2_COST[sp.id][1], "measured_ns_per_orbit_step" => mult * P2_COST[sp.id][2]))
    end
    return arms
end

function p3t_write_manifest(path, p, chunks, h0, ph, arms; status = nothing)
    failed = any(c -> c["status"] == "failed", chunks)
    AX._write_json(path, Dict{String, Any}("study" => P3_SLUG, "schema_version" => 2, "status" => status === nothing ? (failed ? "complete with failed chunks" : "complete") : status, "plan_sha256" => p.plan_sha256,
        "plan_counts" => Dict("items" => sum(length(c["orbits"]) for c in chunks), "chunks" => length(chunks), "steps_before" => 1, "steps_after" => 1),
        "git" => Dict("sha" => "b"^40, "clean" => true), "sessions" => Any[Dict("clean" => true, "git_sha" => "b"^40, "index" => 1, "threads_used" => 4, "wall_seconds" => 10.0)],
        "prereg" => Dict("header" => h0, "sha256_start" => h0, "sha256_end" => h0), "panel" => Dict("sha256_start" => ph, "sha256_end" => ph, "build_hash" => "e"^64),
        "overrides" => Dict{String, Any}(), "control_stops" => Any[], "reasons" => String[], "chunks" => chunks, "arms" => arms, "indicators" => Dict{String, Any}()))
end

# ---------------------------------------------------------------- the inputs of the test side: record, 71-key table, test-side table, preregistration (one set per target pair)
function p3t_model()
    r(p, s) = (p = p, s = s)
    comp = Dict{Symbol, Any}(X => (T = r(0.6, 0.5), RW = r(0.4, 0.3)) for X in P2_X7)
    return (rule = Dict(:T => (p = 0.1, d = 0.9), :RW => (p = 0.05, d = 0.9)), comp = comp, ts = P2_X7, w2 = P2_X7)
end
function p3t_record_nt(q, t, h0, sha71)
    th = (V = String.(P3T_THETA.V), qC = P3T_THETA.qC, qR = P3T_THETA.qR, w = P3T_THETA.w, v = P3T_THETA.v)
    row = (j = 0, rho = 1 // 1, nT = 5, nRW = 6, nW = t.target_w2, allowed = true, successes = 12, runs = q.pp_runs, power = t.power_hb_planned, se = 0.05)
    return (program = (generator = "xoshiro256++", quantile = "AS241", git_commit = "b"^40, git_clean = true),
            fingerprints = (preregistration = h0, p3_parameters = sha71, git_commit = "b"^40, git_clean = true),
            parameters = NamedTuple{Tuple(P3_PARAM_KEYS)}(Tuple(getfield(q, k) for k in P3_PARAM_KEYS)),
            outcome = (selection = "selected", optimism_skipped = 0, pp_hb = "ok", pp_hc = "target reached"),
            relabelling = (dropped = (total = 0, by_reason = Dict{String, Int}()),),
            sets = (cal = (pooled = 36, panel = 12, grid = 24), S2 = (pooled = 8, panel = 2, grid = 6), T_s = (pooled = 6, panel = 2, grid = 4), W2 = (pooled = 13, panel = 5, grid = 8)),
            theta_space = (size = 55,), selection = (feasible = 1, tol = 1, index = 11, star = (index = 11, V = th.V, qC = th.qC, qR = th.qR, w = th.w, v = th.v), J_min = 1 // 20, theta = th),
            pp_hb = (rates = Dict{String, Any}(), model = p3t_model(), path = [row], targets = (outcome = "ok", cap = t.cap_flag, j = 0, rho_star = 1 // 1, n_T = t.target_t_s, n_W = t.target_w2)),
            pp_hc = (auc_i_used = 0.9, auc_p = 0.55, deliverable_expected = t.target_f, deliverable_worst = t.target_f + 10, power_expected = (nf = t.target_f, successes = 33, runs = q.pp_runs, power = t.power_hc_planned, se = 0.05), target = t.target_f,
                     target_power = 0.8, underpowered = false, power_cap = false), r_ab = (), readings = String[])
end

# writes prereg, 71-key table, record, test-side table of one pair of targets into a new folder; returns the paths and the values. `st` is the D93 study (the p2 and P1 inputs)
function p3t_inputs(S, st, tag::AbstractString; target_t_s = 5, target_w2 = 11, target_f = 6, q = S.q, values = (;), a1 = true, record_edit = identity)
    dir = an_path("p3t-in-$tag"); rm(dir; recursive = true, force = true); mkpath(dir)
    q = merge(q, (; p2_manifest_sha256 = p3ref_sha256_file(st.s2.manifest)))
    t0 = p3t_values(; target_t_s, target_w2, target_f, values...)
    text0, h0 = p3t_prereg_text(q, t0; a1 = false)
    pre0 = joinpath(dir, "pre-a1.md"); write(pre0, text0)
    params = p3_write_params(joinpath(dir, "$P3_SLUG-calibration-parameters.toml"), q)
    sha71 = p3ref_sha256_file(params)
    rec = joinpath(dir, "$P3_SLUG-calibration-record.json")
    write(rec, AN.p3_json(record_edit(p3t_record_nt(q, t0, h0, sha71))))
    t = merge(t0, (; calibration_record_sha256 = p3ref_sha256_file(rec), p3_table_sha256 = sha71))
    text, h = p3t_prereg_text(q, t; a1)
    prereg = joinpath(dir, "$P3_SLUG-preregistration.md"); write(prereg, text)
    rm(pre0)
    tparams = p3t_write_values(joinpath(dir, "$P3_SLUG-test-parameters.toml"), t, q)
    return (; dir, prereg, params, tparams, record = rec, q, t, prereg_hash = h, qt = merge(q, t))
end
p3ref_sha256_file(f) = p3ref_sha256(String(read(f)))

# ---------------------------------------------------------------- expected values from the design (the independent reference)
# counts after round r (E62 b); the text of the entry, the chunk hashes of the ref arm up to round r
p3t_label(o) = o.tag === :U ? :U : P3_RHO_NEW[o.tmpl]
function p3t_counts(S, qt, r)
    D = S.design.f
    os = [o for o in S.orbits if o.round <= r && o.tag !== :gate]
    exc = count(o -> o.tag === :excl, os); dr = count(o -> o.tag === :drift, os)
    ok = [o for o in os if o.tag !== :excl && o.tag !== :drift && p3t_label(o) === :C]
    fs = p3tref_rat(qt.f_s); lw = p3tref_rat(qt.lam_w2)
    ts(o) = p3tref_rat(D.lns[o.id]) <= fs * p3tref_rat(D.lref[o.id])
    w2(o) = ts(o) || p3tref_rat(D.lref[o.id]) <= lw
    return Dict("drift" => dr, "excluded" => exc, "orbits" => length(os), "t_s" => count(ts, ok), "target_t_s" => qt.target_t_s, "target_w2" => qt.target_w2, "w2" => count(w2, ok))
end
p3t_ref_hashes(S, r) = [c.hash for c in S.chunks if c.arm == "ref" && c.round <= r]
p3t_digest(S, r) = p3ref_sha256(join(sort(p3t_ref_hashes(S, r)) .* "\n"))
# the entries of the log by the reference writer, up to the close (or r_max); returns the list of texts and the decisions
function p3t_log_entries(S, I; upto = P3T_ROUNDS, prereg = I.prereg_hash, panel = S.panel_hash)
    ents = String[]; decs = String[]; prev = "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855"
    for r in 1:upto
        c = p3t_counts(S, I.qt, r)
        dec = p3tref_decide(c["t_s"], c["w2"], I.qt.target_t_s, I.qt.target_w2, r, Int(I.qt.r_max))
        e = p3tref_entry(; round = r, decision = dec, counts = c, panel, prereg, digest = p3t_digest(S, r), prev)
        push!(ents, e); push!(decs, dec); prev = p3ref_sha256(e)
        dec == "close" && break
    end
    return ents, decs
end

# ---------------------------------------------------------------- an instance of the campaign for one set of inputs: the store with the preregistration hash of the inputs patched into
# the root attributes, the manifest rewritten; `rounds` limits the ref chunks (the state of a campaign after that round), `closure` keeps or drops the closure arms, `clo_rounds` limits the
# closure chunks to the cells of the rounds up to it (the runner releases the closure arms only up to the closing round)
function p3t_instance(S, I, tag::AbstractString; rounds::Int = P3T_ROUNDS, closure::Bool = true, panel_edit = nothing, clo_rounds::Int = P3T_ROUNDS)
    dir = an_path("p3t-inst-$tag"); rm(dir; recursive = true, force = true); mkpath(dir)
    store = joinpath(dir, "store"); mkpath(store)
    keep = Set{Int}()
    for c in S.chunks
        ((c.arm == "ref" && c.round <= rounds) || (c.arm != "ref" && closure && c.round <= clo_rounds)) && push!(keep, c.ordinal)
    end
    panel = joinpath(dir, "panel.toml")
    if panel_edit === nothing
        cp(S.panel, panel)
    else
        d = TOML.parsefile(S.panel); panel_edit(d)
        open(io -> TOML.print(io, d), panel, "w")
    end
    ph = p3ref_sha256(read(panel, String))
    for f in readdir(S.store; join = true)
        m = match(r"chunk-(\d{5})\.h5$", f)
        (m === nothing || !(parse(Int, m.captures[1]) in keep)) && continue
        dst = joinpath(store, basename(f)); cp(f, dst)
        AX.HDF5.h5open(dst, "r+") do h
            AX.HDF5.delete_attribute(h, "prereg_sha256"); AX.HDF5.write_attribute(h, "prereg_sha256", I.prereg_hash)
            AX.HDF5.delete_attribute(h, "panel_sha256"); AX.HDF5.write_attribute(h, "panel_sha256", ph)
        end
    end
    man = AX.read_json(S.manifest)
    man["prereg"] = Dict("header" => I.prereg_hash, "sha256_start" => I.prereg_hash, "sha256_end" => I.prereg_hash)
    man["panel"] = Dict("sha256_start" => ph, "sha256_end" => ph, "build_hash" => "e"^64)
    man["chunks"] = [c for c in man["chunks"] if c["ordinal"] in keep]
    man["status"] = any(c -> c["status"] == "failed", man["chunks"]) ? "complete with failed chunks" : "complete"
    manifest = joinpath(dir, "manifest.json"); AX._write_json(manifest, man)
    return (; dir, store, manifest, panel, panel_hash = ph, log = joinpath(dir, "round-log.txt"), out = joinpath(dir, "out"), keep)
end
p3t_write_log(inst, ents) = (write(inst.log, join(ents)); inst.log)

# the expected quantities of the analysis, from the design (the independent reference): sets, strata, outputs, HA, HB, HC
function p3t_expected(S, qt; noD = Int[])
    D, Dd = S.design.f, S.design.d
    th = (V = P3T_THETA.V, qC = P3T_THETA.qC, qR = P3T_THETA.qR, w = P3T_THETA.w, v = P3T_THETA.v)
    kP = p3ref_kP(qt)
    fs, lw, ls = p3tref_rat(qt.f_s), p3tref_rat(qt.lam_w2), p3tref_rat(qt.lam_s2)
    row(o) = (rho = p3t_label(o); C = rho === :C; lr = p3tref_rat(D.lref[o.id]); ln = p3tref_rat(D.lns[o.id]); Ts = C && ln <= fs * lr;
              (id = o.id, rho = rho, out = p3ref_out(D.V[o.id], kP, th), W2 = C && (lr <= lw || Ts), Ts = Ts, S2 = C && lr >= ls))
    drops = NamedTuple[]
    E = NamedTuple[]
    for o in S.orbits
        o.tag === :gate && continue
        reason = o.tag === :excl ? "excluded" : o.tag === :drift ? "drift" : o.tag === :U ? "unresolved" : o.tag === :nov ? "no verdict data" : ""
        isempty(reason) ? push!(E, row(o)) : push!(drops, (id = o.id, reason = reason, round = o.round, kidx = o.kidx))
    end
    gate = [p3ref_out(D.V[o.id], kP, th) for o in S.orbits if o.tag === :gate]
    ha = p3tref_ha(E, gate, qt)
    byid = Dict(r.id => r for r in E)
    M = p3ref_matched(qt, th.V)
    names = Symbol.(qt.voter_bits)
    pos(X) = X === :LE ? findfirst(==(M.le_read), qt.budgets) : findfirst(==(M.X[X].b), qt.budgets)
    capped = [X for X in names if M.X[X].capped]
    cmp = Dict{Tuple{Symbol, Symbol}, Any}(); cov = Dict{Symbol, Any}(); sizes = (Ts = 0, W2 = 0)
    for (A, key, Xs) in ((:Ts, :Ts, names), (:W2, :W2, vcat(:LE, names)))
        set = sort([r for r in E if getfield(r, A)]; by = r -> r.id)
        cov[key] = P3T_R(count(r -> r.out !== :unresolved, set)) // length(set)
        for X in Xs
            cmp[(A, X)] = (rule = [r.out === :regular_supported for r in set], x = [D.V[r.id][X][pos(X)] === :regular for r in set])
        end
    end
    sizes = (Ts = count(r -> r.Ts, E), W2 = count(r -> r.W2, E))
    hb = p3tref_hb(cmp, cov, sizes, qt; capped)
    hcrows = NamedTuple[]
    for r in E
        (r.id in P3T_D && !(r.id in noD)) || continue
        nC = count(X -> D.V[r.id][X][kP + 1] === :chaotic, th.V); nR = count(X -> D.V[r.id][X][kP + 1] === :regular, th.V)
        sP = P3T_R(count(X -> (a = D.V[r.id][X][kP + 1]; b = Dd.V[r.id][X][kP + 1]; a !== :undecided && b !== :undecided && a != b), th.V)) // length(th.V)
        push!(hcrows, (id = r.id, rho = r.rho, nC = nC, nR = nR, sP = sP))
    end
    hc = p3tref_hc(hcrows, qt)
    vfun(id, Y) = D.V[id][Y]
    return (; E, drops, gate, ha, hb, hc, cmp, cov, sizes, M, th, kP, hcrows, byid, vfun, D = [r.id for r in hcrows])
end

# ---------------------------------------------------------------- shared objects (built once) and small utilities
p3t_files(store) = sort([joinpath(store, f) for f in readdir(store) if occursin(r"^chunk-\d{5}\.h5$", f)])
const P3T_CACHE = Dict{String, Any}()
p3t_cached() = get!(() -> p3t_build("main"), P3T_CACHE, "S")
p3t_cached_inputs() = get!(() -> p3t_inputs(p3t_cached(), p3_cached("d95"), "main"), P3T_CACHE, "I")

# ---------------------------------------------------------------- calls of the two commands, clones, fingerprints, alterations of an instance
# keyword arguments of count_round and analyze (the p2 inputs are those of the D93 study `st`)
p3t_cr_args(I, st, inst; kw...) = (; prereg = I.prereg, params = I.params, test_params = I.tparams, store = inst.store, panel = inst.panel, round_log = inst.log, p2_manifest = st.s2.manifest,
    p2_prereg = st.s2.prereg, p2_params = st.s2.params, calibration_record = I.record, kw...)
p3t_an_args(I, st, inst; kw...) = (; p3t_cr_args(I, st, inst)..., manifest = inst.manifest, p1_prereg = st.p1.prereg, p1_params = st.p1.params, out = inst.out, kw...)
p3t_cli_cr(I, st, inst; extra = String[]) = vcat(["count-round", I.prereg, I.params, I.tparams, "--store", inst.store, "--panel", inst.panel, "--round-log", inst.log, "--p2-manifest", st.s2.manifest,
    "--p2-preregistration", st.s2.prereg, "--p2-parameters", st.s2.params, "--calibration-record", I.record], extra)
p3t_cli_an(I, st, inst; extra = String[]) = vcat(["analyze", I.prereg, I.params, I.tparams, "--store", inst.store, "--panel", inst.panel, "--round-log", inst.log, "--manifest", inst.manifest,
    "--p2-manifest", st.s2.manifest, "--p2-preregistration", st.s2.prereg, "--p2-parameters", st.s2.params, "--calibration-record", I.record, "--p1-preregistration", st.p1.prereg,
    "--p1-parameters", st.p1.params, "--out", inst.out], extra)

# a copy of an instance in a new folder (the log, the store, the manifest, the panel); the inputs of the inputs-set I stay where they are
function p3t_clone(inst, tag::AbstractString)
    dir = an_path("p3t-clone-$tag"); rm(dir; recursive = true, force = true)
    sub(s) = replace(s, inst.dir => dir)
    cp(inst.dir, dir)
    rm(joinpath(dir, "out"); recursive = true, force = true)
    return merge(inst, (; dir, store = sub(inst.store), manifest = sub(inst.manifest), panel = sub(inst.panel), log = sub(inst.log), out = joinpath(dir, "out")))
end
# the bytes of all files of the inputs and of an instance (what must not change): a list of (path, hash of the bytes)
p3t_fp(paths) = [(p, bytes2hex(AN_SHA(read(p)))) for p in paths if isfile(p)]
AN_SHA(b) = hex2bytes(an_sha256(b))
function p3t_fp_all(I, st, inst)
    fs = vcat([I.prereg, I.params, I.tparams, I.record, st.s2.manifest, st.s2.prereg, st.s2.params, st.p1.prereg, st.p1.params, inst.manifest, inst.panel, inst.log], p3t_files(inst.store))
    return [(f, p3ref_sha256(String(read(f)))) for f in fs if isfile(f)]
end

# set or replace a root attribute of a chunk file
function p3t_set_attr(file, name, value)
    AX.HDF5.h5open(file, "r+") do h
        AX.HDF5.haskey(AX.HDF5.attributes(h), name) && AX.HDF5.delete_attribute(h, name)
        AX.HDF5.write_attribute(h, name, value)
    end
end
# the first chunk file of the instance that matches the arm, role and cell prefix
p3t_chunk(S, inst; arm, role, cell = "") = joinpath(inst.store, "chunk-" * lpad(string(first(c for c in S.chunks if c.arm == arm && c.role === role && startswith(c.cell, cell) && c.status == "ok").ordinal), 5, '0') * ".h5")
# rewrite a chunk of the instance from a changed design (same root attributes, other items); `drop` removes orbits from the chunk
function p3t_alter_chunk(S, inst, I, ordinal::Integer; design = nothing, drop = Int[], edit = Dict{String, Any}(), sp_edit = identity)
    c = first(c for c in S.chunks if c.ordinal == ordinal)
    file = joinpath(inst.store, "chunk-" * lpad(string(ordinal), 5, '0') * ".h5")
    attrs = AX.read_chunk_attrs(file)
    sp = sp_edit(first(s for s in p3t_allspecs(S.p) if s.arms == [c.arm] && s.id == c.spec && s.precision == c.precision))
    design === nothing && (design = c.arm in ("d64", "d64_pa1") ? S.design.d : S.design.f)
    recs = AX.OrbitRecord[p3t_record(sp, p3t_info(S.orbits[id]), design, S.p; flags = p3t_flags(S.orbits[id], sp)) for id in c.ids if !(id in drop)]
    rm(file)
    AX.write_chunk(file, merge(attrs, edit), recs)
    return file
end

# ---------------------------------------------------------------- derived instances (Amendment D97): failed and healed chunks, the same campaign with other tags
# chunks of an instance rewritten as failed twice: no orbit record, status failed, attempt 2; the manifest entry gets status failed, the steps and seconds of the attempts, the orbits as excluded
function p3t_fail_chunks(S, inst, ordinals)
    man = AX.read_json(inst.manifest)
    for o in ordinals
        c = first(c for c in S.chunks if c.ordinal == o)
        file = joinpath(inst.store, "chunk-" * lpad(string(o), 5, '0') * ".h5")
        attrs = AX.read_chunk_attrs(file)
        merge!(attrs, Dict{String, Any}("status" => "failed", "attempt" => 2, "failure_reasons" => ["failed twice (fixture)"]))
        rm(file); AX.write_chunk(file, attrs, AX.OrbitRecord[])
        e = first(e for e in man["chunks"] if e["ordinal"] == o)
        merge!(e, Dict{String, Any}("status" => "failed", "attempt" => 2, "orbits" => 0, "steps" => 3 * c.horizon, "thread_seconds" => 7.0))
        e["flags"]["excluded"] = collect(c.ids); e["flags"]["drift"] = Int[]
    end
    man["status"] = "complete with failed chunks"
    AX._write_json(inst.manifest, man)
    return inst
end
# the failed chunks of an instance rewritten as ok chunks with their items (the campaign without failures); the manifest status is `complete` if no entry is left failed
function p3t_heal_chunks(S, inst, I, ordinals)
    man = AX.read_json(inst.manifest)
    for o in ordinals
        c = first(c for c in S.chunks if c.ordinal == o)
        p3t_alter_chunk(S, inst, I, o; edit = Dict{String, Any}("status" => "ok", "attempt" => 1, "failure_reasons" => String[]))
        e = first(e for e in man["chunks"] if e["ordinal"] == o)
        steps = length(c.ids) * c.horizon
        merge!(e, Dict{String, Any}("status" => "ok", "attempt" => 1, "orbits" => length(c.ids), "steps" => steps, "thread_seconds" => steps * P3T_NS_PER_STEP * 1e-9))
        e["flags"]["excluded"] = Int[]
    end
    man["status"] = any(e -> e["status"] == "failed", man["chunks"]) ? "complete with failed chunks" : "complete"
    AX._write_json(inst.manifest, man)
    return inst
end
# the same campaign with the tags of some orbits changed (the design of the verdict histories does not depend on the tag, except :U and :drift)
p3t_retag(S, ids, tag) = merge(S, (; orbits = [o.id in ids ? merge(o, (; tag)) : o for o in S.orbits]))
# the ordinals of the chunks of an arm (and a spec, a cell) in the table of the campaign
p3t_ords(S; arm, spec = nothing, cell = nothing, precision = nothing) = [c.ordinal for c in S.chunks if c.arm == arm && (spec === nothing || c.spec == spec) && (cell === nothing || c.cell == cell) && (precision === nothing || c.precision == precision)]
