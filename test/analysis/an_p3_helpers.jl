# Amendment D93, stage 02: helpers of the p3-abstention suite (prefix p3_). Definitions only: nothing here calls Chaoslab.Analysis at load time, so a missing function
# errors the testsets that need it and nothing else. Uses an_helpers.jl and the p2 helpers (an_p2_*.jl) without changing them.
# Spec: workspaces/module-build/stages/01-spec/output/analysis-spec.md, "Amendment D93". API assumed: analysis-test-plan.md, "Amendment D93 (stage 02)", "Contract proposed by the tests".
#
# The synthetic study "p3-abstention" (scaled): the p2 fixture panel and grid (16 panel orbits, 4 controls, 24 grid starts) with its p2 store of 33 chunk files written by the
# real chunk writer, a P1 reference store whose horizon N_1 equals the p2 horizon (80000, so that the panel relabelling C0 finds its L5 items), a fixture preregistration with
# the SELF hash rule, and a parameter table of 71 keys. Budgets 40 x 2^k, b_P = 2560, N_c = 20480 (the p2 mirror of the fixture). The costs are the frozen costs of the
# preregistration (the fixture keeps b_P / N_c = 1 / 8, so that the cost cap F1 removes the same elements as in the real study: 690 pass). The thresholds of calibration are
# scaled to the 36 orbits of Cal (the fixture owns its own numbers: it is not the study).

p3_seed(letter::Char, sub::Integer = 0) = AN_SEED + 34000 + 50 * (letter - 'a' + 1) + sub      # D93<letter> draws use 34000 + 50 n + sub, sub < 50
const P3_SLUG = "p3-abstention"

# Stage 03, owner decision 2 of the stage 02 checkpoint (D38): the package suite carries 30 of the 52 refusal cases of D93c (the criterion asks for at least 30, one for each rule a to k and each
# sub-cause); verify/analysis-p3.jl sets P3_FULL_RUN = 1 and runs all 52 with the same assertions (test plan, "Rework at stage 03").
const P3_FULL = get(ENV, "P3_FULL_RUN", "") == "1"
const P3_CASE_SAMPLE = ["slug-absent", "slug-invalid", "hash-row-absent", "hash-row-repeated", "hash-wrong", "param-missing", "param-unknown", "param-quote", "param-range", "param-mirror", "param-voters",
    "p2-param-missing", "p2-manifest-hash", "p2-study", "p2-prereg-header", "p2-status", "p2-chunk-failed", "p2-pilot", "p2-nan-series", "p2-drift-unflagged", "p1-prereg", "p1-status", "p1-chunk-hash",
    "p1-study", "id-grid-ic", "id-p1-ic", "item-missing-cal", "item-missing-gate", "item-listed-exists", "step-missing"]
const P3_I7 = ["LE", "SALI", "FLI", "MEGNO", "E1", "E2", "WBA"]
const P3_PARAM_KEYS = [:budgets, :b_primary, :N_c, :N_1, :n_s, :f_s, :lam_w2, :lam_s2, :min_set, :gate_size,
    :indicators_i7, :voter_bits, :w_range, :v_range, :min_voters,
    :cost_l1, :cost_tg, :cost_pa1, :cost_bk1, :cost_fq1, :ratio_d64_tg, :ratio_d64_fq1, :ratio_d64_pa1, :ratio_d64_bk1, :spec_of,
    :e_cal, :e_max, :a_cal, :a_max, :s_cal, :delta_a, :alpha_ci, :g_max, :c_min, :alpha_hb, :rho_hb, :n_t_min, :n_w_min, :p_hb, :p_hc, :delta_c, :delta_tol,
    :b_opt, :rule_rate_floor, :rule_rate_factor, :auc_shrink, :auc_p, :nonfail_per_fail, :pp_runs, :test_seed, :r_max, :rho_path_start, :rho_path_second, :y_t, :y_w,
    :nf_start, :nf_second, :ceil_step, :z_one_sided, :shared_default, :d_expected_total, :d_worst_total,
    :screening_points_per_axis, :d_frac_high, :d_frac_low, :d_b_size, :safety_float, :safety_dd, :wall_cap_hours, :bootstrap_resamples,
    :p2_manifest_sha256]
const P3_INT_KEYS = [:b_primary, :N_c, :N_1, :n_s, :min_set, :gate_size, :min_voters, :n_t_min, :n_w_min, :b_opt, :nonfail_per_fail, :pp_runs, :test_seed, :y_t, :y_w, :nf_start,
                     :nf_second, :ceil_step, :d_expected_total, :d_worst_total, :screening_points_per_axis, :d_b_size, :bootstrap_resamples, :g_max]
const P3_INTVEC_KEYS = [:budgets, :w_range, :v_range]
const P3_STRVEC_KEYS = [:indicators_i7, :voter_bits]
const P3_CARRIED = [:screening_points_per_axis, :d_frac_high, :d_frac_low, :d_b_size, :safety_float, :safety_dd, :wall_cap_hours, :bootstrap_resamples]

# the value table of the fixture (the mirror keys follow the p2 fixture table of an_p2_helpers.jl)
p3_params(; kw...) = merge((; budgets = P2_BUDGETS, b_primary = 2560, N_c = P2_NC, N_1 = P2_N1, n_s = P2_NS, f_s = 0.2, lam_w2 = 0.1, lam_s2 = 1.0, min_set = 5, gate_size = 2,
    indicators_i7 = copy(P3_I7), voter_bits = P3_I7[2:end], w_range = [0, 1, 2], v_range = [0, 1], min_voters = 3,
    cost_l1 = 207.6, cost_tg = 2658.0, cost_pa1 = 1810.5, cost_bk1 = 3000.6, cost_fq1 = 16850.4, ratio_d64_tg = 10.5, ratio_d64_fq1 = 1.38, ratio_d64_pa1 = 7.28, ratio_d64_bk1 = 6.56,
    spec_of = Dict("LE" => "l1", "SALI" => "tg", "FLI" => "tg", "MEGNO" => "tg", "E1" => "fq1", "E2" => "pa1", "WBA" => "bk1"),
    e_cal = 0.2, e_max = 0.3, a_cal = 0.3, a_max = 0.4, s_cal = 0.15, delta_a = 0.1, alpha_ci = 0.05, g_max = 1, c_min = 0.5, alpha_hb = 0.05, rho_hb = 0.5, n_t_min = 100,
    n_w_min = 300, p_hb = 0.8, p_hc = 0.8, delta_c = 0.1, delta_tol = 0.05,
    b_opt = 6, rule_rate_floor = 0.01, rule_rate_factor = 2.0, auc_shrink = 0.5, auc_p = 0.55, nonfail_per_fail = 20, pp_runs = 60, test_seed = 20261030, r_max = 5.0,
    rho_path_start = 1.0, rho_path_second = 1.25, y_t = 57, y_w = 247, nf_start = 20, nf_second = 25, ceil_step = 5, z_one_sided = 1.645, shared_default = 0.5,
    d_expected_total = 1048, d_worst_total = 1309,
    screening_points_per_axis = 512, d_frac_high = 0.3, d_frac_low = 0.05, d_b_size = 100, safety_float = 2.0, safety_dd = 3.0, wall_cap_hours = 12.0, bootstrap_resamples = 10000,
    p2_manifest_sha256 = "0"^64), NamedTuple(kw))

# section headings of the fixture preregistration (the real ones) and, per key, the section and the verbatim quote that contains the value text
const P3_SEC = (REG = "Thresholds register (single home; all set at the checkpoint, 2026-10-30)",
    SETS = "Sets (this study's own definitions; the parameters are in the register, marked as p2 mirrors)", CAL = "Calibration data (read only)", RULE = "The rule R",
    COST = "Cost model and matched budget", PP = "Test draw size and power calculation (procedure PP, frozen in Amendment A1)", CALP = "Calibration procedure (deterministic; its output is frozen by Amendment A1)",
    BUDGET = "Budget", METH = "Methods")
const P3_SEC_ORDER = [:CAL, :SETS, :METH, :RULE, :COST, :CALP, :PP, :BUDGET, :REG]
const P3_KEY_SEC = Dict{Symbol, Symbol}(
    :budgets => :REG, :b_primary => :REG, :N_c => :REG, :N_1 => :REG, :n_s => :SETS, :f_s => :SETS, :lam_w2 => :SETS, :lam_s2 => :SETS, :min_set => :REG, :gate_size => :SETS,
    :indicators_i7 => :RULE, :voter_bits => :RULE, :w_range => :REG, :v_range => :REG, :min_voters => :RULE,
    :cost_l1 => :COST, :cost_tg => :COST, :cost_pa1 => :COST, :cost_bk1 => :COST, :cost_fq1 => :COST, :ratio_d64_tg => :COST, :ratio_d64_fq1 => :COST, :ratio_d64_pa1 => :COST,
    :ratio_d64_bk1 => :COST, :spec_of => :COST, :p2_manifest_sha256 => :CAL, :y_t => :PP, :y_w => :PP, :nf_start => :PP, :nf_second => :PP, :ceil_step => :PP, :z_one_sided => :PP,
    :shared_default => :PP, :d_expected_total => :BUDGET, :d_worst_total => :BUDGET, :rho_path_start => :PP, :rho_path_second => :PP)
p3_sec(k) = getfield(P3_SEC, get(P3_KEY_SEC, k, :REG))

# the text of a value as the section writes it (A2): the field `text`; it differs from Julia's string(value) for several keys (0.10, 0.80, 10^4, ...)
function p3_texts(p; texts = Dict{Symbol, String}())
    T = Dict{Symbol, String}()
    for k in P3_PARAM_KEYS
        v = getfield(p, k)
        T[k] = k === :spec_of ? join(["$x=$(v[x])" for x in p.indicators_i7 if haskey(v, x)], ", ") : v isa AbstractVector ? join(v, ", ") : string(v)
    end
    two(x) = !isfinite(x) ? string(x) : (u = string(round(x; digits = 2)); occursin('.', u) && length(split(u, '.')[2]) < 2 ? u * "0" : u)
    for k in (:a_cal, :a_max, :p_hb, :p_hc, :delta_a, :delta_c, :delta_tol, :e_cal, :e_max, :s_cal, :alpha_ci, :alpha_hb, :rule_rate_floor, :d_frac_low)
        T[k] = two(getfield(p, k))
    end
    T[:bootstrap_resamples] = (e = round(Int, log10(p.bootstrap_resamples)); p.bootstrap_resamples == 10^e ? "10^$e" : string(p.bootstrap_resamples))
    return merge(T, texts)
end

function p3_quotes(p; texts = Dict{Symbol, String}())
    T = p3_texts(p; texts)
    q = Dict{Symbol, Tuple{String, String, String}}()
    for k in P3_PARAM_KEYS
        what = k === :spec_of ? "Spec of each indicator (indicator=spec): $(T[k])" : k === :budgets ? "Budget list B = {b_k, k = 0 to 9} = $(T[k])" :
               k === :p2_manifest_sha256 ? "SHA-256 of the p2 manifest $(T[k])" : "$(k) = $(T[k])"
        q[k] = (p3_sec(k), what, T[k])
    end
    return q
end

# TOML text of the parameter table (one table per key: value, section, quote and text). edit[k] = Dict(field => value or nothing to remove)
function p3_write_params(path::AbstractString, p::NamedTuple; drop = Symbol[], extra = Dict{String, Any}(), edit = Dict{Symbol, Any}(), texts = Dict{Symbol, String}())
    q = p3_quotes(p; texts)
    d = Dict{String, Any}()
    for k in P3_PARAM_KEYS
        k in drop && continue
        v = getfield(p, k)
        (k in P3_INT_KEYS && v isa Real) && (v = Int(v))
        (k in P3_INTVEC_KEYS) && (v = Int.(v))
        e = Dict{String, Any}("value" => v, "section" => q[k][1], "quote" => q[k][2], "text" => q[k][3])
        haskey(edit, k) && for (kk, vv) in edit[k]; vv === nothing ? delete!(e, kk) : (e[kk] = vv); end
        d[string(k)] = e
    end
    merge!(d, extra)
    open(io -> TOML.print(io, d), path, "w")
    return path
end

# the fixture preregistration text with the SELF hash rule; returns (text, header hash)
function p3_prereg_text(p; slug = P3_SLUG, hashrow = true, drop_sections = Symbol[], extra = "", texts = Dict{Symbol, String}())
    q = p3_quotes(p; texts)
    body = IOBuffer()
    println(body, "# Preregistration: ", slug, "\n\n## Header\n")
    println(body, "| Field | Value |\n|-------|-------|\n| Study slug | ", slug, " |")
    hashrow && println(body, "| SHA-256 of this file (current version) | ", "0"^64, " (rule as p2: fixture) |")
    println(body, "| Amendment log | none yet |\n\nFixture. Not a study.\n")
    for sk in P3_SEC_ORDER
        sk in drop_sections && continue
        println(body, "## ", getfield(P3_SEC, sk), "\n")
        for k in P3_PARAM_KEYS
            q[k][1] == getfield(P3_SEC, sk) && println(body, "- ", q[k][2])
        end
        println(body)
    end
    println(body, extra)
    text = String(take!(body))
    prep(t) = replace(replace(t, "\r\n" => "\n"), r"(?m)^(\| SHA-256 of this file \(current version\) \|)[^\n|]*\|" => s"\1 SELF |")
    h = an_sha256(prep(text))
    hashrow && (text = replace(text, "0"^64 => h; count = 1))
    return text, h
end

# ---------------------------------------------------------------- the designed orbits (verdict histories of the seven indicators at the ten budgets)
# Sets of Cal by design, after the relabelling C0 (rho_new): see the test plan, "Synthetic study", for the hand table.
const P3_DESIGN_NOTE = "S2 = 11 12 33..38; T_s = 5 6 23 24 25 39; W2 minus T_s = 1 3 7 8 21 22 26; other chaotic = 2 4 40; regular = 9 10 15 16 27..31 41 42 43; unresolved = 13 14 32 44"
function p3_design()
    prof = Dict{Int, Symbol}()
    for i in (11, 12, 33, 34, 35, 36, 37, 38); prof[i] = :S; end
    for i in (5, 6, 23, 24, 25, 39); prof[i] = :T; end
    for i in (1, 3, 7, 8, 21, 22, 26); prof[i] = :W; end
    for i in (2, 4, 40, 13, 14, 32, 44, 17, 18, 19, 20); prof[i] = :N; end
    for i in (9, 10, 15, 16, 27:31..., 41, 42, 43); prof[i] = :R; end
    # overrides: false-chaotic and unstable votes of regular orbits, split votes of weak orbits, undecided votes (the rule has something to decide about)
    ov = Dict{Tuple{Int, Symbol}, String}(
        (28, :MEGNO) => "rrrrrrcrrr",       # chaotic at b_P only (not stable at w >= 1), counts in n_Cany
        (29, :FLI) => "rrrrrrcrrr", (29, :E2) => "rrrrrrcrrr", (42, :SALI) => "rrrrrrcrrr", (30, :SALI) => "rrrrrrrrcc", (31, :E2) => "rrrrcrrrrr", (41, :WBA) => "uuuuuuuuuu", (42, :FLI) => "rrrrrrurrr", (43, :E1) => "rrrrrrcccc",
        (9, :E1) => "cccrrrrrrr", (10, :SALI) => "uurrrrrrrr", (15, :E1) => "uuuuuuuuuu", (16, :LE) => "rrrrrrrrrc",
        (3, :FLI) => "rrrrrrcrrr", (7, :WBA) => "uuuuuuuuuu", (21, :SALI) => "rrrrrrrccc", (26, :MEGNO) => "rrrrrrrrcc",
        # the tangent-map indicators of three strong orbits are undecided at every budget: every V of tg alone abstains on S2 (F3)
        (33, :SALI) => "uuuuuuuuuu", (33, :FLI) => "uuuuuuuuuu", (33, :MEGNO) => "uuuuuuuuuu", (34, :SALI) => "uuuuuuuuuu", (34, :FLI) => "uuuuuuuuuu", (34, :MEGNO) => "uuuuuuuuuu",
        (35, :SALI) => "uuuuuuuuuu", (35, :FLI) => "uuuuuuuuuu", (35, :MEGNO) => "uuuuuuuuuu",
        # a gate orbit that two indicators call chaotic at b_P only: a rule with w = 0 that trusts both is excluded by F4
        (10, :FLI) => "rrrrrrcrrr", (10, :E2) => "rrrrrrcrrr", (9, :WBA) => "rrrrrrcrrr",
        # a sticky orbit that WBA calls chaotic at b_P, and a regular orbit that SALI and WBA call chaotic: error to trade against coverage
        (24, :WBA) => "rrrrrrcccc", (27, :SALI) => "rrrrrrcrrr", (27, :WBA) => "rrrrrrcccc",
        # the sticky orbit 6 is regular to every voter at every budget (a chaotic orbit that looks regular): every rule that can call it regular makes a false-regular call (J > 0)
        (6, :SALI) => "rrrrrrrrrr", (6, :FLI) => "rrrrrrrrrr", (6, :MEGNO) => "rrrrrrrrrr", (6, :E1) => "rrrrrrrrrr", (6, :E2) => "rrrrrrrrrr", (6, :WBA) => "rrrrrrrrrr")
    nan_nu = Dict{Int, Vector{Int}}()
    rho = Dict{Int, Symbol}()
    for i in (1:8..., 11, 12, 21:26..., 33:40...); rho[i] = :C; end
    for i in (9, 10, 15, 16, 27:31..., 41, 42, 43); rho[i] = :R; end
    for i in (13, 14, 32, 44, 17, 18, 19, 20); rho[i] = :U; end     # the design class of the p2 grid is what the reduced rule gives; panel entries are used only by the p2 writer
    lref = Dict{Int, Float64}(i => (rho[i] === :R ? 1.4e-4 : 0.5) for i in 1:44)
    lns = copy(lref)
    for (j, i) in enumerate((33, 34, 35, 36, 37, 38, 11, 12)); lref[i] = lns[i] = 1.1 + 0.02 * (j - 1); end
    for i in (5, 6, 23, 24, 39); lns[i] = 0.05; end
    lns[25] = 0.09
    for (i, a) in ((1, 0.05), (3, 0.08), (7, 0.07), (8, 0.06), (21, 0.09), (22, 0.02), (26, 0.04)); lref[i] = lns[i] = a; end
    seqs = Dict{Int, Dict{Symbol, Vector{Symbol}}}()
    for id in 1:44
        pr = prof[id]
        d = Dict{Symbol, Vector{Symbol}}()
        for X in P2_X8
            s = pr === :R ? "rrrrrrrrrr" : pr === :T ? P2_TSTR[X] : (o = getfield(P2_ONSET[pr === :S ? :S : pr], X); repeat("r", o) * repeat("c", 10 - o))
            d[X] = p2_sym(get(ov, (id, X), s))
        end
        seqs[id] = d
    end
    haveRLI(id) = id in P2_SL || id >= 21
    for id in 1:44
        haveRLI(id) || delete!(seqs[id], :RLI)
        id == 20 && (delete!(seqs[id], :E1); delete!(seqs[id], :WBA))
    end
    d64 = Dict{Tuple{Int, Symbol}, Symbol}((5, :E2) => :regular, (11, :LE) => :regular)
    return (variant = :p3, prof = prof, rho = rho, lref = lref, lns = lns, V = seqs, d64 = d64, dir = Dict{Tuple{Int, Symbol, Int}, Symbol}(), drift = Int[], nan_nu = nan_nu)
end

# the expected rho_new of the fixture by the reduced rule C0 (E49): panel orbits from the P1 items, grid orbits from the p2 items
const P3_RHO_NEW = merge(Dict{Int, Symbol}(i => :C for i in (1:8..., 11, 12, 21:26..., 33:40...)), Dict{Int, Symbol}(i => :R for i in (9, 10, 15, 16, 27:31..., 41, 42, 43)),
                         Dict{Int, Symbol}(i => :U for i in (13, 14, 32, 44)))

# ---------------------------------------------------------------- the P1 reference study with N_1 equal to the p2 horizon
# orbit 13: the L5 item read for C0 carries "diagnostics exceeded" (the diagnostics condition of E49 c); orbit 14: the P1 L1 item differs from its L5 item by 1e-6 (the
# L1 against L5 test of E49 c, with the P1 L1 item as in P45); orbit 2 carries the P1 ladder defect of the p2 fixture (an old label U that the reduced rule does not give).
function p3_p1_study(dir; p = an_fx_params(N_1 = P2_N1))
    prim(sp) = sp.method == "L1" && "main" in sp.arms && sp.mode == "lifted" && sp.precision == "Float64" && sp.horizon == p.N_1 && sp.tau == 1 && sp.direction == "pi/4"
    return p2_p1_study(dir; p, cover = Dict(2 => Dict("L9" => 100.0)),
        flags = (id, sp) -> (id == 13 && sp.method == "L5" && "main" in sp.arms) ? ["diagnostics exceeded"] : String[],
        lam = (id, sp, steps, vals) -> (id == 14 && prim(sp)) ? vals .+ 1e-6 : vals,
        attrs = (id, sp) -> id <= 20 ? Dict{String, Any}("ic" => "$(AN_ORBITS[id]["x0_dec"]),$(AN_ORBITS[id]["p0_dec"])") : Dict{String, Any}())     # the stored start text equals the panel file (E48 g)
end
const P3_P1 = Ref{Any}(nothing)
p3_p1_cached() = (P3_P1[] === nothing && (P3_P1[] = p3_p1_study(an_path("p3-p1-shared"))); P3_P1[])

# the orbit rows of the calibration straight from a design (the verdicts as designed, rho_new as designed): what the p2 store holds when it realises the design (checked by
# the fixture self-test); for tests that need no store
function p3_design_rows(design; q = p3_params(), rho = P3_RHO_NEW)
    h4 = p2_h4_tables(p2_params())
    rows = NamedTuple[]
    for id in vcat(1:16, 21:44)
        r = rho[id]
        C = r === :C
        lref, lns = design.lref[id], design.lns[id]
        Ts = C && lns <= q.f_s * lref
        push!(rows, (id = id, source = id <= 20 ? :panel : :grid, kindex = id <= 20 ? Int(AN_ORBITS[id]["k_index"]) : h4[id - 20].kindex, in_cal = r in (:C, :R), rho = r,
                     S2 = C && lref >= q.lam_s2, Ts = Ts, W2 = C && (lref <= q.lam_w2 || Ts), gate = id in (9, 10), V = Dict(X => design.V[id][X] for X in P2_X7)))
    end
    return rows
end

# ---------------------------------------------------------------- the p2 store of the study (the body of p2_study with a given design), and the p3 inputs
function p3_make_p2(tag::AbstractString, design, p1; p = nothing)
    dir = an_path("p3-$tag")
    rm(dir; recursive = true, force = true)
    mkpath(dir)
    p1dir = joinpath(dir, "p1")
    cp(p1.dir, p1dir)
    sub(s) = replace(s, p1.dir => p1dir)
    p1c = merge(p1, (; dir = p1dir, store = sub(p1.store), prereg = sub(p1.prereg), panel = sub(p1.panel), manifest = sub(p1.manifest), params = sub(p1.params)))
    panel = joinpath(dir, "$P2_SLUG-panel.toml")
    cp(p1c.panel, panel)
    ph = an_sha256(read(panel))
    p = p === nothing ? p2_params(; panel_sha256 = ph) : merge(p, (; panel_sha256 = ph))
    text, h = p2_prereg_text(p)
    prereg = joinpath(dir, "$P2_SLUG-preregistration.md"); write(prereg, text)
    params = p2_write_params(joinpath(dir, "$P2_SLUG-analysis-parameters.toml"), p)
    h4 = p2_h4_tables(p)
    store = joinpath(dir, "store"); mkpath(store)
    specs = p2_specs(p)
    chunks = Dict{String, Any}[]
    total = 0
    for (ord, sp) in enumerate(specs)
        recs = AX.OrbitRecord[]
        for id in sp.orbits
            info = p2_info(id, h4)
            ex = Dict{String, Any}()
            if "lyap_p" in sp.arms && sp.role === :L1 && id in (17, 18, 20)
                ex["control_scope"] = id == 18 ? "info" : "bound"
                ex["control_deviation"] = abs(design.lref[id] - Float64(AN_ORBITS[id]["exponent"]))
            end
            push!(recs, p2_orbit_record(sp, info, design, p; flags = id in design.drift ? ["drift"] : String[], extra = ex))
        end
        total += length(recs)
        key = "method=$(sp.id) kind=$(sp.kind) precision=$(sp.precision) horizon=$(sp.horizon) tau=1 direction=$(sp.direction) mode=wrapped arm=$(sp.arms[1])"
        ch = "c"^60 * lpad(string(ord), 4, '0')
        root = Dict{String, Any}("schema_version" => 2, "study" => P2_SLUG, "chunk" => ord, "chunk_key" => key, "chunk_hash" => ch, "plan_sha256" => p.plan_sha256,
            "panel_sha256" => ph, "prereg_sha256" => h, "git_sha" => "a"^40, "code_version" => "an-code-1", "backend_request" => "cpu", "backend" => "cpu",
            "device" => "cpu", "backend_fallback_reason" => "", "precision" => sp.precision, "method" => sp.id, "horizon" => sp.horizon, "tau" => 1,
            "vectors" => sp.role in (:tg, :L5) ? 2 : 1, "coordinate_mode" => "wrapped", "grid_sha256" => "g"^64,
            "status" => "ok", "attempt" => 1, "setup_check" => "not applicable", "first_attempt_reasons" => String[], "failure_reasons" => String[],
            "pilot" => false, "arms" => sp.arms, "session" => 1, "control_stop" => false, "dirty" => false, "created_utc" => "2026-10-29T00:00:00Z")
        sp.role in (:L1, :L5, :tg, :st) && (root["direction"] = sp.direction)
        if sp.kind != ""
            kd, mem = P2_SPEC_INFO[sp.id]
            st_text = "indicator.id = \"$(sp.id)\"\nindicator.kind = \"$kd\"\n"
            merge!(root, Dict{String, Any}("kind" => kd, "spec" => st_text, "spec_sha256" => an_sha256(st_text), "chain" => [kd], "members" => mem))
            kd in ("frequency", "phase_advance", "birkhoff") && merge!(root, Dict{String, Any}("windows" => p.windows, "run_length" => sp.horizon))
            kd == "birkhoff" && (root["observable_names"] = P2_BKROWS)
        end
        AX.write_chunk(joinpath(store, "chunk-" * lpad(string(ord), 5, '0') * ".h5"), root, recs)
        ent = Dict{String, Any}("ordinal" => ord, "key" => key, "hash" => ch, "arms" => sp.arms, "status" => "ok", "attempt" => 1, "orbits" => length(recs), "backend" => "cpu",
            "steps" => sp.horizon * length(recs), "wall_seconds" => 1.0, "thread_seconds" => 2.0, "fallback_reason" => "", "first_attempt_reasons" => String[],
            "setup_check" => "not applicable", "control_stop" => false, "control_info" => Any[],
            "flags" => Dict{String, Any}("diagnostics exceeded" => Int[], "drift" => design.drift, "excluded" => Int[], "breakdown_tau" => Int[]))
        if sp.kind != ""
            ent["indicator_flags"] = Dict{String, Any}(k => v for (k, v) in p2_flag_counts(recs))
            sp.kind == "frequency" && (ent["transform_seconds"] = 1.5)
            sp.role in (:fq1, :bk1) && sp.arms[1] == "win_p" && (ent["not_applicable"] = Dict{String, Any}(sp.id => [20]))
        end
        push!(chunks, ent)
    end
    arms = Dict{String, Any}()
    for sp in specs, a in sp.arms
        e = get!(arms, a, Dict{String, Any}("status" => "complete", "steps" => 0, "cost_ns_per_orbit_step" => 100.0, "budget_ns" => 415.2, "specs" => Dict{String, Any}()))
        haskey(P2_COST, sp.id) && (e["specs"][sp.id] = Dict{String, Any}("budget_line_ns" => P2_COST[sp.id][1], "measured_ns_per_orbit_step" => P2_COST[sp.id][2]))
    end
    inds = Dict{String, Any}(id => Dict{String, Any}("kind" => P2_SPEC_INFO[id][1], "windows" => p.windows,
                                                    "bytes_per_orbit" => Dict("B" => 0, "M" => 80, "R" => 24, "S" => 9128, "total" => 9232 + 100 * i),
                                                    "budget" => Dict("Float64" => Dict("ns_cpu_thread" => P2_COST[id][1] / 2, "safety" => 2.0)))
                             for (i, id) in enumerate(sort(collect(keys(P2_SPEC_INFO)))))
    manifest = joinpath(dir, "$P2_SLUG-manifest.json")
    AX._write_json(manifest, Dict{String, Any}("study" => P2_SLUG, "schema_version" => 2, "status" => "complete", "plan_sha256" => p.plan_sha256,
        "plan_counts" => Dict("items" => total, "chunks" => length(chunks), "steps_before" => 1, "steps_after" => 1),
        "git" => Dict("sha" => "b"^40, "clean" => true), "sessions" => Any[Dict("clean" => true, "git_sha" => "b"^40, "index" => 1, "threads_used" => 4, "wall_seconds" => 10.0)],
        "prereg" => Dict("header" => h, "sha256_start" => h, "sha256_end" => h),
        "panel" => Dict("sha256_start" => ph, "sha256_end" => ph, "build_hash" => p.panel_build_sha256), "overrides" => Dict{String, Any}(), "control_stops" => Any[],
        "reasons" => String[], "chunks" => chunks, "arms" => arms, "indicators" => inds))
    return (; dir, store, prereg, panel, manifest, params, p, design, specs, chunks = length(chunks), items = total, study = P2_SLUG, prereg_hash = h, panel_hash = ph,
            out = joinpath(dir, "out"), h4, p1 = p1c, variant = :p3)
end

# the whole study: the p2 study (a p2-like NamedTuple for the p2 reference helpers), the p3 preregistration and the p3 parameter table
function p3_study(tag::AbstractString; design = p3_design(), q = p3_params(), p1 = p3_p1_cached(), p = nothing)
    s2 = p3_make_p2(tag, design, p1; p)
    q = merge(q, (; p2_manifest_sha256 = an_sha256(read(s2.manifest))))
    text, h = p3_prereg_text(q)
    prereg = joinpath(s2.dir, "$P3_SLUG-preregistration.md"); write(prereg, text)
    params = p3_write_params(joinpath(s2.dir, "$P3_SLUG-calibration-parameters.toml"), q)
    return (; s2, q, dir = s2.dir, prereg, params, prereg_hash = h, out = joinpath(s2.dir, "out"), store = s2.store, p1 = s2.p1, design)
end
const P3_STUDIES = Dict{String, Any}()
p3_cached(tag; kw...) = get!(() -> p3_study(tag; kw...), P3_STUDIES, tag)

# copy of a study into a new directory, with the paths rewritten
function p3_copy(st, tag::AbstractString)
    dir = an_path("p3copy-$tag")
    rm(dir; recursive = true, force = true)
    cp(st.dir, dir)
    rm(joinpath(dir, "out"); recursive = true, force = true)
    sub(s) = replace(s, st.dir => dir)
    p1 = merge(st.p1, (; dir = sub(st.p1.dir), store = sub(st.p1.store), prereg = sub(st.p1.prereg), panel = sub(st.p1.panel), manifest = sub(st.p1.manifest), params = sub(st.p1.params)))
    s2 = merge(st.s2, (; dir, store = sub(st.s2.store), prereg = sub(st.s2.prereg), panel = sub(st.s2.panel), manifest = sub(st.s2.manifest), params = sub(st.s2.params),
                       out = joinpath(dir, "out"), p1))
    return merge(st, (; s2, dir, prereg = sub(st.prereg), params = sub(st.params), out = joinpath(dir, "out"), store = s2.store, p1))
end

# the call of the command in the form of keywords and in the form of arguments
p3_run(st; params = st.params, prereg = st.prereg, out = st.out, kw...) = AN.calibrate(; prereg, params, out, p2_store = st.s2.store, p2_manifest = st.s2.manifest,
    p2_prereg = st.s2.prereg, p2_params = st.s2.params, panel = st.s2.panel, p1_store = st.p1.store, p1_manifest = st.p1.manifest, p1_prereg = st.p1.prereg,
    p1_params = st.p1.params, kw...)
p3_cli_args(st; extra = String[]) = vcat(["calibrate", st.prereg, st.params, "--out", st.out, "--p2-store", st.s2.store, "--p2-manifest", st.s2.manifest, "--p2-preregistration", st.s2.prereg,
    "--p2-parameters", st.s2.params, "--panel", st.s2.panel, "--p1-store", st.p1.store, "--p1-manifest", st.p1.manifest, "--p1-preregistration", st.p1.prereg,
    "--p1-parameters", st.p1.params], extra)
p3_chunk(st, i) = p2_chunk(st.s2, i)
p3_edit_manifest(st, f) = p2_edit_manifest(st.s2, f)

# a hash of every byte of every input and store file (the read-only checks, I30); both stores, all manifests, tables and preregistrations
p3_fingerprint(st) = [hash(String(read(f))) for f in vcat([st.prereg, st.params, st.s2.manifest, st.s2.prereg, st.s2.panel, st.s2.params, st.p1.manifest, st.p1.prereg, st.p1.params],
    sort(filter(f -> endswith(f, ".h5"), readdir(st.s2.store; join = true))), sort(filter(f -> endswith(f, ".h5"), readdir(st.p1.store; join = true))))]

# ---------------------------------------------------------------- pure inputs: orbit rows for the selection (E50, E54) built by hand
# rows: (id, source, kindex, in_cal, rho, S2, Ts, W2, gate, V); `vs` maps an indicator to the verdict string of the ten budgets ('c', 'r', 'u'); an indicator that is not
# listed gets the string `fill`
function p3_row(id; source = :panel, kindex = 2, in_cal = true, rho = :C, S2 = false, Ts = false, W2 = Ts, gate = false, fill = "rrrrrrrrrr", vs = Dict{Symbol, String}())
    V = Dict{Symbol, Vector{Symbol}}(Symbol(X) => p2_sym(get(vs, Symbol(X), fill)) for X in P3_I7)
    return (id = id, source = source, kindex = kindex, in_cal = in_cal, rho = rho, S2 = S2, Ts = Ts, W2 = W2, gate = gate, V = V)
end
# a verdict history with the same verdict `c` at the budgets 4 to 9 (b_P is index 6): all voters agree
p3_hist(c::Char) = repeat(string(c), 10)

# the model of one simulated draw (E56 c): rule and comparator rates per stratum (:T and :RW); a comparator that is not listed has rate 0
function p3_model(; rule_T = (p = 0.1, d = 1.0), rule_RW = (p = 0.1, d = 1.0), comp = Dict{Symbol, Any}(), ts = Symbol[], w2 = Symbol[])
    z = (p = 0.0, s = 0.0)
    cm = Dict{Symbol, Any}(X => get(comp, X, (T = z, RW = z)) for X in P2_X7)
    return (rule = Dict(:T => rule_T, :RW => rule_RW), comp = cm, ts = ts, w2 = w2)
end

# ---------------------------------------------------------------- prescribed tables of the selection tests (O18), shared by an_p3_select.jl and an_p3_pp.jl
# five voters, |V| >= 4, one depth and one v: 55 elements, of which the 8 of V = SALI, FLI, MEGNO, E2 have no E1 (the cost cap F1 removes the 47 others)
p3_toy_q(; kw...) = merge(p3_params(), (; voter_bits = ["SALI", "FLI", "MEGNO", "E1", "E2"], indicators_i7 = ["LE", "SALI", "FLI", "MEGNO", "E1", "E2"], min_voters = 4, w_range = [0], v_range = [0],
    min_set = 3, e_cal = 0.4, s_cal = 0.15, a_cal = 0.3, c_min = 0.5, delta_tol = 0.05), NamedTuple(kw))
# an orbit votes with its first nc voters chaotic and the next nr regular at every budget; the others (and E1, LE, WBA) are undecided
p3_votes(nc, nr; names = [:SALI, :FLI, :MEGNO, :E2]) = Dict{Symbol, String}(X => (i <= nc ? "c"^10 : i <= nc + nr ? "r"^10 : "u"^10) for (i, X) in enumerate(names))
function p3_toy_rows(spec)         # spec: vector of (nc, nr, kw...)
    rows = NamedTuple[]
    for (i, (nc, nr, kw)) in enumerate(spec)
        push!(rows, p3_row(i; fill = "u"^10, vs = p3_votes(nc, nr), kw...))
    end
    return rows
end

# the prescribed table of O18 (a): the element E* = (V = SALI, FLI, MEGNO, E2; q_C = 2, q_R = 3) is the only feasible one, by construction:
#  - S2 has four orbits, one of them with two chaotic votes only: q_C = 3 or 4 leaves it unresolved (u(S2) = 1/4 > s_cal): F3
#  - two regular orbits and the gate orbit have one stray chaotic vote: q_C = 1 calls them chaotic: error 7/17 > e_cal (F2) and the gate fires (F4)
#  - the W2 orbits outside T_s have three regular votes only: q_R = 4 leaves them unresolved, the decided share on W2 is 2/7 < 1/2 (F6); the share of unresolved on Cal exceeds a_cal (F5)
#  - every V that contains E1 costs more than the cap (F1)
const P3_TOY_A = [(4, 0, (S2 = true,)), (4, 0, (S2 = true,)), (4, 0, (S2 = true,)), (2, 0, (S2 = true,)),
                  (3, 0, (Ts = true,)), (3, 0, (Ts = true,)), (0, 3, (Ts = true,)),
                  (0, 3, (W2 = true,)), (0, 3, (W2 = true,)), (0, 3, (W2 = true,)), (0, 3, (W2 = true,)),
                  (0, 4, (rho = :R,)), (0, 4, (rho = :R,)), (0, 4, (rho = :R,)), (0, 4, (rho = :R,)), (1, 0, (rho = :R,)), (1, 0, (rho = :R,)), (1, 0, (rho = :R, gate = true))]

same_sel(a, b) = a.outcome === b.outcome && a.star == b.star && a.feasible == b.feasible && a.tol == b.tol && a.first_fail == b.first_fail && a.alone == b.alone && a.J == b.J &&
                 a.J_min == b.J_min && all(i -> a.flags[i] == b.flags[i] && a.stats[i] == b.stats[i], eachindex(b.thetas)) && length(a.thetas) == length(b.thetas)


# the same design with three sticky orbits made ordinary (lambda_ns = lambda_ref): T_s has 3 < min_set = 5 orbits in Cal, the outcome is "sets too small"
function p3_design_small()
    d = p3_design()
    for i in (24, 25, 39); d.lns[i] = d.lref[i]; end
    return d
end

# the table of a parameter set that makes F empty (no rule is feasible): every threshold is as strict as the designed study can meet only partly
p3_q_empty(q) = merge(q, (; e_cal = 0.0, a_cal = 0.0, s_cal = 0.0, c_min = 1.0))

# the orbit rows of the calibration for the reference implementation: the designed verdicts read back from the chunk files by the independent p2 reference (never the module)
function p3_ref_rows(st; rho = P3_RHO_NEW)
    ref = p2_ref_all(st.s2)
    d, q = st.design, st.q
    rows = NamedTuple[]
    for id in vcat(1:16, 21:44)
        r = rho[id]
        C = r === :C
        lref, lns = d.lref[id], d.lns[id]
        Ts = C && lns <= q.f_s * lref
        push!(rows, (id = id, source = id <= 20 ? :panel : :grid, kindex = id <= 20 ? Int(AN_ORBITS[id]["k_index"]) : st.s2.h4[id - 20].kindex, in_cal = r in (:C, :R), rho = r,
                     S2 = C && lref >= q.lam_s2, Ts = Ts, W2 = C && (lref <= q.lam_w2 || Ts), gate = id in (9, 10), V = Dict(X => ref.V1[(id, X)] for X in P2_X7)))
    end
    return rows
end


# ---------------------------------------------------------------- editing helpers for the malformed inputs and the amendment of the preregistration (used by an_p3_refusals.jl and an_p3_pipeline.jl)
# the preregistration with its hash cell rewritten by the SELF rule, so that the edit is the only defect
function p3_edit_prereg(c, f)
    text = f(read(c.prereg, String))
    prep(t) = replace(replace(t, "\r\n" => "\n"), r"(?m)^(\| SHA-256 of this file \(current version\) \|)[^\n|]*\|" => s"\1 SELF |")
    h = an_sha256(prep(text))
    write(c.prereg, replace(text, r"(?m)^(\| SHA-256 of this file \(current version\) \| )[0-9a-f]{64}" => SubstitutionString("\\g<1>" * h)))
    return h
end
# the p3 preregistration and table follow the SHA-256 of the (edited) p2 manifest, as the real table would after a regeneration of the p2 manifest
function p3_resync!(c)
    q = merge(c.q, (; p2_manifest_sha256 = an_sha256(read(c.s2.manifest))))
    text, _ = p3_prereg_text(q); write(c.prereg, text)
    p3_write_params(c.params, q)
    return q
end
function p3_rewrite(c, i, f)                    # a chunk file of the p2 store rewritten with the real writer after a change of its orbit records
    path = p3_chunk(c, i)
    ch = AX.read_chunk(path)
    AX.write_chunk(path, Dict{String, Any}(ch.attrs), AX.OrbitRecord[f(o) for o in ch.orbits if f(o) !== nothing])
end
p3_p1chunk(c, i) = joinpath(c.p1.store, "chunk-" * lpad(string(i), 5, '0') * ".h5")


# the keys of the outermost object of a JSON text in file order (a small scanner: strings with escapes, nesting depth), so that the test does not depend on the layout
function p3_top_keys(text::AbstractString)
    keys_ = String[]
    depth = 0; i = firstindex(text); n = lastindex(text)
    while i <= n
        c = text[i]
        if c == '"'
            j = nextind(text, i); buf = IOBuffer()
            while j <= n && text[j] != '"'
                text[j] == '\\' && (j = nextind(text, j))
                j <= n && print(buf, text[j])
                j = nextind(text, j)
            end
            k = nextind(text, j)
            while k <= n && isspace(text[k]); k = nextind(text, k); end
            depth == 1 && k <= n && text[k] == ':' && push!(keys_, String(take!(buf)))
            i = nextind(text, j)
            continue
        end
        (c == '{' || c == '[') && (depth += 1)
        (c == '}' || c == ']') && (depth -= 1)
        i = nextind(text, i)
    end
    return keys_
end

# a parameter set written as the preregistration and the table of the study (the quotes follow the values)
function p3_set_q!(c, q)
    text, _ = p3_prereg_text(q)
    write(c.prereg, text)
    p3_write_params(c.params, q)
    return q
end
