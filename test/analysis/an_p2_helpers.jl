# Amendment D89, stage 02: helpers of the P2 suite (prefix p2_). Definitions only: nothing here calls Chaoslab.Analysis at load time, so a missing
# function errors the testsets that need it and nothing else. Uses an_helpers.jl (AN, AX, an_sha256, an_path, an_study pieces) without changing it.
# Spec: workspaces/module-build/stages/01-spec/output/analysis-spec.md, "Amendment D89". API assumed: analysis-test-plan.md, "API proposed by the tests".
#
# The synthetic study "p2-indicators" (scaled): 16 panel orbits (ids 1..16) and 4 controls (17..20, Y5 = 20) of the P1 fixture panel, a grid G of 24 starts
# (ids 21..44: 2 K values x 3 x 4 lattice), 33 chunk files written by the real chunk writer, a manifest, a parameter table of 58 keys, a fixture
# preregistration with run blocks, and a P1 reference study (slug p1-le-benchmark) for the labels. Budgets 40 x 2^k, windows 20 x 2^k, N_c = 20480,
# N_1 = 80000, N_dir = 40000, n_s = 800, min_set = 5. All thresholds are the frozen study values (rules R1, R2), except label_lambda_min (the P1 fixture value).

p2_seed(letter::Char, sub::Integer = 0) = AN_SEED + 31000 + 50 * (letter - 'a' + 1) + sub      # D89<letter> draws use 31000 + 50 n + sub, sub < 50
const P2_SLUG, P1_SLUG = "p2-indicators", "p1-le-benchmark"
const P2_BUDGETS = [40 * 2^k for k in 0:9]
const P2_WINDOWS = [20 * 2^k for k in 0:9]
const P2_NC, P2_N1, P2_NDIR, P2_NS = 20_480, 80_000, 40_000, 800
const P2_DIRS = ["pi/4", "0.3", "0.3 + pi/5", "0.3 + 2pi/5", "0.3 + 3pi/5"]
const P2_KTEXT = ["0.5", "6.283185307179586"]
const P2_KP = 7                                   # 1-based position of b_primary = 2560 in the budget list (budget index 6)
const P2_SD = [1, 2, 5, 7, 9, 11, 12, 15, 17, 18, 19, 20]
const P2_SL = [1, 5, 9, 11]

# recorded grids (E26). GA: LE and tangent items of horizon N_1; GAD: direction and Double64 items (horizon N_dir); GB: pair, stretching and window items
# (horizon N_c); GB1: phase advance (N_c + 1). n(b_k) on GA = 44 x 2^k, n' = 4, 9, 18, 35, 70, 141, 282, 563, 1126, 2253 (hand table in the plan).
const P2_NK = [44 * 2^k for k in 0:9]
const P2_NP = [4, 9, 18, 35, 70, 141, 282, 563, 1126, 2253]
const P2_GA = sort(unique(vcat(1:10, P2_NP, P2_NK, [P2_NS, 8_000, 40_000, 80_000])))
const P2_GAD = filter(<=(P2_NDIR), P2_GA)
const P2_GB = sort(unique(vcat(filter(<=(5632), P2_GA), [2048, 10_240, 20_480])))
const P2_GB1 = vcat(P2_GB, 20_481)
const P2_NKB = vcat(P2_NK[1:8], [10_240, 20_480])  # n(b_k) on GB (pair and stretching)

# ---------------------------------------------------------------- parameter table (58 keys, E22)
const P2_PARAM_KEYS = [:budgets, :b_primary, :N_c, :N_1, :N_dir, :n_s, :f_s, :windows, :directions,
    :ratio_threshold, :le_floor, :rli_threshold, :e1_d_threshold, :e2_d_threshold, :wba_digits_threshold,
    :r2_le_factor, :r2_sali_threshold, :r2_fli_factor, :r2_megno_threshold, :persist_divisor, :censor_factor,
    :s2_lambda_min, :w2_lambda_max, :gate_k_index, :gate_classes, :gate_size, :drift_bound,
    :label_lambda_min, :label_growth_ratio, :label_rotation_tol, :label_estimator_bound,
    :h1_u_min, :h1_u_lower_min, :h1_amin_max, :h2_q_min, :h2_q_lower_gt, :h3_r_min, :h3_r_lower_gt, :h3_count_min,
    :h4_dphi_min, :h4_dphi_lower_min, :h5_fmed_min, :h5_fmed_lower_min, :h6_flip_min,
    :bootstrap_resamples, :bootstrap_level, :bootstrap_seed, :min_set, :range_lower_percentile, :range_upper_percentile, :gate_max_share,
    :grid_k_values, :grid_i_count, :grid_j_count, :grid_orbit_count, :panel_sha256, :panel_build_sha256, :plan_sha256]
const P2_MACHINE = [:N_1, :N_c, :N_dir, :windows, :directions, :drift_bound, :grid_k_values, :grid_i_count, :grid_j_count, :grid_orbit_count, :panel_build_sha256]
const P2_INT_KEYS = [:b_primary, :N_c, :N_1, :N_dir, :n_s, :persist_divisor, :censor_factor, :gate_k_index, :gate_size, :h3_count_min, :bootstrap_resamples,
                     :bootstrap_seed, :min_set, :grid_i_count, :grid_j_count, :grid_orbit_count]

p2_params(; kw...) = merge((; budgets = P2_BUDGETS, b_primary = 2560, N_c = P2_NC, N_1 = P2_N1, N_dir = P2_NDIR, n_s = P2_NS, f_s = 0.2, windows = P2_WINDOWS,
    directions = P2_DIRS, ratio_threshold = 5.0, le_floor = 1e-4, rli_threshold = 1e-7, e1_d_threshold = -4.0, e2_d_threshold = -4.0, wba_digits_threshold = 6.0,
    r2_le_factor = 3.0, r2_sali_threshold = 1e-8, r2_fli_factor = 3.0, r2_megno_threshold = 2.5, persist_divisor = 4, censor_factor = 2,
    s2_lambda_min = 1.0, w2_lambda_max = 0.1, gate_k_index = 1, gate_classes = ["T"], gate_size = 2, drift_bound = 1e6,
    label_lambda_min = 1e-2, label_growth_ratio = 5.0, label_rotation_tol = 1e-7, label_estimator_bound = 1e-8,
    h1_u_min = 0.2, h1_u_lower_min = 0.1, h1_amin_max = 0.8, h2_q_min = 10.0, h2_q_lower_gt = 3.0, h3_r_min = 10.0, h3_r_lower_gt = 3.0, h3_count_min = 5,
    h4_dphi_min = 0.5, h4_dphi_lower_min = 0.25, h5_fmed_min = 0.1, h5_fmed_lower_min = 0.05, h6_flip_min = 0.05,
    bootstrap_resamples = 1000, bootstrap_level = 0.95, bootstrap_seed = AN_SEED, min_set = 5, range_lower_percentile = 5.0, range_upper_percentile = 95.0,
    gate_max_share = 0.1, grid_k_values = P2_KTEXT, grid_i_count = 3, grid_j_count = 4, grid_orbit_count = 24,
    panel_sha256 = "0"^64, panel_build_sha256 = "e"^64, plan_sha256 = "p"^64), NamedTuple(kw))

# section headings of the fixture preregistration (the real ones) and, per key, the verbatim quote that contains the value
const P2_SEC = (M1 = "M1 Budgets and horizons", M2 = "M2 Verdicts: score and class rule of each indicator", M4 = "M4 Detection time and outcomes",
    M8 = "M8 Sets, bootstrap, undecidable sets", M10 = "M10 Regular range of sticky orbits (naff-fma R4, weighted-birkhoff R4)",
    H1 = "H1 (claim C3, classification divergence at a fixed budget)", H2 = "H2 (claim C3, detection-time divergence)",
    H3 = "H3 (detection latency, regime dependence)", H4 = "H4 (claim C3, sticky orbits)", H5 = "H5 (classification stability over the budget)",
    H6 = "H6 (classification stability under direction and precision)", HYP = "Hypotheses", PANEL = "Panel P (reused by reference)",
    RHO = "Reference class rho(o)", GRID = "Grid G (new starts, K x initial-condition)", SETS = "Weak-chaos set W2, strong-chaos set S2, sticky ensemble T_s (R4)",
    ELIG = "Eligible sets", RUN = "Run block (draft)", PLAN = "Analysis plan", OPEN = "Open points for the owner at the checkpoint",
    A1 = "Amendment A1 (2026-10-29)", A2 = "Amendment A2 (2026-10-29)")
const P2_SEC_ORDER = [:HYP, :H1, :H2, :H3, :H4, :H5, :H6, :PANEL, :RHO, :GRID, :SETS, :ELIG, :RUN, :M1, :M2, :M4, :M8, :M10, :PLAN, :OPEN, :A1, :A2]

# The text of a value as the section writes it (A2): the field `text` of the table, which must occur in the quote. It differs from Julia's string(value) for many keys
# (1e-4, 0.20, -4, 10^3, 95 percent, 5th), so that a reader that matches string(value) is caught.
function p2_texts(p; style = :section, texts = Dict{Symbol, String}())
    pow10(x) = (e = round(Int, log10(x)); x == 10^e ? "10^$e" : string(x))
    T = Dict{Symbol, String}()
    for k in P2_PARAM_KEYS
        v = getfield(p, k)
        T[k] = v isa AbstractVector ? join(v, k === :gate_classes ? " and " : ", ") : string(v)
    end
    style === :julia && return merge(T, texts)                      # the text is Julia's string(value) everywhere (a variant the reader must accept as well)
    sci = Dict(:le_floor => "1e-4", :rli_threshold => "1e-7", :r2_sali_threshold => "1e-8", :label_rotation_tol => "1e-7", :label_estimator_bound => "1e-8")
    for (k, t) in sci
        parse(Float64, t) == getfield(p, k) && (T[k] = t)
    end
    two(x) = (u = string(round(x; digits = 2)); length(split(u, '.')[2]) < 2 ? u * "0" : u)
    for k in (:h1_u_min, :h1_u_lower_min, :h1_amin_max, :h4_dphi_min, :h5_fmed_min, :h5_fmed_lower_min, :h6_flip_min, :gate_max_share, :h4_dphi_lower_min)
        T[k] = two(getfield(p, k))
    end
    for k in (:e1_d_threshold, :e2_d_threshold, :wba_digits_threshold, :ratio_threshold, :label_growth_ratio, :r2_le_factor, :r2_fli_factor, :h2_q_min, :h2_q_lower_gt, :h3_r_min, :h3_r_lower_gt)
        v = getfield(p, k)
        isinteger(v) && (T[k] = string(Int(v)))
    end
    T[:bootstrap_resamples] = pow10(p.bootstrap_resamples)
    T[:bootstrap_level] = p.bootstrap_level == 0.95 ? "95 percent" : string(p.bootstrap_level)
    T[:range_lower_percentile] = isinteger(p.range_lower_percentile) ? "$(Int(p.range_lower_percentile))th" : string(p.range_lower_percentile)
    T[:range_upper_percentile] = isinteger(p.range_upper_percentile) ? "$(Int(p.range_upper_percentile))th" : string(p.range_upper_percentile)
    return merge(T, texts)                            # `texts` overrides single texts (the quote follows the text)
end

function p2_quotes(p; style = :section, texts = Dict{Symbol, String}())
    s = P2_SEC
    T = p2_texts(p; style, texts)
    q = Dict{Symbol, Tuple{String, Any, String}}()
    for k in P2_MACHINE
        q[k] = (k in (:N_1, :N_c, :N_dir, :windows, :directions) ? s.RUN : k === :drift_bound ? s.RUN : k === :panel_build_sha256 ? s.PANEL : s.A1, :machine, "")
    end
    put(k, sec, text) = (q[k] = (sec, text, T[k]))
    put(:budgets, s.M1, "Budget list B = {b_k = 40 x 2^k, k = 0 to 9} = " * T[:budgets])
    put(:b_primary, s.M1, "Primary budget b_P = b_6 = $(T[:b_primary])")
    put(:n_s, s.SETS, "T_s (sticky): lambda_1(n_s) <= f_s lambda_ref, with n_s = $(T[:n_s])")
    put(:f_s, s.SETS, "sticky if lambda_1(n_s) <= f_s lambda_ref with f_s = $(T[:f_s])")
    put(:s2_lambda_min, s.SETS, "S2 (strong chaos): lambda_ref >= $(T[:s2_lambda_min])")
    put(:w2_lambda_max, s.SETS, "W2 (weak chaos): lambda_ref <= $(T[:w2_lambda_max]), or the orbit is in T_s")
    put(:gate_k_index, s.ELIG, "K = 0 gate set: the $(T[:gate_size]) orbits of the panel cells (K index $(T[:gate_k_index]), classes $(T[:gate_classes]))")
    put(:gate_classes, s.ELIG, "K = 0 gate set: the $(T[:gate_size]) orbits of the panel cells (K index $(T[:gate_k_index]), classes $(T[:gate_classes]))")
    put(:gate_size, s.ELIG, "K = 0 gate set: the $(T[:gate_size]) orbits of the panel cells (K index $(T[:gate_k_index]), classes $(T[:gate_classes]))")
    put(:ratio_threshold, s.M2, "SALI | r_SALI | r_SALI >= $(T[:ratio_threshold])")
    put(:le_floor, s.M2, "LE | lambda_1(n), r_LE | lambda_1(n) >= $(T[:le_floor]) and r_LE >= $(T[:ratio_threshold])")
    put(:rli_threshold, s.M2, "RLI | column RLI | RLI >= $(T[:rli_threshold])")
    put(:e1_d_threshold, s.M2, "E1 | D_E1 | D_E1 >= $(T[:e1_d_threshold])")
    put(:e2_d_threshold, s.M2, "E2 | D_E2 | D_E2 >= $(T[:e2_d_threshold])")
    put(:wba_digits_threshold, s.M2, "WBA | dig_min | dig_min < $(T[:wba_digits_threshold])")
    put(:r2_le_factor, s.M2, "LE: n lambda_1(n) >= $(T[:r2_le_factor]) ln n")
    put(:r2_sali_threshold, s.M2, "SALI: SALI(n) <= $(T[:r2_sali_threshold])")
    put(:r2_fli_factor, s.M2, "FLI: fli(n) >= $(T[:r2_fli_factor]) ln n")
    put(:r2_megno_threshold, s.M2, "MEGNO: Ybar20_1(n) >= $(T[:r2_megno_threshold])")
    put(:persist_divisor, s.M4, "and b_k <= N_c / $(T[:persist_divisor])")
    put(:censor_factor, s.M4, "t_det for correct, and $(T[:censor_factor]) N_c for wrong and undetected")
    put(:label_lambda_min, s.RHO, "instability axis: lambda_L5 >= $(T[:label_lambda_min])")
    put(:label_growth_ratio, s.RHO, "growth ratio r >= $(T[:label_growth_ratio])")
    put(:label_rotation_tol, s.RHO, "rotation numbers differ by at most $(T[:label_rotation_tol])")
    put(:label_estimator_bound, s.RHO, "|lambda_L1 - lambda_L5| above $(T[:label_estimator_bound]) gives U")
    put(:h1_u_min, s.H1, "U(b_P) >= $(T[:h1_u_min])")
    put(:h1_u_lower_min, s.H1, "the lower end of its 95 percent bootstrap interval is >= $(T[:h1_u_lower_min])")
    put(:h1_amin_max, s.H1, "A_min(b_P) <= $(T[:h1_amin_max])")
    put(:h2_q_min, s.H2, "H2 is SUPPORTED if Q >= $(T[:h2_q_min])")
    put(:h2_q_lower_gt, s.H2, "the lower end of the interval of Q is > $(T[:h2_q_lower_gt])")
    put(:h3_r_min, s.H3, "R_X >= $(T[:h3_r_min])")
    put(:h3_r_lower_gt, s.H3, "the lower end of the interval of R_X is > $(T[:h3_r_lower_gt])")
    put(:h3_count_min, s.H3, "for at least $(T[:h3_count_min]) of the 7 indicators")
    put(:h4_dphi_min, s.H4, "Delta_phi >= $(T[:h4_dphi_min])")
    put(:h4_dphi_lower_min, s.H4, "the lower end of its interval is >= $(T[:h4_dphi_lower_min])")
    put(:h5_fmed_min, s.H5, "F_med >= $(T[:h5_fmed_min])")
    put(:h5_fmed_lower_min, s.H5, "the lower end of the interval is >= $(T[:h5_fmed_lower_min])")
    put(:h6_flip_min, s.H6, "max_X delta_dir(X) >= $(T[:h6_flip_min])")
    put(:bootstrap_resamples, s.HYP, "Every statistic with a bootstrap interval uses $(T[:bootstrap_resamples]) resamples of orbits")
    put(:bootstrap_level, s.HYP, "with a $(T[:bootstrap_level]) level interval")
    put(:bootstrap_seed, s.OPEN, "the bootstrap settings ($(T[:bootstrap_resamples]) resamples, seed $(T[:bootstrap_seed]))")
    put(:min_set, s.M8, "A set with fewer than $(T[:min_set]) orbits is not decidable")
    put(:range_lower_percentile, s.M10, "the interval from the $(T[:range_lower_percentile]) to the $(T[:range_upper_percentile]) percentile")
    put(:range_upper_percentile, s.M10, "the interval from the $(T[:range_lower_percentile]) to the $(T[:range_upper_percentile]) percentile")
    put(:gate_max_share, s.PLAN, "chaotic must be at most $(T[:gate_max_share])")
    put(:panel_sha256, s.PANEL, "SHA-256 $(T[:panel_sha256])")
    put(:plan_sha256, s.A2, "plan hash $(T[:plan_sha256])")
    return q
end

# TOML text of the parameter table (one table per key: value, section, and machine or quote and text)
function p2_write_params(path::AbstractString, p::NamedTuple; drop = Symbol[], extra = Dict{String, Any}(), edit = Dict{Symbol, Any}(), style = :section, texts = Dict{Symbol, String}())
    q = p2_quotes(p; style, texts)
    d = Dict{String, Any}()
    for k in P2_PARAM_KEYS
        k in drop && continue
        v = getfield(p, k)
        (k in P2_INT_KEYS && v isa Real) && (v = Int(v))
        e = Dict{String, Any}("value" => v, "section" => q[k][1])
        q[k][2] === :machine ? (e["machine"] = true) : (e["quote"] = q[k][2]; e["text"] = q[k][3])
        haskey(edit, k) && for (kk, vv) in edit[k]; vv === nothing ? delete!(e, kk) : (e[kk] = vv); end        # a value `nothing` removes the field
        d[string(k)] = e
    end
    merge!(d, extra)
    open(io -> TOML.print(io, d), path, "w")
    return path
end

# ---------------------------------------------------------------- the grid starts (Amendment A1) and the fixture preregistration
function p2_h4_tables(p)
    out = NamedTuple[]
    id = 21
    for (ki, kt) in enumerate(p.grid_k_values), i in 0:(p.grid_i_count - 1), j in 0:(p.grid_j_count - 1)
        x0 = 2 * pi * (i + 0.5) / p.grid_i_count
        p0 = pi * (j + 0.5) / p.grid_j_count
        push!(out, (id = id, K = kt, x0 = string(x0), p0 = string(p0), label = "grid-k$(ki)-i$(lpad(i, 2, '0'))-j$(lpad(j, 2, '0'))", kindex = ki, i = i, j = j))
        id += 1
    end
    return out
end

function p2_run_block(p)
    wl = "[" * join(p.windows, ", ") * "]"
    N1, Nc, Nd = p.N_1, p.N_c, p.N_dir
    dirs = "[" * join(["\"$d\"" for d in p.directions], ", ") * "]"
    spec(id, kind, extra, f64 = 1000, d64 = nothing) = "[[run.indicator]]\nid = \"$id\"\nkind = \"$kind\"\n$extra\n[run.indicator.budget.Float64]\nns_cpu_thread = $f64\nsafety = 2\n" *
        (d64 === nothing ? "" : "[run.indicator.budget.Double64]\nns_cpu_thread = $d64\nsafety = 3\n") * "\n"
    freq(id, obs, mr) = spec(id, "frequency", "observable = \"$obs\"\nwindows = $wl\norder = 1\npositions = \"half\"\npad = 8\nmean_removal = $mr\nmirror_ratio = 0.9\nfloor_multiplier = 32.0", 13400,
                             id == "fq1" ? 28100 : nothing)
    arm(name, orbits, ind, prec, hor; methods = nothing, direction = nothing, h1 = false) = "[[run.arm]]\nname = \"$name\"\norbits = \"$orbits\"\n" *
        (methods === nothing ? "" : "methods = [$(join(["\"$m\"" for m in methods], ", "))]\n") * "indicators = [$(join(["\"$i\"" for i in ind], ", "))]\n" *
        "precisions = [\"$prec\"]\nhorizons = [$hor]\n" * (direction === nothing ? "" : "direction = $direction\n") * (h1 ? "h1_check = true\n" : "") * "\n"
    return """
```toml
[run]
panel_prereg_sha256_accepted = ["$(p.panel_build_sha256)"]

[run.primary]
tau = 1
eps = 1e-8
direction = "pi/4"
coordinate_mode = "wrapped"

[run.grid]
points_per_decade = 20
first_steps = 100
extra_fractions = [0.1, 0.5]

[run.backend]
Float64 = "auto"
Double64 = "cpu"

[run.health]
fail_fraction = 0.05
control_bound = 1e-8
control_methods = ["L1"]
control_scope = "exact_start"
setup_check_steps = { T = 100, R = 100, control_elliptic = 100, control_hyperbolic = 5, C = 5, S = 5, X = 5, h4 = 5 }
setup_check_tol = 1e-10
drift_bound = $(p.drift_bound)
rerun_attempts = 1
stop_failed_fraction = 0.05
first_chunk_cost_factor = 2.0
wall_cap_hours = 12

[run.budget.Float64]
ns_cpu_thread = 207.6
ns_gpu = 2.468
safety = 2

[run.budget.Double64]
ns_cpu_thread = 6634.6
safety = 3

[run.k19]
threads = 23

[run.verify]
backend_agreement = { per_orbit_sigma = 6.0, mean_standard_errors = 4.0 }
rerun_sample = 20

[run.pilot]
orbits_per_k = 4
horizon = 10000
seed_label = "p2-indicators-pilot"

$(spec("tg", "tangent", "members = [\"FLI\", \"MEGNO\", \"SALI\"]", 1000, 6000))$(spec("pr", "pair", "members = [\"RLI\", \"MEGNO\"]\nseparation = 1e-12\nframe_source = \"seeded\"", 13000))\
$(spec("st", "stretching", "range = \"derived\"\nbins = 1000\ninterval = 100\nmin_per_cell = 100\ncounts = \"S_D\"", 1100))$(freq("fq1", "exp_i_q", true))\
$(freq("fq1_raw", "exp_i_q", false))$(freq("fq1_qp", "q_plus_i_p", true))$(spec("pa1", "phase_advance", "observable = \"lifted\"\nwindows = $wl\nbirkhoff_order = 1", 1000, 4200))\
$(spec("pa_planar", "phase_advance", "observable = \"planar\"\nwindows = $wl\nbirkhoff_order = 1", 1000))\
$(spec("bk1", "birkhoff", "observables = [\n  { kind = \"cos\", coordinate = 1, frequency = 1 },\n  { kind = \"sin\", coordinate = 1, frequency = 1 },\n  { kind = \"cos\", coordinate = 2, frequency = 1 },\n  { kind = \"sin\", coordinate = 2, frequency = 1 },\n]\nwindows = $wl\nbirkhoff_order = 1", 1000, 5000))\
$(arm("lyap_p", "P", ["tg"], "Float64", N1; methods = ["L1"], h1 = true))$(arm("lyap_g", "h4", ["tg"], "Float64", N1; methods = ["L1", "L5"]))\
$(arm("stretch_p", "P", ["st"], "Float64", Nc))$(arm("stretch_g", "h4", ["st"], "Float64", Nc))$(arm("pair_l", "S_L", ["pr"], "Float64", Nc))\
$(arm("pair_g", "h4", ["pr"], "Float64", Nc))$(arm("win_p", "P", ["fq1", "bk1"], "Float64", Nc))$(arm("win_g", "h4", ["fq1", "bk1"], "Float64", Nc))\
$(arm("phase_p", "P", ["pa1"], "Float64", Nc + 1))$(arm("phase_g", "h4", ["pa1"], "Float64", Nc + 1))\
$(arm("dir_tg", "S_D", ["tg"], "Float64", Nd; methods = ["L1"], direction = dirs))$(arm("sens_win", "S_L", ["fq1_raw", "fq1_qp"], "Float64", Nc))\
$(arm("sens_phase", "S_L", ["pa_planar"], "Float64", Nc + 1))$(arm("d64_lyap", "S_L", ["tg"], "Double64", Nd; methods = ["L1"]))\
$(arm("d64_win", "S_L", ["fq1", "bk1"], "Double64", Nc))$(arm("d64_phase", "S_L", ["pa1"], "Double64", Nc + 1))```
"""
end

# the fixture preregistration text with the SELF hash rule; returns (text, header hash)
function p2_prereg_text(p; slug = P2_SLUG, hashrow = true, drop_sections = Symbol[], extra = "", amend_hash_row = false, style = :section, texts = Dict{Symbol, String}())
    q = p2_quotes(p; style, texts)
    h4 = p2_h4_tables(p)
    body = IOBuffer()
    println(body, "# Preregistration: ", slug, "\n\n## Header\n")
    println(body, "| Field | Value |\n|-------|-------|\n| Study slug | ", slug, " |")
    hashrow && println(body, "| SHA-256 of this file (current version) | ", "0"^64, " (rule: fixture) |")
    println(body, "| Amendment log | A1 (2026-10-29, fixture); A2 (2026-10-29, fixture) |\n\nFixture. Not a study.\n")
    for sk in P2_SEC_ORDER
        sk in drop_sections && continue
        println(body, "## ", getfield(P2_SEC, sk), "\n")
        for k in P2_PARAM_KEYS
            (q[k][1] == getfield(P2_SEC, sk) && q[k][2] !== :machine) && println(body, "- ", q[k][2])
        end
        sk === :RUN && println(body, "\n", p2_run_block(p))
        if sk === :A1
            println(body, "\nThe starts of grid G as `[[run.h4]]` tables.\n\n```toml\n[run]\n")
            for g in h4
                println(body, "[[run.h4]]\nK = \"", g.K, "\"\nx0 = \"", g.x0, "\"\np0 = \"", g.p0, "\"\nlabel = \"", g.label, "\"\n")
            end
            println(body, "```")
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
