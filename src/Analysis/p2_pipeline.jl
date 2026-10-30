# Study p2-indicators (Amendment D89): the pipeline of the Algorithm of the amendment, in the order of its steps. Nothing is written before every check has passed.

const P2_READINGS = [
    "The parameter table holds a quote of the preregistration for every key that the run blocks do not hold, and the field text (the value as the section writes it) must occur in the quote and the quote in the section (P9, A2)",
    "Output of study p1-le-benchmark is unchanged except the header lines module version and git commit (P10)",
    "The figure drawers of study p2-indicators are a follow-up to the figures extension; the command writes the tables and exits with code 2 (figures skipped) (P11)",
    "The P1 labels of the panel orbits are computed in memory from the P1 reference inputs, designated by the options --p1-* or by the study layout (P12)",
    "lambda_ref and lambda_ns of every orbit are read from the P2 wrapped L1 item at N_1 and at n_s (P13)",
    "A grid orbit is identified by the order of the [[run.h4]] tables with the checks of the stored texts (P14)",
    "n(b) is the first recorded step not before b and n'(n) the recorded step nearest to n / 10, a tie going to the larger step; the step used is printed in T2 (P15)",
    "g_SALI = -L (column SALI.L is ln GALI_2 and negative) so that a non-positive denominator reads alike for the four ratio indicators (P16)",
    "WBA is undecided at a window if any of its four rows has the flag zero or a non-finite dig; E1 and E2 are undecided on their flag or a non-finite D (P17)",
    "W2 and S2 may overlap through T_s; the overlap is counted in T1 (P18)",
    "A grid orbit's structure axis uses the lifted first coordinate (wrapped value plus 2 pi times the winding); only the L1 against L5 test applies (P19)",
    "Bootstrap: one set of draws per label shared by all statistics of the label, 5 percent dropped replicates make the interval undefined, intervals for the pooled statistics only (P20)",
    "H3 keeps h3_count_min for eight and for six indicators; a second version is made per violating indicator (P21)",
    "The stretching items of the other directions, the pooled counts and the counts matrices are not read (P22)",
    "An orbit undecided in some direction or precision is dropped from H6 and counted; RLI has no Double64 item (P23)",
    "Packages (T15) are cross-tabulated with rho; no combination rule is applied (P24)",
    "Regular range: percentile of the same budget, strictly above the upper percentile, the chaotic side is the larger score, K values with fewer than min_set regular orbits have no range (P25)",
    "The gate set is not intersected with E_all; orbits undecided at the last budget leave the denominator (P26)",
    "Statistics with an empty denominator are undefined and make their decision not decidable; the number decided is printed with every statistic (P27, P33)",
    "Controls are the orbits with the control flag; K index 0 belongs to the Y5 control only (P28)",
    "AUC per K value, panel and grid pooled and by source, STR included, minus infinity kept for RLI = 0 (P32)",
    "Pooling by K (AUC, regular range, leaving budget, share below the median, the rows by K) is by the binary64 K value, panel and grid together; the Y5 control has no K and is in no such row (E24, E31)",
    "The share 0.05 of dropped bootstrap replicates (E34) and the design class S of T_s^design (E29) are constants of the module, not keys of the 58 (a gap of the spec)",
    "The 5 percent rule of dropped replicates is applied jointly to the statistics of a label, not per statistic; with the real sets it cannot be reached (P20)",
    "In the E8 extension T_s^design is intersected with E8 like the other sets (owner decision, stage 04 checkpoint)",
    "The verdict series of table R2 give their own detection times, flip counts and direction and precision verdicts in the R2 version (E33)"]

# The P1 reference inputs by option or by the study layout next to the P2 inputs (P12).
function _p1_inputs(prereg, manifest, params, p1_prereg, p1_manifest, p1_params, p1_store)
    pre = p1_prereg === nothing ? joinpath(dirname(abspath(prereg)), P1_SLUG * "-preregistration.md") : p1_prereg
    man = p1_manifest === nothing ? joinpath(dirname(abspath(manifest)), P1_SLUG * "-manifest.json") : p1_manifest
    par = p1_params === nothing ? joinpath(dirname(abspath(params)), P1_SLUG * "-analysis-parameters.toml") : p1_params
    isfile(pre) || _fail("the P1 preregistration $pre does not exist (option --p1-preregistration)")
    return (prereg = pre, manifest = man, params = par, store = p1_store === nothing ? default_store(pre) : p1_store)
end

# E36 m and E23 a: the label of every panel orbit as the base pipeline computes it from the P1 inputs, in memory, writing nothing.
function p2_p1_labels(prereg, manifest, panel, params, store)
    pre = _p1_prepare(prereg, manifest, panel, params, store)
    pre.ctrl.ok || _fail("the P1 reference run stops at a control defect: " * join(["control $(d.id) method $(d.method)" for d in pre.ctrl.defects], "; "))
    c = pre.c
    flags, dirs, bf = _flag_map(pre.items), _dirs_map(pre.items, c.p), _bf_label(pre.items)
    labels = Dict{Int, Any}()
    for o in pre.orbits
        r = orbit_record(pre.idx, o, c.p, flags, dirs, bf)
        labels[r.id] = (instability = r.label.instability, structure = r.label.structure, forced = r.label.forced)
    end
    return (labels = labels, hashes = c.hashes, raw = c.raw)
end

# The K index of a grid start is the position of its K text in grid_k_values (E24).
function _check_grid(p, grid)
    length(grid) == p.grid_orbit_count || _fail("the preregistration has $(length(grid)) h4 starts but grid_orbit_count is $(p.grid_orbit_count)")
    for g in grid
        k = findfirst(==(g.K), p.grid_k_values)
        (k !== nothing && k == g.kindex) || _fail("grid orbit $(g.id) ($(g.label)): the K index of the label is not the position of K = $(g.K) in grid_k_values")
    end
end

# The data tree of the repository: an output folder inside it refuses (I12).
function _check_out(out)
    root = lowercase(replace(abspath(joinpath(pkgdir(parentmodule(@__MODULE__))::String, "data")), '\\' => '/'))
    o = lowercase(replace(abspath(out), '\\' => '/'))
    (o == root || startswith(o, root * "/")) && _fail("the output folder $out is inside the data folder; nothing is written into data/")
end

function _deviations(man, p)
    ov = get(man, "overrides", Dict())
    return String["Overrides recorded in the manifest: " * (isempty(ov) ? "none" : join(sort!(collect(keys(ov))), ", ")::String),
                  "The reference class of the grid orbits uses the reduced L1/L5 rule of the preregistration, section Reference class rho(o), and every grid orbit carries the flag ref_weak (E23 b)",
                  "lambda_ref and lambda_ns of the panel orbits come from the P2 wrapped L1 items; the preregistration says they equal the P1 L1 values, which they do not (P13)",
                  "SALI: the growth function is g_SALI = -L; read literally with a negative L every SALI verdict would be undecided (P16)",
                  "An orbit undecided for E1 or WBA at every budget is undetected, with the censored time, and stays in the regular set (P17)",
                  "The budget step n(b) is the first recorded step not before b, so the indicators are read at slightly different steps (n(b) up to 1.13 b) (P15)",
                  "The panel labels are those of the P1 pipeline at the P1 horizon N(o); the L1 item of table R1 is wrapped, the P1 labels were made with lifted L5 items (P29)"]
end

# Everything the tables need, computed before anything is written.
function _p2_compute(ctx, c, lab, byc, threads)
    p = ctx.p
    orbits, grid_rows = p2_orbit_table(ctx, lab.labels)
    S = p2_sets(orbits, p)
    vt = p2_verdict_tables(ctx, orbits, S)
    hyp = p2_hypotheses(ctx, vt, S, threads)
    orb = Dict(o.id => o for o in orbits)
    auc = p2_auc_table(vt, orbits, S, p)
    range, leaving, below = p2_range_tables(vt, orbits, S, p)
    rep = (omissions = ctx.omissions, auc = auc, sensitivity = p2_sensitivity(ctx, S), flag_counts = p2_flag_counts(byc, c.man), controls = p2_controls(ctx, vt),
           k0line = p2_k0line(ctx, S), cost = p2_cost(c.man))
    packages, byrho = p2_packages(vt, orbits, S, p)
    return (orbits = orbits, grid_rows = grid_rows, S = S, vt = vt, hyp = hyp, orb = orb, rep = rep, range = range, leaving = leaving, below = below, packages = packages, byrho = byrho)
end

# The items T1 to T17 and F1 to F6, the second versions.
function _p2_items(ctx, d)
    p, S, vt, h, r = ctx.p, d.S, d.vt, d.hyp, d.rep
    f4 = p2_f4(r.auc)
    its = p2_hypothesis_items(p, d.orbits, S, vt, h.decisions, P2_I7)
    t = merge((T1 = p2_t1(d.orbits, S, r), T2 = p2_t2(d.orbits, S, vt, p)), its,
              (T9 = p2_t9(f4), T10 = p2_stacked(h.r2, P2_I7), T11 = p2_t11(r.sensitivity, p), T12 = p2_t12(h.gate), T13 = p2_stacked(h.extension, P2_I8; robust = h.robust),
               T14 = p2_t14(d.range, d.leaving, d.below, d.orb), T15 = p2_t15(d.byrho, p), T16 = p2_t16(r.flag_counts, r.omissions, r.controls, r.k0line), T17 = p2_t17(r.cost),
               F1 = p2_f1(S, vt, p), F2 = p2_f2(S, vt, p), F3 = p2_f3(S, vt, p), F4 = f4, F5 = p2_f5(d.grid_rows, vt, p), F6 = p2_f6(d.orbits, S, vt, ctx)))
    second = Dict{Symbol, Dict{String, Any}}()
    for (X, v) in h.second
        i7 = [Y for Y in P2_I7 if Y !== X]
        sub = p2_hypothesis_items(p, d.orbits, S, vt, v.decisions, i7)
        second[X] = Dict{String, Any}(string(k) => getfield(sub, k) for k in keys(sub))
        second[X]["T13"] = p2_stacked(v.extension, [Y for Y in P2_I8 if Y !== X])
    end
    return t, second
end

"""
    analyze_p2(; prereg, manifest, panel, params, store, out, figures, threads, p1_prereg, p1_manifest, p1_params, p1_store) -> NamedTuple

The analysis of study `p2-indicators` (Amendment D89). Every input is checked before anything is written; `<out>/p2-indicators-tables.md` is the only file written. Status `:ok`, or
`:figures_skipped` when figures are wanted (no drawer exists for the study, P11).
"""
function analyze_p2(; prereg, manifest, panel, params, store, out, figures::Bool = false, threads::Integer = Threads.nthreads(), p1_prereg = nothing, p1_manifest = nothing,
                    p1_params = nothing, p1_store = nothing)
    _check_out(out)
    d1 = _p1_inputs(prereg, manifest, params, p1_prereg, p1_manifest, p1_params, p1_store)
    c = p2_check_inputs(prereg, manifest, panel, params, store, d1.params)
    p, man = c.p, c.man
    panel_orbits = TOML.parsefile(panel)["orbit"]
    grid = p2_grid_starts(prereg; panel = panel)
    _check_grid(p, grid)
    items = p2_read_store(manifest, store)
    length(items) == man["plan_counts"]["items"] || _fail("$(length(items)) orbit records read, but plan_counts.items of the manifest is $(man["plan_counts"]["items"])")
    byc = p2_check_items(items, c, panel_orbits, grid)
    idx = p2_index(items)
    omissions = p2_check_coverage(idx, man, p2_requirements(p, panel_orbits, grid))
    _chk_flag_counts(byc, man)
    ctx = (p = p, man = man, panel = panel_orbits, grid = grid, items = items, idx = idx, omissions = omissions)
    p2_controls(ctx)
    lab = p2_p1_labels(d1.prereg, d1.manifest, panel, d1.params, d1.store)
    d = _p2_compute(ctx, c, lab, byc, threads)
    its, second = _p2_items(ctx, d)
    digests = Dict{String, String}("F$i" => p2_digest(getproperty(its, Symbol("F$i"))) for i in 1:6)
    h = merge(c.hashes, (p1_prereg_header_sha256 = lab.hashes.prereg_header_sha256, p1_manifest_sha256 = lab.hashes.manifest_sha256, p1_params_sha256 = lab.hashes.params_sha256,
                         chunks = length(man["chunks"]), items = length(items), seed = p.bootstrap_seed, module_version = string(package_version()),
                         git_sha = git_state(pkgdir(parentmodule(@__MODULE__))::String).sha))
    res = (status = figures ? :figures_skipped : :ok, study = P2_SLUG, tables = nothing, header = h, params = p, items = its, second_items = second, digests = digests,
           rho = Dict(o.id => o.rho for o in d.orbits), lambda_ref = Dict(o.id => o.lambda_ref for o in d.orbits), lambda_ns = Dict(o.id => o.lambda_ns for o in d.orbits),
           p1_labels = lab.labels, grid = d.grid_rows, sets = d.S, verdicts = d.vt.verdicts, scores = d.vt.scores, steps = d.vt.steps, detection = d.vt.detection,
           decisions = d.hyp.decisions, extension = d.hyp.extension, r2 = d.hyp.r2, r2_extension = d.hyp.r2_extension, robust = d.hyp.robust, gate = d.hyp.gate, second = d.hyp.second,
           auc = d.rep.auc, range = d.range, leaving = d.leaving, below_median = d.below, sensitivity = d.rep.sensitivity, packages = d.packages, flag_counts = d.rep.flag_counts,
           omissions = omissions, controls_table = d.rep.controls, k0line = d.rep.k0line, cost = d.rep.cost, readings = P2_READINGS, deviations = _deviations(man, p))
    recheck_inputs(c.raw, prereg, manifest, panel, params)
    recheck_inputs(lab.raw, d1.prereg, d1.manifest, d1.params)
    return merge(res, (tables = p2_write_tables(out, res),))
end
