# Included by verify/experiments-indicators.jl (section D86W): D86w (V). The module holds no study value, no threshold and no detection logic.
# There is no p2 preregistration yet, so the values searched for are those of the reduced block of this verification (verify/experiments-indicators-block.toml: chosen for the check,
# not study values) and the key names of a threshold, a class or a detection; the search also lists every assignment of a literal number to a spec key in src/Experiments.
const SRC = joinpath(ROOT, "src")
const SRC_FILES = sort([joinpath(d, f) for (d, _, fs) in walkdir(SRC) for f in fs if endswith(f, ".jl")])
srcrel(p) = replace(relpath(p, ROOT), '\\' => '/')
function grep(rx; files = SRC_FILES, skip_comments = true)
    hits = Tuple{String, Int, String}[]
    for f in files
        for (i, l) in enumerate(eachline(f))
            s = skip_comments ? replace(l, r"#.*$" => "") : l
            occursin(rx, s) && push!(hits, (srcrel(f), i, strip(l)))
        end
    end
    return hits
end
const EXP_FILES = filter(f -> occursin("/src/Experiments/", replace(f, '\\' => '/')), SRC_FILES)
let
    run = block_run()
    # distinctive numbers of the reduced block, as numeric tokens (not part of a longer number, not in a comment or a string)
    vals = Dict("250 (window)" => r"(?<![\w.])250(?![\w.])", "1000 (window)" => r"(?<![\w.])1000(?![\w.])", "2000 / 2001 (run length)" => r"(?<![\w.])200[01](?![\w.])", "4000 (horizon)" => r"(?<![\w.])4000(?![\w.])",
                "1e-12 (separation)" => r"(?<![\w.])1(\.0)?e-12(?![\w])", "8.0 (floor_multiplier)" => r"(?<![\w.])8\.0(?![\w])", "400.0 / 800.0 / 600.0 (budget lines)" => r"(?<![\w.])(400|800|600)\.0(?![\w])")
    for (nm, rx) in sort(collect(vals); by = first)
        h = grep(rx; files = EXP_FILES)
        info("D86w", "src/Experiments, numeric token $nm: $(length(h)) hits" * (isempty(h) ? "" : "; " * join(("$(f):$(i)" for (f, i, _) in first(h, 6)), ", ")))
        for (f, i, l) in first(h, 6)
            say("   hit | $f:$i | $(first(l, 140))")
        end
    end
    h = grep(r"(?<![\w.])(250|1000|2001|4000)(?![\w.])"; files = EXP_FILES)
    h_nonzero = [x for x in h if !occursin(r"^\s*(#|\"\"\")", x[3])]
    explained = [("PILOT_MAX_DRAWS", "limit of the draws of a pilot start (D53, not an indicator value)"), ("WARM_STEPS", "steps of the warm-up run of a run type (D53)"), ("PAR_STEPS_PER_TASK", "steps per task of the parallel scheduler (D53)"),
                 ("_estimator(method, 8.0)", "placeholder orthogonality tolerance of the hot-loop allocation probe of an exponent method (D53)"),
                 ("StretchingNumbers(; range = (-1.0, 1.0), bins = 1, interval = 1, min_per_cell = 0)", "placeholder reducer built only to call check_frame, which reads none of these parameters"),
                 ("det_threshold", "the monitor exceedance threshold of R5 (determinant and symplecticity residual), not an indicator threshold")]
    why(l) = (i = findfirst(e -> occursin(e[1], l), explained); i === nothing ? nothing : explained[i][2])
    h = grep(r"(?<![\w.])(250|1000|2001|4000)(?![\w.])|(?<![\w.])8\.0(?![\w])"; files = EXP_FILES)
    unexplained = [x for x in h if why(x[3]) === nothing]
    emit("D86w", "src/Experiments: numeric tokens 250, 1000, 2001, 4000 and 8.0 outside comments: hits / hits without an explanation as an engineering constant of D53 or a placeholder (each hit is listed above with its line)", "$(length(h)) / $(length(unexplained))", "0 unexplained", pf(isempty(unexplained)))
    foreach(x -> say("   explained | $(x[1]):$(x[2]) | $(why(x[3]))"), [x for x in h if why(x[3]) !== nothing])
    # assignments of a literal number to a spec key in src/Experiments (a default would be one)
    keys_rx = r"\b(bins|interval|min_per_cell|order|pad|mirror_ratio|floor_multiplier|separation|birkhoff_order|windows|window|frequency)\b\s*(=|::\w+\s*=)\s*[0-9\[(]"
    hk = grep(keys_rx; files = EXP_FILES)
    hk_un = [x for x in hk if why(x[3]) === nothing]
    emit("D86w/I27", "src/Experiments: assignments of a literal number to a spec key name (bins, interval, min_per_cell, order, pad, mirror_ratio, floor_multiplier, separation, birkhoff_order, windows): hits / hits that are a default of a spec key", "$(length(hk)) / $(length(hk_un))", "0 defaults", pf(isempty(hk_un)))
    foreach(x -> say("   hit | $(x[1]):$(x[2]) | $(first(x[3], 160)) | " * (why(x[3]) === nothing ? "UNEXPLAINED" : why(x[3]))), hk)
    # key names of a threshold, a class or a detection in the grammar and the recorders of the indicators
    ind_files = filter(f -> occursin("run_indicator", f) || occursin("run_block.jl", f), EXP_FILES)
    ht = grep(r"(?i)threshold|detect|crossing|first_cross|sticky"; files = ind_files, skip_comments = true)
    ht_un = [x for x in ht if why(x[3]) === nothing]
    emit("D86w", "run_indicator*.jl and run_block.jl: code lines (comments excluded) naming threshold, detect, crossing, first_cross, sticky: lines / lines that are not the R5 monitor threshold", "$(length(ht)) / $(length(ht_un))", "0 (the grammar has no such key)", pf(isempty(ht_un)))
    foreach(x -> say("   hit | $(x[1]):$(x[2]) | $(first(x[3], 160)) | " * (why(x[3]) === nothing ? "UNEXPLAINED" : why(x[3]))), ht)
    ht2 = grep(r"(?i)threshold|detect|crossing|first_cross|sticky"; files = filter(f -> !(f in ind_files), EXP_FILES), skip_comments = true)
    info("D86w", "rest of src/Experiments, code lines with those words (the monitor exceedance rule of R5 and the panel builder; not indicator logic): $(length(ht2)) lines in files $(sort(unique(first.(ht2))))")
    # the grammar refuses a key that could hold a threshold
    base = block_run()
    function refused(mod!)
        d = deepcopy(base)
        mod!(d)
        try
            EX.run_block(d)
            return "accepted"
        catch e
            return e isa ArgumentError ? e.msg : sprint(showerror, e)
        end
    end
    tg(d) = first(s for s in d["indicator"] if s["id"] == "tg")
    for key in ("threshold", "class", "detect", "first_crossing")
        msg = refused(d -> tg(d)[key] = 0.5)
        emit("D86w/D86a", "a spec carrying the key `$key` is refused before any computation", first(msg, 110), "refused naming the key", pf(msg != "accepted" && occursin(key, msg)))
    end
    msg = refused(d -> d["health"]["control_methods"] = ["L1", "tg"])
    emit("D86w/D86a", "a spec id in [run.health].control_methods is refused", first(msg, 110), "refused", pf(msg != "accepted" && occursin("tg", msg)))
    msg = refused(d -> delete!(first(s for s in d["indicator"] if s["id"] == "st"), "bins"))
    emit("D86w/I27", "a spec without one of its kind keys (no default in the module) is refused naming the key", first(msg, 110), "refused naming `bins`", pf(msg != "accepted" && occursin("bins", msg)))
    # names of datasets and attributes written to the files (the independent reader lists every name): none holds a threshold or a first-crossing step
    names = unique(sort(String[m.captures[1] for l in eachline(joinpath(SCR, "reader-B.txt")) for m in (match(r"^N\|([^|]+)\|", l),) if m !== nothing]))
    sus = [n for n in names if occursin(r"(?i)thresh|detect|cross|class|sticky", n)]
    info("D86w", "names of root attributes, orbit attributes, provenance keys and datasets in the files of campaign B ($(length(names)) distinct): matching threshold, detect, cross, class, sticky: $sus")
    allowed = Set(["diagnostics@threshold", "orbit@class_label", "orbit@design_class"])
    emit("D86w", "no attribute or dataset of the new kinds names a threshold, a detection or a crossing; the only matches are the monitor exceedance threshold of R5 and the panel class labels (D53)", sus, "subset of $(sort(collect(allowed)))", pf(Set(sus) <= allowed))
end
