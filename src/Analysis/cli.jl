# The commands `analyze` and `calibrate` (entry point: `bin/chaoslab.jl` routes them here). Exit codes: 0 done, 1 refusal or usage, 2 (analyze) tables written but figures skipped.

const CLI_USAGE = "usage: chaoslab analyze <preregistration> <manifest> <panel> <parameters> [--store <dir>] [--out <dir>] [--threads N]" *
                  " [--p1-store <dir>] [--p1-preregistration <file>] [--p1-manifest <file>] [--p1-parameters <file>]  (the --p1-* options are for study p2-indicators only)"
const CLI_P1_OPTIONS = ("--p1-store", "--p1-preregistration", "--p1-manifest", "--p1-parameters")
const CLI_VALUE_OPTIONS = ("--store", "--out", "--threads", CLI_P1_OPTIONS...)
const CLI_CALIBRATE_OPTIONS = ("--out", "--threads", "--p2-store", "--p2-manifest", "--p2-preregistration", "--p2-parameters", "--panel", "--p1-store", "--p1-manifest", "--p1-preregistration",
                               "--p1-parameters")
const CLI_CALIBRATE_USAGE = "usage: chaoslab calibrate <preregistration> <parameters> [--out <dir>] [--threads N] [--p2-store <dir>] [--p2-manifest <file>] [--p2-preregistration <file>]" *
                            " [--p2-parameters <file>] [--panel <file>] [--p1-store <dir>] [--p1-manifest <file>] [--p1-preregistration <file>] [--p1-parameters <file>]" *
                            "  (study p3-abstention only; the files of the other studies default to the study layout)"

const CLI_COUNT_OPTIONS = ("--store", "--panel", "--round-log", "--p2-manifest", "--p2-preregistration", "--p2-parameters", "--calibration-record", "--threads")
const CLI_P3T_OPTIONS = (CLI_COUNT_OPTIONS..., "--manifest", "--p1-preregistration", "--p1-parameters", "--out")
const CLI_COUNT_USAGE = "usage: chaoslab count-round <preregistration> <parameters> <test-parameters> [--store <dir>] [--panel <file>] [--round-log <file>] [--p2-manifest <file>]" *
                        " [--p2-preregistration <file>] [--p2-parameters <file>] [--calibration-record <file>] [--threads N]  (study p3-abstention only)"
const CLI_P3T_USAGE = "usage: chaoslab analyze <preregistration> <parameters> <test-parameters> [--store <dir>] [--panel <file>] [--round-log <file>] [--manifest <file>] [--p2-manifest <file>]" *
                      " [--p2-preregistration <file>] [--p2-parameters <file>] [--calibration-record <file>] [--p1-preregistration <file>] [--p1-parameters <file>] [--out <dir>] [--threads N]" *
                      "  (study p3-abstention: the test side)"

function _parse(args, valid = CLI_VALUE_OPTIONS)
    pos, opt = String[], Dict{String, String}()
    i = 1
    while i <= length(args)
        a = String(args[i])
        if startswith(a, "--")
            a in valid || throw(ArgumentError("unknown option $a"))
            i < length(args) || throw(ArgumentError("option $a needs a value"))
            opt[a] = String(args[i + 1])
            i += 2
        else
            push!(pos, a)
            i += 1
        end
    end
    return pos, opt
end

# The --p1-* options designate the P1 reference inputs of study p2-indicators; for any other study they are a usage error (E21).
function _p1_keywords(prereg, opt)
    given = [o for o in CLI_P1_OPTIONS if haskey(opt, o)]
    slug = study_slug(prereg)
    isempty(given) || slug == P2_SLUG || throw(ArgumentError("option $(first(given)) is for study $P2_SLUG only; this preregistration is of study $slug"))
    optval(o) = haskey(opt, o) ? opt[o] : nothing
    return (p1_store = optval("--p1-store"), p1_prereg = optval("--p1-preregistration"), p1_manifest = optval("--p1-manifest"), p1_params = optval("--p1-parameters"))
end

# Whether `analyze` is called for the preregistration of study p3-abstention (its first positional argument names a file with that slug).
function _is_p3t(a)
    pos, i = String[], 2
    while i <= length(a)
        known = a[i] in CLI_P3T_OPTIONS || a[i] in CLI_VALUE_OPTIONS      # an option of either form takes the next argument as its value
        startswith(a[i], "--") ? (i += known ? 2 : 1) : (push!(pos, String(a[i])); i += 1)
    end
    isempty(pos) && return false
    return try
        study_slug(pos[1]) == P3_SLUG
    catch
        false
    end
end

# `--threads N` (the thread count does not change a number, I13); the default is the thread count of the session.
function _threads(opt)
    haskey(opt, "--threads") || return Threads.nthreads()
    n = tryparse(Int, opt["--threads"])
    n === nothing && throw(ArgumentError("--threads needs an integer"))
    return n
end

function _run(pos, opt, io)
    length(pos) == 4 || throw(ArgumentError("analyze needs 4 arguments (preregistration, manifest, panel, parameters), got $(length(pos))"))
    prereg, manifest, panel, params = pos
    kw = _p1_keywords(prereg, opt)
    store = get(() -> Experiments.default_store(prereg), opt, "--store")
    out = get(() -> joinpath(dirname(manifest), "..", "..", "04-analyze", "output"), opt, "--out")
    nt = _threads(opt)
    res = analyze(; prereg, manifest, panel, params, store, out, figures = true, threads = nt, kw...)
    res.status === :control_defect && (println(io, "control defect: ", join(["control $(d.id) $(d.method)" for d in res.defects], "; "), "; report ", res.tables); return 1)
    println(io, "tables written: ", res.tables)
    res.status === :figures_skipped && (println(io, "figures skipped: the plotting package is not loaded (run under the verify environment)"); return 2)
    res.status === :figures_failed && (println(io, res.study == P2_SLUG ? "figures failed: $(res.figures_error)" : "figures failed: the tables are written, see the warning above"); return 2)
    println(io, "figures written in ", out)
    return 0
end

# `calibrate <preregistration> <parameters> [options]` (E46): the two positional arguments, the eleven options as keywords; the thread count does not change a number (I31).
function _run_calibrate(pos, opt, io)
    length(pos) == 2 || throw(ArgumentError("calibrate needs 2 arguments (preregistration, parameters), got $(length(pos))"))
    key(o) = Symbol(replace(replace(o[3:end], "preregistration" => "prereg", "parameters" => "params"), '-' => '_'))
    kw = Dict{Symbol, Any}(key(o) => v for (o, v) in opt if o != "--out" && o != "--threads")
    res = calibrate(; prereg = pos[1], params = pos[2], out = get(opt, "--out", nothing), threads = _threads(opt), kw...)
    println(io, "record written: ", res.record, "
tables written: ", res.tables, "
outcome: ", res.outcome)
    return 0
end

# The files of the study layout for the commands of the test side (E59, D93 E46): the siblings of the preregistration.
function _p3t_keywords(prereg, opt)
    stages = joinpath(dirname(abspath(prereg)), "..", "..")
    st(sub, name) = normpath(joinpath(stages, sub, name))
    p2 = joinpath(dirname(abspath(prereg)), P2_SLUG * "-preregistration.md")
    pick(o, default) = get(opt, o, default)
    return (store = pick("--store", default_store(prereg)), panel = pick("--panel", st("02-panel/output", P3_SLUG * "-panel.toml")),
            round_log = pick("--round-log", st("03-run/output", P3_SLUG * "-round-log.txt")), p2_manifest = pick("--p2-manifest", st("03-run/output", P2_SLUG * "-manifest.json")),
            p2_prereg = pick("--p2-preregistration", p2), p2_params = pick("--p2-parameters", st("04-analyze/references", P2_SLUG * "-analysis-parameters.toml")),
            calibration_record = pick("--calibration-record", st("04-analyze/output", P3_SLUG * "-calibration-record.json")), threads = _threads(opt))
end

function _run_count(pos, opt, io)
    length(pos) == 3 || throw(ArgumentError("count-round needs 3 arguments (preregistration, parameters, test-parameters), got $(length(pos))"))
    r = count_round(; prereg = pos[1], params = pos[2], test_params = pos[3], _p3t_keywords(pos[1], opt)...)
    r.action === :closed && println(io, "round log closed: verified, nothing appended")
    r.action === :nothing && println(io, "nothing to append: the round log is verified and has an entry for every finished round")
    r.action === :appended && println(io, "appended the entry of round ", r.round, ":")
    isempty(r.entry) || print(io, r.entry)
    r.counts === nothing || println(io, "counts: ", join(["$k=$v" for (k, v) in pairs(r.counts)], " "))
    println(io, "targets: t_s ", r.targets.t_s, ", w2 ", r.targets.w2, "; decision: ", r.decision)
    return 0
end

function _run_p3t_analyze(pos, opt, io)
    length(pos) == 3 || throw(ArgumentError("analyze (study $P3_SLUG) needs 3 arguments (preregistration, parameters, test-parameters), got $(length(pos))"))
    kw = _p3t_keywords(pos[1], opt)
    stages = joinpath(dirname(abspath(pos[1])), "..", "..")
    p1(o, name) = get(opt, o, normpath(joinpath(dirname(abspath(pos[1])), name)))
    res = analyze(; prereg = pos[1], params = pos[2], test_params = pos[3], kw..., manifest = get(opt, "--manifest", normpath(joinpath(stages, "03-run/output", P3_SLUG * "-manifest.json"))),
                  p1_prereg = p1("--p1-preregistration", P1_SLUG * "-preregistration.md"), p1_params = get(opt, "--p1-parameters", normpath(joinpath(stages, "04-analyze/references", P1_SLUG * "-analysis-parameters.toml"))),
                  out = get(opt, "--out", normpath(joinpath(stages, "04-analyze/output"))))
    for h in ("HA", "HB", "HC")
        println(io, h, ": ", getproperty(res.decisions, Symbol(h)))
    end
    println(io, "results written: ", res.results, "
tables written: ", res.tables)
    return 0
end

"""
    cli(args; io = stdout) -> Int

`analyze <preregistration> <manifest> <panel> <parameters> [--store dir] [--out dir] [--threads N] [--p1-store dir] [--p1-preregistration file] [--p1-manifest file]
[--p1-parameters file]`, `analyze --help`, and `calibrate <preregistration> <parameters> [options]`, `calibrate --help` (Amendment D93). Returns the exit code.
"""
function cli(args; io::IO = stdout)
    a = String.(args)
    if !isempty(a) && a[1] == "calibrate"
        ("--help" in a || "-h" in a) && (println(io, CLI_CALIBRATE_USAGE); return 0)
        try
            pos, opt = _parse(a[2:end], CLI_CALIBRATE_OPTIONS)
            return _run_calibrate(pos, opt, io)
        catch e
            e isa ArgumentError || rethrow()
            println(io, "error: ", e.msg, "
", CLI_CALIBRATE_USAGE)
            return 1
        end
    end
    if !isempty(a) && a[1] == "count-round"
        ("--help" in a || "-h" in a) && (println(io, CLI_COUNT_USAGE); return 0)
        try
            pos, opt = _parse(a[2:end], CLI_COUNT_OPTIONS)
            return _run_count(pos, opt, io)
        catch e
            e isa ArgumentError || rethrow()
            println(io, "error: ", e.msg, "
", CLI_COUNT_USAGE)
            return 1
        end
    end
    if isempty(a) || a[1] != "analyze"
        println(io, CLI_USAGE, "
", CLI_CALIBRATE_USAGE)
        return 1
    end
    ("--help" in a || "-h" in a) && (println(io, _is_p3t(a) ? CLI_P3T_USAGE : CLI_USAGE); return 0)
    try
        if _is_p3t(a)
            pos, opt = _parse(a[2:end], CLI_P3T_OPTIONS)
            return _run_p3t_analyze(pos, opt, io)
        end
        pos, opt = _parse(a[2:end])
        return _run(pos, opt, io)
    catch e
        e isa ArgumentError || rethrow()
        println(io, "error: ", e.msg, "\n", CLI_USAGE)
        return 1
    end
end
