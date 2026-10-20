# The command `analyze` (entry point: `bin/chaoslab.jl` routes it here). Exit codes: 0 done, 1 refusal or usage, 2 tables written but figures skipped.

const CLI_USAGE = "usage: chaoslab analyze <preregistration> <manifest> <panel> <parameters> [--store <dir>] [--out <dir>] [--threads N]"
const CLI_VALUE_OPTIONS = ("--store", "--out", "--threads")

function _parse(args)
    pos, opt = String[], Dict{String, String}()
    i = 1
    while i <= length(args)
        a = String(args[i])
        if startswith(a, "--")
            a in CLI_VALUE_OPTIONS || throw(ArgumentError("unknown option $a"))
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

function _run(pos, opt, io)
    length(pos) == 4 || throw(ArgumentError("analyze needs 4 arguments (preregistration, manifest, panel, parameters), got $(length(pos))"))
    prereg, manifest, panel, params = pos
    store = get(() -> Experiments.default_store(prereg), opt, "--store")
    out = get(() -> joinpath(dirname(manifest), "..", "..", "04-analyze", "output"), opt, "--out")
    nt = haskey(opt, "--threads") ? something(tryparse(Int, opt["--threads"]), throw(ArgumentError("--threads needs an integer"))) : Threads.nthreads()
    res = analyze(; prereg, manifest, panel, params, store, out, figures = true, threads = nt)
    res.status === :control_defect && (println(io, "control defect: ", join(["control $(d.id) $(d.method)" for d in res.defects], "; "), "; report ", res.tables); return 1)
    println(io, "tables written: ", res.tables)
    res.status === :figures_skipped && (println(io, "figures skipped: the plotting package is not loaded (run under the verify environment)"); return 2)
    res.status === :figures_failed && (println(io, "figures failed: the tables are written, see the warning above"); return 2)
    println(io, "figures written in ", out)
    return 0
end

"""
    cli(args; io = stdout) -> Int

`analyze <preregistration> <manifest> <panel> <parameters> [--store dir] [--out dir] [--threads N]`, `analyze --help`. Returns the exit code.
"""
function cli(args; io::IO = stdout)
    a = String.(args)
    if isempty(a) || a[1] != "analyze"
        println(io, CLI_USAGE)
        return 1
    end
    ("--help" in a || "-h" in a) && (println(io, CLI_USAGE); return 0)
    try
        pos, opt = _parse(a[2:end])
        return _run(pos, opt, io)
    catch e
        e isa ArgumentError || rethrow()
        println(io, "error: ", e.msg, "\n", CLI_USAGE)
        return 1
    end
end
