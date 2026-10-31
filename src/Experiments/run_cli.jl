# The command-line entry point (D16, D49, D53): argument parsing only; every command calls the module. `bin/chaoslab.jl` passes ARGS.

const CLI_USAGE = "usage: chaoslab panel <prereg> --out <file> | run <prereg> <panel> [--store <dir>] [--manifest <file>] | pilot <prereg> <panel> --out <dir>; " *
                  "options: --threads N --thread-override --allow-dirty --repo <dir> --code-version <text> --capacity N --cap-override a,b --stale move --retry-failed --retry-control-stop --backend-override --manifest-reset --prior-record <file> --round-log <file>"

const CLI_FLAGS = ("--thread-override", "--allow-dirty", "--retry-failed", "--retry-control-stop", "--backend-override", "--manifest-reset")
const CLI_VALUES = ("--out", "--store", "--manifest", "--threads", "--repo", "--code-version", "--capacity", "--cap-override", "--stale", "--prior-record", "--round-log")

function _cli_parse(args)
    pos, opt = String[], Dict{String, Any}()
    i = 1
    while i <= length(args)
        a = String(args[i])
        if startswith(a, "--")
            a in CLI_FLAGS && (opt[a] = true; i += 1; continue)
            a in CLI_VALUES || throw(ArgumentError("unknown option $a"))
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

_cli_int(opt, key) = haskey(opt, key) ? (n = tryparse(Int, opt[key]); n === nothing ? throw(ArgumentError("$key needs an integer, got $(opt[key])")) : n) : nothing

_device_name() = Foundation.probe(Foundation.DEVICE[], Float64).device

# Keywords shared by `run` and `pilot`: only what the user gave, the module holds the defaults.
function _cli_common(opt; device::AbstractString = _device_name())
    kw = Dict{Symbol, Any}(:thread_override => haskey(opt, "--thread-override"), :allow_dirty => haskey(opt, "--allow-dirty"))
    haskey(opt, "--repo") && (kw[:repo] = opt["--repo"])
    haskey(opt, "--code-version") && (kw[:code_version] = opt["--code-version"])
    (n = _cli_int(opt, "--threads")) === nothing || (kw[:threads] = n)
    (n = _cli_int(opt, "--capacity")) === nothing || (kw[:capacity] = n)
    if haskey(opt, "--prior-record")      # D58: the record of the device in force from an older record file (A4 K19)
        isempty(device) && throw(ArgumentError("--prior-record needs a device in force; run in the verify environment (julia --project=verify), which loads CUDA"))
        kw[:prior_record] = read_gpu_record(opt["--prior-record"], device)
    end
    return kw
end

function _need(pos, n, what)
    length(pos) == n || throw(ArgumentError("$what needs $n arguments, got $(length(pos)); $CLI_USAGE"))
    return pos
end

function _cli_run(pos, opt, io)
    prereg, panel = _need(pos, 2, "run")
    kw = _cli_common(opt)
    kw[:store] = get(() -> default_store(prereg), opt, "--store")
    kw[:manifest] = get(opt, "--manifest") do
        joinpath(_root(), "workspaces", "study", "stages", "03-run", "output", _study_slug(_lf(read(prereg, String))) * "-manifest.json")
    end
    haskey(opt, "--cap-override") && (kw[:cap_override] = String.(split(opt["--cap-override"], ',')))
    haskey(opt, "--stale") && (kw[:stale] = opt["--stale"] == "move" ? :move : throw(ArgumentError("--stale takes the value move")))
    kw[:retry_failed] = haskey(opt, "--retry-failed")
    kw[:retry_control_stop] = haskey(opt, "--retry-control-stop")
    kw[:backend_override] = haskey(opt, "--backend-override")
    kw[:manifest_reset] = haskey(opt, "--manifest-reset")
    haskey(opt, "--round-log") && (kw[:round_log] = opt["--round-log"])
    res = run_campaign(prereg, panel; kw...)
    println(io, "status: ", res.status)
    foreach(r -> println(io, "  ", r), res.reasons)
    println(io, "chunks run: ", length(res.ran), ", skipped: ", length(res.skipped), ", failed twice: ", length(res.failed), "; manifest: ", res.manifest)
    return exit_status(res.status)
end

function _cli_pilot(pos, opt, io)
    prereg, panel = _need(pos, 2, "pilot")
    haskey(opt, "--out") || throw(ArgumentError("pilot needs --out <dir>"))
    res = run_pilot(prereg, panel; out = opt["--out"], _cli_common(opt)...)
    println(io, "pilot report: ", res.path, " (", res.report["not_study_data"], ")")
    return 0
end

function _cli_panel(pos, opt, io)
    prereg, = _need(pos, 1, "panel")
    haskey(opt, "--out") || throw(ArgumentError("panel needs --out <file>"))
    println(io, "panel file: ", build_panel(prereg; out = opt["--out"]))
    return 0
end

"""
    cli(args; io = stdout) -> Int

The commands `panel`, `run` and `pilot` of `bin/chaoslab.jl`. Exit code: 0 complete, 1 refusal or usage error (message on `io`),
2 complete with failed chunks, 3 stopped, 4 paused arms.
"""
function cli(args::AbstractVector{<:AbstractString}; io::IO = stdout)
    if isempty(args) || !(args[1] in ("panel", "run", "pilot"))
        println(io, isempty(args) ? "no command; " : "unknown command $(args[1]); ", CLI_USAGE)
        return 1
    end
    return try
        pos, opt = _cli_parse(args[2:end])
        (args[1] == "run" ? _cli_run : args[1] == "pilot" ? _cli_pilot : _cli_panel)(pos, opt, io)
    catch e
        e isa InterruptException && rethrow()
        println(io, "chaoslab ", args[1], ": ", e isa ArgumentError ? e.msg : sprint(showerror, e))
        1
    end
end
