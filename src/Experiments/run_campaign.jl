# The campaign (R3, R5, R7, R8): session state, start checks, hooks. The group loop is in run_groups.jl, the manifest in
# run_manifest.jl, the store scan in run_resume.jl, and the entry points `run_campaign`, `run_pilot` and `cli` use them.

import Dates

"""
    RunHooks(; clock, kill, mutate, observe)

Test hooks of the runner; every field is optional and the module calls a hook only through this struct. `clock(event, info)` is the
only source of time (seconds since the Unix epoch); `kill(point, info)` is called at group and writer points and its exceptions
propagate unchanged; `mutate(info, orbits)` runs after a chunk attempt is computed, before the checks; `observe(info, attrs, orbits)`
runs after the checks, just before the write.
"""
Base.@kwdef struct RunHooks
    clock::Any = nothing
    kill::Any = nothing
    mutate::Any = nothing
    observe::Any = nothing
end

const NS_PER_S = 10^9          # nanoseconds per second
const STRICT_SLACK = 2.0^-40   # a count equal to a fraction of the chunk passes whatever the rounding of that product

# Seconds since the Unix epoch at the resolution of the monotonic clock (`time()` ticks in milliseconds on some systems).
const CLOCK_ORIGIN = Ref((0.0, UInt64(0)))
function real_now()
    CLOCK_ORIGIN[][1] == 0.0 && (CLOCK_ORIGIN[] = (time(), time_ns()))
    t0, n0 = CLOCK_ORIGIN[]
    return t0 + (time_ns() - n0) / NS_PER_S
end

_now(h::RunHooks, event, info = nothing) = h.clock === nothing ? real_now() : Float64(h.clock(event, info))
_kill(h::RunHooks, point, info) = (h.kill === nothing || h.kill(point, info); nothing)

_utc(t::Real) = Dates.format(Dates.unix2datetime(t), "yyyy-mm-ddTHH:MM:SS") * "Z"

"Everything a work item needs besides itself: the plan, the grids by (horizon, tau) and the Float64 reference items of the setup check."
struct RunCtx
    plan::RunPlan
    df::Dict{String, Any}
    grids::Dict{Tuple{Int, Int}, Vector{Int}}
    ref::Dict{Tuple, NamedTuple}
    frames::Dict{Tuple{String, String, Int}, Any}
    warmed::Set{Tuple}
end

function make_ctx(plan::RunPlan, items = plan.items)
    gr = plan.block.raw["grid"]
    grids = Dict{Tuple{Int, Int}, Vector{Int}}()
    ref = Dict{Tuple, NamedTuple}()
    for it in items
        get!(grids, (it.horizon, it.tau)) do
            record_grid(it.horizon, it.tau; points_per_decade = gr["points_per_decade"], first_steps = gr["first_steps"],
                        extra_fractions = gr["extra_fractions"])
        end
        _label(it.precision) == "Float64" && get!(ref, (it.orbit, it.method, it.tau, it.direction, it.coordinate_mode), it)
    end
    return RunCtx(plan, plan.defaults, grids, ref, Dict{Tuple{String, String, Int}, Any}(), Set{Tuple}())
end

# ---------------------------------------------------------------- start checks (R7)
function _git(repo, args...)
    return try
        String(read(pipeline(`git -C $repo $args`; stderr = devnull)))
    catch
        nothing
    end
end

"Git state of `repo`: head commit, tracked changes, untracked files that are not ignored, and a hash of the difference."
function git_status(repo::AbstractString)
    sha = _git(repo, "rev-parse", "HEAD")
    sha === nothing && return (sha = "unknown", known = false, tracked = String[], untracked = String[], diff_sha256 = "")
    st = _git(repo, "status", "--porcelain", "--untracked-files=no")
    un = _git(repo, "ls-files", "--others", "--exclude-standard")
    tracked = [String(strip(l[4:end])) for l in split(st === nothing ? "" : st, '\n') if length(l) > 3]
    untracked = [String(l) for l in split(un === nothing ? "" : un, '\n') if !isempty(l)]
    diff = _git(repo, "diff", "HEAD")
    h = bytes2hex(sha256((diff === nothing ? "" : diff) * "\n" * join(sort(untracked), "\n")))
    return (sha = String(strip(sha)), known = true, tracked = tracked, untracked = untracked, diff_sha256 = h)
end

# Refuse on a dirty tree unless allowed; returns the git record and the dirty-tree override (or nothing).
function _check_tree(repo, allow_dirty::Bool)
    g = git_status(repo)
    dirty = !g.known || !isempty(g.tracked) || !isempty(g.untracked)
    if dirty && !allow_dirty
        parts = String[]
        g.known || push!(parts, "the git state cannot be read")
        isempty(g.tracked) || push!(parts, "the working tree has tracked changes: " * join(g.tracked, ", "))
        isempty(g.untracked) || push!(parts, "the working tree has untracked files that are not ignored: " * join(g.untracked, ", "))
        throw(ArgumentError(join(parts, "; ") * " (commit them, or pass allow_dirty = true to record an override)"))
    end
    return g, dirty ? (paths = vcat(g.tracked, g.untracked), diff_sha256 = g.diff_sha256) : nothing
end

# The thread count of the run against [run.k19].threads (R4 item 2) and against the threads Julia was started with (M7): a difference
# is a refusal unless overridden. The k19 mismatch is recorded in the overrides, the used and available counts in the manifest.
function _check_threads(plan::RunPlan, threads::Integer, override::Bool)
    k19 = plan.block.raw["k19"]["threads"]
    differs = threads != k19
    differs && !override && throw(ArgumentError("the thread count $threads differs from [run.k19].threads = $k19 (pass thread_override = true to record an override)"))
    threads == Threads.nthreads() || override || throw(ArgumentError("the thread count $threads differs from Threads.nthreads() = $(Threads.nthreads()): " *
        "start Julia with -t $threads (the runner cannot use more or fewer workers than Julia has), or pass thread_override = true"))
    return differs ? (used = Int(threads), k19 = k19) : nothing
end

# Hashes the run is bound to: the preregistration by the E18 rule and the panel file by its bytes.
_input_hashes(prereg, panel) = (prereg = prereg_hash(read(prereg, String)), panel = bytes2hex(sha256(read(panel))))

# D54: `[<device>."Float64".recorder.<method>]` as a Dict of NamedTuples; anything that is not a table of tables is an empty table.
function _recorder_table(t)
    t isa AbstractDict || return Dict{String, Any}()
    return Dict{String, Any}(String(m) => (threads = get(r, "threads", 0), cpu_ns_per_step = get(r, "cpu_ns_per_step", nothing),
        gpu_ns_per_step = get(r, "gpu_ns_per_step", nothing), gpu_faster = get(r, "gpu_faster", false), date = get(r, "date", ""),
        launch_steps = get(r, "launch_steps", nothing), orbits = get(r, "orbits", nothing), steps = get(r, "steps", nothing))
        for (m, r) in t if r isa AbstractDict)
end

"""
    read_gpu_record(path, device) -> NamedTuple

Float64 record of `device` from a record file in the format of `data/gpu-benchmark.toml`, with its per-method `recorder` table (D54, D58).
A missing or unreadable file, a device that is not in the file and a device without a Float64 table are an `ArgumentError` naming the path or the device.
"""
function read_gpu_record(path::AbstractString, device::AbstractString)
    isfile(path) || throw(ArgumentError("record file $path does not exist"))
    d = try
        TOML.parsefile(path)
    catch e
        throw(ArgumentError("record file $path is not valid TOML: $(first(sprint(showerror, e), 200))"))
    end
    haskey(d, device) || throw(ArgumentError("record file $path has no table for the device $device"))
    r = d[device] isa AbstractDict ? get(d[device], "Float64", nothing) : nothing
    r isa AbstractDict || throw(ArgumentError("record file $path has no Float64 table for the device $device"))
    return (threads = get(r, "threads", 0), cpu_ns_per_step = get(r, "cpu_ns_per_step", nothing), gpu_ns_per_step = get(r, "gpu_ns_per_step", nothing),
            gpu_faster = get(r, "gpu_faster", false), date = get(r, "date", ""), recorder = _recorder_table(get(r, "recorder", nothing)))
end

"GPU benchmark record (data/gpu-benchmark.toml) of the device in force for Float64, with its per-method `recorder` table (D54), or `nothing`."
function gpu_record(root::AbstractString = _root())
    dev = Foundation.probe(Foundation.DEVICE[], Float64).device
    path = joinpath(root, "data", "gpu-benchmark.toml")
    (isempty(dev) || !isfile(path)) && return nothing
    return try
        read_gpu_record(path, dev)
    catch
        nothing
    end
end

# Orbits per launch from the free memory (R3 item 3): the most demanding run type of the plan bounds it.
function _launch_capacity(plan::RunPlan, free::Integer, f)
    best = typemax(Int)
    for c in plan.chunks
        it = plan.items[first(c.items)]
        T = _element_type(_label(it.precision))
        wrapped = it.coordinate_mode == "wrapped"
        cols = 3 + (wrapped ? 2 : 0) + it.vectors + length(diagnostic_names(it.method))
        M = Tangent.bytes_per_orbit(Val(2), Val(it.vectors), T, wrapped ? 4 : 2)
        best = min(best, launch_capacity(free; f = f, M = M, reducer = 256, columns = cols, grid_steps = length(record_grid(it.horizon, it.tau;
            points_per_decade = plan.block.raw["grid"]["points_per_decade"], first_steps = plan.block.raw["grid"]["first_steps"],
            extra_fractions = plan.block.raw["grid"]["extra_fractions"]))))
    end
    return max(1, best)
end

"Status of a finished campaign to the exit status of the command line: 0 only for `complete`."
function exit_status(status::AbstractString)
    status == "complete" && return 0
    status == "complete with failed chunks" && return 2
    startswith(status, "stopped") && return 3
    return 4
end
