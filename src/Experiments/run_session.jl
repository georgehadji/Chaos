# `run_campaign`: refusals first (nothing written), then the store scan, the group loop and the final manifest.

_prev_manifest(path) = isfile(path) ? (try read_json(path) catch; nothing end) : nothing

function _prior_state(prev, cap_override)
    prev === nothing && return Any[], Set{String}(), String[]
    paused = Set{String}(n for (n, a) in get(prev, "arms", Dict()) if get(a, "status", "") == "paused: cost cap")
    cap = String.(get(get(prev, "overrides", Dict()), "cap", String[]))
    return Vector{Any}(get(prev, "sessions", Any[])), setdiff(paused, cap_override), cap
end

function _overrides(dirty, over_threads, cap)
    o = Dict{String, Any}()
    dirty === nothing || (o["dirty"] = Dict{String, Any}("paths" => dirty.paths, "diff_sha256" => dirty.diff_sha256))
    over_threads === nothing || (o["threads"] = Dict{String, Any}("used" => over_threads.used, "k19" => over_threads.k19))
    isempty(cap) || (o["cap"] = sort!(collect(String, cap)))
    return o
end

function _open_session(plan, prereg, panel; store, manifest, threads, capacity, free_bytes, hooks, git, dirty, over_threads, rec, prior, cap_override)
    ctx = make_ctx(plan)
    free = free_bytes === nothing ? Int(Sys.free_memory()) : Int(free_bytes)
    cap = capacity === nothing ? _launch_capacity(plan, free, 0.5) : max(1, Int(capacity))
    table = backend_table(plan; threads = threads, record = rec)
    sessions, paused, prior_cap = _prior_state(_prev_manifest(manifest), cap_override)
    capset = union(prior_cap, cap_override)
    return Session(; ctx = ctx, hooks = hooks, store = String(store), manifest = String(manifest), threads = Int(threads), capacity = cap,
                   free_bytes = free, table = table, git = git, overrides = _overrides(dirty, over_threads, capset), records = (rec, prior),
                   sha_start = _input_hashes(prereg, panel), prereg_path = String(prereg), panel_path = String(panel),
                   prev_sessions = sessions, prev_wall = sum((Float64(x["wall_seconds"]) for x in sessions); init = 0.0), paused = paused,
                   cap_override = String.(cap_override))
end

function _begin!(s::Session, stale_mode, retry)
    plan = s.ctx.plan
    keep = Dict{Int, String}()
    if isdir(s.store)
        keep, stale = classify_store(plan, s.store)
        handle_stale!(s, stale, stale_mode)
    end
    mkpath(s.store)
    remove_temporaries!(s)
    load_entries!(s, keep, retry)
    for arm in plan.block.raw["arm"]
        any(e -> arm["name"] in e["arms"] && e["status"] in ("ok", "failed"), s.entries) && push!(s.measured, arm["name"])
    end
    s.t0 = _now(s.hooks, :session_start)
    s.now = s.t0
    clean = s.git.known && isempty(s.git.tracked) && isempty(s.git.untracked)
    s.session = Dict{String, Any}("git_sha" => s.git.sha, "clean" => clean, "start_utc" => _utc(s.t0), "threads_used" => s.threads)
    _refresh_paused!(s)
    write_manifest(s)
end

function _final_status!(s::Session)
    s.status == "running" || return
    cur = _input_hashes(s.prereg_path, s.panel_path)
    s.sha_end = cur
    if cur.prereg != s.sha_start.prereg
        s.status = "stopped: preregistration changed"
    elseif cur.panel != s.sha_start.panel
        s.status = "stopped: panel changed"
    elseif any(e -> e["status"] == "paused", s.entries)
        s.status = "paused arms"
    elseif any(e -> e["status"] == "failed", s.entries)
        s.status = "complete with failed chunks"
    else
        s.status = "complete"
    end
end

"""
    run_campaign(prereg, panel; store, manifest, threads, ...) -> (status, reasons, ran, skipped, failed, paused, stale, manifest, plan)

Run the study campaign: refusals at start throw `ArgumentError` and write nothing; stops and pauses return a status. See the
specification (Amendment D53, R1 to R8) and `RunHooks` for the test hooks.
"""
function run_campaign(prereg::AbstractString, panel::AbstractString; store::AbstractString, manifest::AbstractString,
                      threads::Integer = Threads.nthreads(), thread_override::Bool = false, allow_dirty::Bool = false, repo = nothing,
                      code_version = nothing, defaults = nothing, capacity = nothing, free_bytes = nothing, order::Symbol = :plan,
                      stale::Symbol = :stop, retry_failed::Bool = false, cap_override = String[], device_record = nothing,
                      prior_record = nothing, hooks::RunHooks = RunHooks())
    plan = plan_run(prereg, panel; code_version = code_version, defaults = defaults)
    over_threads = _check_threads(plan, threads, thread_override)
    git, dirty = _check_tree(repo === nothing ? _root() : repo, allow_dirty)
    rec = device_record === nothing ? gpu_record() : device_record
    s = _open_session(plan, prereg, panel; store = store, manifest = manifest, threads = threads, capacity = capacity, free_bytes = free_bytes,
                      hooks = hooks, git = git, dirty = dirty, over_threads = over_threads, rec = rec, prior = prior_record,
                      cap_override = String.(cap_override))
    Threads.nthreads() > 1 && LinearAlgebra.BLAS.set_num_threads(1)
    _begin!(s, stale, retry_failed)
    groups = collect(enumerate(launch_groups(plan; capacity = s.capacity)))
    order === :reverse && reverse!(groups)
    for (gi, g) in groups
        run_group!(s, gi, order === :reverse ? reverse(g.chunks) : g.chunks) && break
    end
    s.now = _now(hooks, :session_end)
    _final_status!(s)
    _refresh_paused!(s)
    write_manifest(s)
    return (status = s.status, reasons = s.reasons, ran = s.ran, skipped = s.skipped, failed = s.failed, paused = sort!(collect(s.paused)),
            stale = s.stale, manifest = s.manifest, plan = plan)
end

"Default store of a preregistration: `<root>/<study slug>`, the root being `<repository>/data`."
function default_store(prereg::AbstractString; root::AbstractString = joinpath(_root(), "data"))
    return joinpath(root, _study_slug(_lf(read(prereg, String))))
end
