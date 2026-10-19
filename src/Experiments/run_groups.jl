# The launch-group loop (R3, R5): compute a group, finish its chunks one by one (checks, one rerun on CPU, write), then apply the
# campaign stop rules at the group boundary.

# Root attributes of a chunk file (schema/run-store.md, version 2).
function _impl(method)
    method == "L1" && return ("none", "benettin")
    method == "L2" && return ("none", "two-trajectory")
    method == "L9" && return ("RGS", "svd")
    method == "L6" && return ("V-H", "library")
    return (Lyapunov.orthonormalizer_name(_alg(method)), "kernel")
end

function _attrs(plan, it, row, ordinal, key, hash, arms, norbits; git, threads, capacity, utc, status, attempt, first_reasons, failure_reasons,
                setup, pilot, prov = cpu_prov(row.reason), dirty = false, session = 1, control_stop = false)
    gpu = prov.backend === :gpu
    orth, impl = _impl(it.method)
    a = Dict{String, Any}("schema_version" => 2, "study" => plan.study, "chunk" => ordinal, "chunk_key" => key, "chunk_hash" => hash,
        "plan_sha256" => plan.hash, "panel_sha256" => plan.panel.sha256, "git_sha" => git, "code_version" => plan.code_version,
        "prereg_sha256" => plan.prereg_header, "config_sha256" => block_hash(plan.block), "julia_version" => string(VERSION),
        "packages" => packages_hash(), "cpu" => Sys.CPU_NAME, "threads" => threads, "backend_request" => String(row.request),
        "backend" => String(prov.backend), "device" => prov.device, "backend_fallback_reason" => prov.reason, "precision" => it.precision, "method" => it.method,
        "horizon" => it.horizon, "tau" => it.tau, "vectors" => it.vectors, "coordinate_mode" => it.coordinate_mode,
        "direction" => it.direction, "grid_sha256" => it.grid_sha256, "status" => status, "attempt" => attempt,
        "first_attempt_reasons" => first_reasons, "failure_reasons" => failure_reasons, "pilot" => pilot, "arms" => arms,
        "renorm_interval_tau" => it.tau, "deviation_vectors_k" => it.vectors, "seed" => parse(Int, it.run_seed),
        "chunk_size" => gpu ? prov.chunk_size : capacity, "gpu_chunks" => prov.gpu_chunks, "steps_per_launch" => gpu ? prov.launch_steps : it.horizon,
        "launches_per_chunk" => gpu ? prov.launches : cld(norbits, capacity), "orthonormalizer" => orth,
        "implementation" => impl, "setup_check" => setup, "dirty" => dirty, "session" => session, "control_stop" => control_stop,
        "created_utc" => utc, "complete" => true)
    it.eps === nothing || (a["perturbation"] = it.eps)
    return a
end

# Provenance of a chunk: the group's; the rerun of attempt 2 is on the CPU, after a device attempt for the fifth reason.
function _chunk_prov(s, chunk, row, attempt)
    prov = get(s.provs, chunk.ordinal, cpu_prov(row.reason))
    return attempt >= 2 && prov.backend === :gpu ? cpu_prov(Foundation.FALLBACK_REASONS[5]) : prov
end

function chunk_attrs(s, chunk, status, attempt, first_reasons, failure_reasons, setup; control_stop = false)
    plan = s.ctx.plan
    it = plan.items[first(chunk.items)]
    row = _row(s, it)
    return _attrs(plan, it, row, chunk.ordinal, chunk.key, chunk.hash, chunk.arms, length(chunk.items);
                  git = s.git.sha, threads = s.threads, capacity = s.capacity, utc = _utc(s.now), status = status, attempt = attempt,
                  first_reasons = first_reasons, failure_reasons = failure_reasons, setup = setup, pilot = false,
                  prov = _chunk_prov(s, chunk, row, attempt), dirty = !s.session["clean"], session = s.session_index, control_stop = control_stop)
end

# ---------------------------------------------------------------- computing
# Tasks for a group: a task is worth its start-up cost only with enough orbit-steps to run (an engineering constant, not a study
# value); BigFloat precision is scoped per task only from Julia 1.11 on.
const PAR_STEPS_PER_TASK = 4000
const SMALL_GROUP_STEPS = 50000   # orbit-steps below which the first group of an arm is run once untimed before it is measured

function tasks_for(prec, nt::Integer, steps::Integer)
    (_label(prec) == "BigFloat" && VERSION < v"1.11") && return 1
    return clamp(steps ÷ PAR_STEPS_PER_TASK, 1, max(1, Int(nt)))
end

@noinline function compute_records(ctx::RunCtx, items, nt::Integer; times = nothing)
    recs = Vector{OrbitRecord}(undef, length(items))
    isempty(items) && return recs
    prime_frames!(ctx, items)
    function one(_, j)
        t0 = time_ns()
        recs[j] = run_item(items[j], ctx)
        times === nothing || (times[j] = (time_ns() - t0) / NS_PER_S)
    end
    _par(one, () -> nothing, length(items), nt)
    return recs
end

function _fill_costs!(recs, items, wall, workers, total_steps)
    for (o, it) in zip(recs, items)
        c = o.cost[it.method]
        share = it.horizon / total_steps
        c[3] = wall * share
        c[4] = wall * workers * share
    end
end

# The first attempt's wall and thread seconds carry over to the rerun (same inputs, same orbit).
function _copy_costs!(new, old)
    for (n, o) in zip(new, old), (m, c) in n.cost
        c[3], c[4] = o.cost[m][3], o.cost[m][4]
    end
end

# ---------------------------------------------------------------- one chunk
function _attempt!(s::Session, chunk, attempt, recs)
    plan = s.ctx.plan
    backend = attempt == 1 ? get(s.provs, chunk.ordinal, cpu_prov()).backend : :cpu
    s.hooks.mutate === nothing || s.hooks.mutate((ordinal = chunk.ordinal, attempt = attempt, backend = backend, chunk = chunk,
                                                  items = plan.items[chunk.items], arms = chunk.arms), recs)
    return check_chunk!(s.ctx, chunk, recs)
end

# `final_failed = false` is a control-bound stop: the file records it (`control_stop`), no orbit is excluded, a resume refuses it (H1).
function _write_result!(s::Session, chunk, recs, status, attempt, first_reasons, failure_reasons, setup; final_failed = true)
    status == "failed" && final_failed && mark_excluded!(recs)
    attrs = chunk_attrs(s, chunk, status, attempt, first_reasons, failure_reasons, setup; control_stop = !final_failed)
    s.hooks.observe === nothing || s.hooks.observe((ordinal = chunk.ordinal, attempt = attempt), attrs, recs)
    path = chunk_path(s.store, chunk.ordinal)
    ispath(path) && move_aside!(s.store, path, "aside")
    write_chunk(path, attrs, recs; hook = pt -> _kill(s.hooks, pt, (ordinal = chunk.ordinal, path = path, chunks_done = length(s.ran))))
    push!(s.ran, chunk.ordinal)
    status == "failed" && final_failed && push!(s.failed, chunk.ordinal)
    s.entries[chunk.ordinal] = chunk_entry(chunk, attrs, recs)
    _refresh_paused!(s)
    write_manifest(s)
    _kill(s.hooks, :chunk_done, (ordinal = chunk.ordinal, path = path, chunks_done = length(s.ran)))
    return nothing
end

# Finish one chunk of a group; true when the campaign must stop at once (a control over its bound).
function finish_chunk!(s::Session, chunk, recs)
    r1 = _attempt!(s, chunk, 1, recs)
    if !isempty(r1.fatal)
        push!(s.reasons, "control bound exceeded in chunk $(chunk.ordinal): " * join(r1.fatal, "; "))
        push!(s.control_stops, Dict{String, Any}("ordinal" => chunk.ordinal, "code_version" => s.ctx.plan.code_version, "session" => s.session_index,
                                                 "reasons" => r1.fatal, "retried" => false))
        s.status = "stopped: control bound"
        _write_result!(s, chunk, recs, "failed", 1, String[], r1.fatal, r1.setup; final_failed = false)
        return true
    end
    if isempty(r1.reasons)
        _write_result!(s, chunk, recs, "ok", 1, String[], String[], r1.setup)
    elseif s.ctx.plan.block.raw["health"]["rerun_attempts"] >= 1
        items = s.ctx.plan.items[chunk.items]
        recs2 = compute_records(s.ctx, items, tasks_for(first(items).precision, s.threads, sum(it.horizon for it in items)))
        _copy_costs!(recs2, recs)
        r2 = _attempt!(s, chunk, 2, recs2)
        bad = vcat(r2.reasons, r2.fatal)
        _write_result!(s, chunk, recs2, isempty(bad) ? "ok" : "failed", 2, r1.reasons, bad, r2.setup)
    else
        _write_result!(s, chunk, recs, "failed", 1, String[], r1.reasons, r1.setup)
    end
    return false
end

# ---------------------------------------------------------------- one group
_row(s::Session, it) = s.table[(_label(it.precision), it.method)]

function group_info(s::Session, gi, todo)
    plan = s.ctx.plan
    its = [plan.items[i] for o in todo for i in plan.chunks[o].items]
    row = _row(s, first(its))
    armset = Set(a for o in todo for a in plan.chunks[o].arms)
    arms = [a["name"] for a in plan.block.raw["arm"] if a["name"] in armset]
    return (group = gi, chunks = todo, arms = arms, precision = _label(first(its).precision), orbits = length(its), steps = sum(it.horizon for it in its; init = 0),
            backend = row.backend === :gpu && !s.device_lost[] ? :gpu : :cpu)
end

function _refresh_paused!(s::Session)
    for (e, c) in zip(s.entries, s.ctx.plan.chunks)
        e["status"] in ("pending", "paused") && (e["status"] = all(in(s.paused), c.arms) ? "paused" : "pending")
    end
end

runnable(s::Session, ord) = s.entries[ord]["status"] in ("pending", "paused") && !all(in(s.paused), s.ctx.plan.chunks[ord].arms)

# A stop decided before a group: an input changed, or the wall cap is reached (R5 item 4). `nothing` to go on.
function pre_group_stop(s::Session, t)
    h = s.ctx.plan.block.raw["health"]
    cur = _input_hashes(s.prereg_path, s.panel_path)
    s.sha_end = cur
    cur.prereg == s.sha_start.prereg || return "stopped: preregistration changed", "the preregistration hash changed from $(s.sha_start.prereg) to $(cur.prereg)"
    cur.panel == s.sha_start.panel || return "stopped: panel changed", "the panel file hash changed from $(s.sha_start.panel) to $(cur.panel)"
    total = s.prev_wall + (t - s.t0)
    total >= 3600 * h["wall_cap_hours"] && return "stopped: wall cap", "the campaign wall time $(total) s reached the cap of $(h["wall_cap_hours"]) h"
    return nothing
end

# First completed group of an arm against its budget line (R5 item 5): over first_chunk_cost_factor times the line pauses the arm.
function measure_arms!(s::Session, info, wall, workers, backend::Symbol = :cpu)   # the line of the group's own precision
    rb = s.ctx.plan.block
    factor = rb.raw["health"]["first_chunk_cost_factor"]
    cost_ns = NS_PER_S * wall * workers / info.steps
    for arm in rb.raw["arm"]
        name = arm["name"]
        (name in info.arms && !(name in s.measured)) || continue
        push!(s.measured, name)
        name in s.cap_override && continue
        line = budget_line_ns(rb, arm, backend, info.precision)
        if cost_ns > factor * line
            push!(s.paused, name)
            push!(s.reasons, "arm $name paused: measured $(round(cost_ns; digits = 3)) ns per orbit-step above $factor times the budget line $line ns")
        end
    end
end

# The group starts at its clock event, after the checks (they are not part of its cost). A real-clock group this small is dominated
# by cold-start effects and timer jitter: its first measurement is repeated and the shortest wall time taken (same results).
function timed_compute!(s::Session, items, ntask, info)
    small = s.hooks.clock === nothing && info.steps <= SMALL_GROUP_STEPS && any(a -> !(a in s.measured || a in s.cap_override), info.arms)
    best, recs, prov = Inf, OrbitRecord[], cpu_prov()
    for _ in 1:(small ? 5 : 1)
        t = _now(s.hooks, :group_start, info)
        s.now = t
        recs, prov = compute_group(s.ctx, items, ntask, _row(s, first(items)), s.device_lost)
        s.now = _now(s.hooks, :group_end, info)
        best = min(best, s.now - t)
    end
    return recs, best, prov
end

# Run one launch group; true when the campaign stops here.
function run_group!(s::Session, gi, ordinals)
    todo = filter(o -> runnable(s, o), ordinals)
    isempty(todo) && return false
    plan = s.ctx.plan
    info = group_info(s, gi, todo)
    _kill(s.hooks, :group_start, (group = gi, chunks = todo))
    items = NamedTuple[plan.items[i] for o in todo for i in plan.chunks[o].items]
    prime_frames!(s.ctx, items)
    ntask = tasks_for(first(items).precision, s.threads, info.steps)
    warm_up!(s.ctx, items, ntask; dev = info.backend === :gpu)
    t = _now(s.hooks, :group_start, info)
    s.now = t
    stop = pre_group_stop(s, t)
    if stop !== nothing
        s.status = stop[1]
        push!(s.reasons, stop[2])
        return true
    end
    recs, wall, prov = timed_compute!(s, items, ntask, info)
    for o in todo
        s.provs[o] = prov
    end
    workers = prov.backend === :gpu ? 1 : min(ntask, length(items))
    _fill_costs!(recs, items, wall, workers, info.steps)
    npaused = length(s.paused)
    measure_arms!(s, info, wall, workers, prov.backend)
    if length(s.paused) != npaused   # a pause reaches the disk before the first chunk file of its group (M8)
        _refresh_paused!(s)
        write_manifest(s)
    end
    lo = 0
    for o in todo
        n = length(plan.chunks[o].items)
        finish_chunk!(s, plan.chunks[o], recs[(lo + 1):(lo + n)]) && return true
        lo += n
    end
    return _failed_stop!(s)
end

# Strictly more than stop_failed_fraction of the chunks of the plan failed twice: stop, resumable (R5 item 4b).
function _failed_stop!(s::Session)
    h = s.ctx.plan.block.raw["health"]
    M = s.ctx.plan.counts.chunks
    nf = count(e -> e["status"] == "failed", s.entries)
    nf > h["stop_failed_fraction"] * M * (1 + STRICT_SLACK) || return false
    s.status = "stopped: failed chunks"
    push!(s.reasons, "$nf of $M chunks failed twice, more than stop_failed_fraction $(h["stop_failed_fraction"])")
    return true
end
