# The study manifest (R7): JSON, written atomically after every chunk and at every session end. Chunk entries are rebuilt from the
# chunk files, so a resumed run trusts the files and keeps only the session list of the previous manifest (R8).

Base.@kwdef mutable struct Session
    ctx::RunCtx
    hooks::RunHooks = RunHooks()
    store::String = ""
    manifest::String = ""
    threads::Int = 1
    capacity::Int = 1
    free_bytes::Int = 0
    f::Float64 = 0.5
    table::Dict{Tuple{String, String}, NamedTuple} = Dict{Tuple{String, String}, NamedTuple}()
    git::NamedTuple = (sha = "unknown", known = false, tracked = String[], untracked = String[], diff_sha256 = "")
    overrides::Dict{String, Any} = Dict{String, Any}()
    records::Tuple = (nothing, nothing)
    sha_start::NamedTuple = (prereg = "", panel = "")
    removed::Vector{String} = String[]
    entries::Vector{Dict{String, Any}} = Dict{String, Any}[]
    prev_sessions::Vector{Any} = Any[]
    session::Dict{String, Any} = Dict{String, Any}()
    prereg_path::String = ""
    panel_path::String = ""
    prev_wall::Float64 = 0.0
    t0::Float64 = 0.0
    now::Float64 = 0.0
    paused::Set{String} = Set{String}()
    measured::Set{String} = Set{String}()
    cap_override::Vector{String} = String[]
    ran::Vector{Int} = Int[]
    skipped::Vector{Int} = Int[]
    failed::Vector{Int} = Int[]
    stale::Vector{String} = String[]
    status::String = "running"
    reasons::Vector{String} = String[]
    sha_end::NamedTuple = (prereg = "", panel = "")
    provs::Dict{Int, NamedTuple} = Dict{Int, NamedTuple}()   # provenance of the launch group of each chunk computed in this session
    device_lost::Base.RefValue{Bool} = Ref(false)
    prev::Any = nothing                                       # the manifest of the previous session (or its numbered backup), if any
    control_stops::Vector{Any} = Any[]                        # control-bound stops of all sessions (H1)
    backup_failures::Vector{String} = String[]                # numbered backups skipped in this session (D59)
    session_index::Int = 1                                    # ordinal of this session in the manifest's session list
end

_record_dict(r::Nothing) = "absent"
function _record_dict(r)
    d = Dict{String, Any}("threads" => r.threads, "cpu_ns_per_step" => r.cpu_ns_per_step, "gpu_ns_per_step" => r.gpu_ns_per_step,
                          "gpu_faster" => r.gpu_faster, "date" => r.date)
    hasproperty(r, :recorder) && (d["recorder"] = Dict{String, Any}(m => _record_dict(x) for (m, x) in r.recorder))
    return d
end

function _backend_rows(s::Session)
    return [Dict{String, Any}("precision" => p, "method" => m, "request" => String(r.request), "backend" => String(r.backend),
                              "device" => r.device, "reason" => r.reason, "note" => r.note) for ((p, m), r) in sort!(collect(s.table); by = first)]
end

function _static(s::Session)
    plan = s.ctx.plan
    clean = s.git.known && isempty(s.git.tracked) && isempty(s.git.untracked)
    dev = Foundation.probe(Foundation.DEVICE[], Float64).device
    return Dict{String, Any}("study" => plan.study, "schema_version" => 2, "plan_sha256" => plan.hash,
        "plan_counts" => Dict{String, Any}("items" => plan.counts.items, "chunks" => plan.counts.chunks,
                                           "steps_before" => plan.steps.total_before, "steps_after" => plan.steps.total_after),
        "git" => Dict{String, Any}("sha" => s.git.sha, "clean" => clean),
        "versions" => Dict{String, Any}("language" => string(VERSION), "packages_sha256" => packages_hash(), "hdf5" => string(HDF5.libversion),
                                        "module" => string(package_version())),
        "host" => Dict{String, Any}("cpu" => Sys.CPU_NAME, "cores" => Sys.CPU_THREADS),
        "threads" => Dict{String, Any}("used" => s.threads, "available" => Threads.nthreads()),
        "backend_table" => _backend_rows(s), "device" => isempty(dev) ? "cpu" : dev, "device_record" => _record_dict(s.records[1]),
        "prior_record" => _record_dict(s.records[2]),
        "block" => Dict{String, Any}("text" => block_text(plan.block), "sha256" => block_hash(plan.block)), "defaults" => plan.defaults,
        "launch" => Dict{String, Any}("free_bytes" => s.free_bytes, "f" => s.f, "capacity" => s.capacity),
        "removed_temporaries" => s.removed)
end

"Manifest entry of a chunk from its attributes and orbit records (a file read back, or the buffers just written)."
function chunk_entry(chunk, attrs, recs)
    cost = [c for o in recs for c in values(o.cost)]
    flags = Dict{String, Any}(f => flagged_ids(recs, f) for f in ("diagnostics exceeded", "drift", "excluded", "breakdown_tau"))
    return Dict{String, Any}("ordinal" => chunk.ordinal, "key" => chunk.key, "hash" => chunk.hash, "arms" => chunk.arms,
        "status" => attrs["status"], "attempt" => attrs["attempt"], "backend" => attrs["backend"],
        "fallback_reason" => attrs["backend_fallback_reason"], "orbits" => length(recs), "steps" => round(Int, sum(c[1] for c in cost; init = 0.0)),
        "wall_seconds" => sum(c[3] for c in cost; init = 0.0), "thread_seconds" => sum(c[4] for c in cost; init = 0.0),
        "first_attempt_reasons" => attrs["first_attempt_reasons"], "flags" => flags, "setup_check" => get(attrs, "setup_check", "not applicable"),
        "control_stop" => get(attrs, "control_stop", false) === true, "control_info" => control_info(recs))
end

function pending_entry(chunk, status)
    return Dict{String, Any}("ordinal" => chunk.ordinal, "key" => chunk.key, "hash" => chunk.hash, "arms" => chunk.arms, "status" => status,
        "attempt" => 0, "backend" => "", "fallback_reason" => "", "orbits" => length(chunk.items), "steps" => 0, "wall_seconds" => 0.0,
        "thread_seconds" => 0.0, "first_attempt_reasons" => String[], "flags" => Dict{String, Any}("diagnostics exceeded" => Int[],
        "drift" => Int[], "excluded" => Int[], "breakdown_tau" => Int[]), "setup_check" => "", "control_info" => Any[])
end

# Budget line in ns per orbit-step (R5 item 5): ns of a precision (the arm's first unless given) times its safety factor, for the CPU or the device.
function budget_line_ns(rb::RunBlock, arm, backend::Symbol = :cpu, precision = first(arm["precisions"]))
    b = rb.raw["budget"][precision]
    return (backend === :gpu ? get(b, "ns_gpu", b["ns_cpu_thread"]) : b["ns_cpu_thread"]) * b["safety"]
end

function _arm_entries(s::Session)
    rb = s.ctx.plan.block
    out = Dict{String, Any}()
    for arm in rb.raw["arm"]
        es = [e for e in s.entries if arm["name"] in e["arms"]]
        steps = sum(e["steps"] for e in es; init = 0)
        thr = sum(e["thread_seconds"] for e in es; init = 0.0)
        status = arm["name"] in s.paused ? "paused: cost cap" : all(e -> e["status"] in ("ok", "failed"), es) ? "complete" : "pending"
        out[arm["name"]] = Dict{String, Any}("steps" => steps, "wall_seconds" => sum(e["wall_seconds"] for e in es; init = 0.0),
            "thread_seconds" => thr, "cost_ns_per_orbit_step" => steps > 0 ? NS_PER_S * thr / steps : 0.0, "budget_ns" => budget_line_ns(rb, arm),
            "status" => status)
    end
    return out
end

function manifest_dict(s::Session)
    m = _static(s)
    cur = merge(s.session, Dict{String, Any}("end_utc" => _utc(s.now), "wall_seconds" => s.now - s.t0))
    sessions = vcat(s.prev_sessions, Any[cur])
    m["status"] = s.status
    m["reasons"] = s.reasons
    m["sessions"] = sessions
    m["git"] = Dict{String, Any}("sha" => cur["git_sha"], "clean" => cur["clean"])
    m["prereg"] = Dict{String, Any}("sha256_start" => s.sha_start.prereg, "sha256_end" => isempty(s.sha_end.prereg) ? s.sha_start.prereg : s.sha_end.prereg,
                                    "header" => s.ctx.plan.prereg_header)
    m["panel"] = Dict{String, Any}("sha256_start" => s.sha_start.panel, "sha256_end" => isempty(s.sha_end.panel) ? s.sha_start.panel : s.sha_end.panel,
                                   "build_hash" => s.ctx.plan.panel.build_hash)
    m["chunks"] = s.entries
    m["arms"] = _arm_entries(s)
    m["overrides"] = s.overrides
    m["control_stops"] = s.control_stops
    isempty(s.backup_failures) || (m["backup_failures"] = s.backup_failures)
    m["wall"] = Dict{String, Any}("total_seconds" => sum(x["wall_seconds"] for x in sessions), "cap_hours" => s.ctx.plan.block.raw["health"]["wall_cap_hours"])
    return m
end

const MANIFEST_BACKUPS = 3   # numbered backups manifest.json.1 .. .3 (newest first)

# The session clock is read at every write, so a killed session counts its time up to the last manifest (M3).
function write_manifest(s::Session)
    s.now = max(s.now, _now(s.hooks, :manifest_write))
    return _write_json(s.manifest, manifest_dict(s); keep = MANIFEST_BACKUPS, failures = s.backup_failures)
end
