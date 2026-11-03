# Manifest and pilot additions of the indicator specs (Amendment D86, R19, R20): the table by spec id, the flag counts of a chunk, the omitted orbits, the
# transform seconds, the cost lines per spec. Flags are counted, never read by a health rule (I25).

plan_has_specs(plan::RunPlan) = haskey(plan.block.raw, "indicator")

"Number of orbits of a record list that carry each flag of spec `id` (only the flags that occur)."
function ind_flag_counts(recs, id::AbstractString)
    counts = Dict{String, Int}()
    for o in recs
        t = get(o.indicators, id, nothing)
        (t === nothing || t.flags === nothing) && continue
        for (i, nm) in enumerate(t.bits)
            any(w -> (w >> (i - 1)) & 1 == 1, t.flags) && (counts[nm] = get(counts, nm, 0) + 1)
        end
    end
    return counts
end

function indicator_entry!(e::AbstractDict, chunk, recs, plan)
    id = chunk.run_type.method
    e["indicator_flags"] = ind_flag_counts(recs, id)
    plan !== nothing && block_spec(plan.block, id)["kind"] == "frequency" && (e["transform_seconds"] = sum((o.cost[id][6] for o in recs if length(o.cost[id]) >= 6); init = 0.0))
    if plan !== nothing
        na = [o for (o, cell) in get(plan.omitted_cells, id, Tuple{Int, String}[]) if cell == chunk.cell]
        isempty(na) || (e["not_applicable"] = Dict{String, Any}(id => sort!(na)))
    end
    return e
end

# The table of R19: kind, canonical text, hash, chain record, windows, observables, budget lines and the terms of bytes(orbit) of the first item of the spec.
# The table is a function of the plan alone and `write_manifest` asks for it at every chunk, so the last one is kept (R9 of the D86 review).
const _IND_MANIFEST = Ref{Tuple{UInt, String, Dict{String, Any}}}((UInt(0), "", Dict{String, Any}()))
const _IND_MANIFEST_LOCK = ReentrantLock()
function indicator_manifest(plan::RunPlan)
    key = (objectid(plan), plan.hash)
    lock(_IND_MANIFEST_LOCK) do
        c = _IND_MANIFEST[]
        (c[1], c[2]) == key && return c[3]
        t = _indicator_manifest(plan)
        _IND_MANIFEST[] = (key..., t)
        return t
    end
end

function _indicator_manifest(plan::RunPlan)
    out = Dict{String, Any}()
    for spec in plan.block.raw["indicator"]
        id = spec["id"]
        k = findfirst(i -> i.method == id, plan.items)
        k === nothing && continue   # R28: a spec held only by closure arms has no item before closure and is not listed
        it = plan.items[k]
        b = orbit_bytes(plan, it)
        m = Dict{String, Any}("kind" => spec["kind"], "spec" => spec_text(spec), "spec_sha256" => spec_hash(spec), "chain" => chain_record(spec), "budget" => spec["budget"],
                              "bytes_per_orbit" => Dict{String, Any}("M" => b.M, "R" => b.R, "S" => b.S, "B" => b.B, "total" => b.total))
        ind_window_kind(spec["kind"]) && (m["windows"] = spec["windows"]; m["run_length"] = ind_run_length(spec))
        spec["kind"] == "pair" && (m["frame_source"] = spec["frame_source"]; m["direction_rule"] = IND_DIRECTION_RULE)
        spec["kind"] == "birkhoff" && (m["observable_names"] = ind_observable_names(spec, it))
        out[id] = m
    end
    return out
end

# Measured cost per orbit-step of each spec of an arm against the budget line of the spec (the spec's own, not the block's).
function arm_spec_costs(s, arm, entries)
    plan = s.ctx.plan
    out = Dict{String, Any}()
    for id in arm["indicators"]
        es = [e for e in entries if plan.chunks[e["ordinal"]].run_type.method == id]
        steps, thr = sum((e["steps"] for e in es); init = 0), sum((e["thread_seconds"] for e in es); init = 0.0)
        spec = s.ctx.specs[id]
        prec = first(p for p in arm["precisions"] if haskey(spec["budget"], p))
        out[id] = Dict{String, Any}("measured_ns_per_orbit_step" => steps > 0 ? NS_PER_S * thr / steps : 0.0,
                                    "budget_line_ns" => budget_line_ns(plan.block, arm, :cpu, prec; spec = spec))
    end
    return out
end
