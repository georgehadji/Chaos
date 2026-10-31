# The plan with rounds (Amendment D94: E39, E41, R25, R26): which rounds and arms the log releases, the orbit subsets D and D_B,
# the chunk cells with the round, the digest of the ref arm and the requests a run may make.

# An arm without a release key in a block with release keys is planned, listed and counted in the ref digest as each_round (M2: one rule).
_release(a) = get(a, "release", "each_round")

"Round of a chunk cell `k<i>:<class>:r<round>`; the controls (no round in the name) count in round 1."
_cell_round(cell) = (m = match(r":r(\d+)$", cell); m === nothing ? 1 : parse(Int, m.captures[1]))

# Subsets D and D_B need the flag `in_D` of the panel, whatever the log says.
function _check_d_columns(rb, pnl)
    any(o -> haskey(o, "in_D"), pnl.orbits) && return
    for a in rb.raw["arm"]
        a["orbits"] in ("D", "D_B") &&
            throw(ArgumentError("run block `arm.orbits`: arm `$(a["name"])` uses the subset $(a["orbits"]) but the panel has no `in_D` flag"))
    end
end

# Orbit ids of the subsets of one round: controls and h4 starts (no round) belong to round 1.
function _round_subsets(bases, r, dbset)
    pick(f) = [b.id for b in bases if (b.round == r || (b.round == 0 && r == 1)) && f(b)]
    return Dict{String, Vector{Int}}("P" => pick(b -> b.design_class != "h4"), "S_D" => pick(b -> b.in_SD), "S_L" => pick(b -> b.in_SL),
                                     "h4" => pick(b -> b.design_class == "h4"), "D" => pick(b -> b.in_D), "D_B" => pick(b -> b.id in dbset))
end

# One pass of the item builder per release phase and round: the each_round arms for the rounds 1 to R, then the closure arms when closed (E41).
function _round_passes(rb, bases, R, closed, dbset)
    passes = NamedTuple[]
    for (phase, on) in (("each_round", true), ("closure", closed)), r in 1:R
        on || continue
        push!(passes, (arms = Set(a["name"] for a in rb.raw["arm"] if _release(a) == phase), subsets = _round_subsets(bases, r, dbset)))
    end
    return passes
end

# D_R, D_B and the allocation of D_B per K index at closure (E38); `R` the number of rounds run.
function _subsets_at(bases, R, closed, db_size)
    DR = [(id = b.id, k_index = b.k_index) for b in bases if b.in_D && 1 <= b.round <= R]
    isempty(bases) && return (DR = DR, dbids = Int[], alloc = Int[])
    i0 = findfirst(b -> b.round > 0, bases)
    i0 === nothing && throw(ArgumentError("run block `arm.release`: the panel has no orbit with a round (a panel with rounds is needed for release keys)"))
    seed = parse(BigInt, bases[i0].run_seed)
    nK = maximum(b.k_index for b in bases)
    alloc = closed ? db_allocation([count(x -> x.k_index == i, DR) for i in 1:nK], db_size) : Int[]
    return (DR = DR, dbids = closed ? db_subset(seed, R, DR, db_size) : Int[], alloc = alloc)
end

"""
Validate a block with release keys against the panel and the log and prepare the passes; `nothing` for a block without release keys (the
log is then ignored). With `check` and no log the plan is a dry plan that holds no round (R25 item 2).
"""
function _round_setup(rb, pnl, bases, round_log, check)
    _check_d_columns(rb, pnl)
    arms = rb.raw["arm"]
    r_max = Int(get(pnl.man, "r_max", 1))
    if !any(a -> haskey(a, "release"), arms)
        for a in arms
            a["orbits"] == "D_B" && throw(ArgumentError("run block `arm.release`: arm `$(a["name"])` uses the subset D_B, which needs a release key (closure)"))
        end
        r_max > 1 && throw(ArgumentError("run block `arm.release`: the panel has $r_max rounds (manifest r_max), so the arms need release keys (R25 item 2)"))
        round_log === nothing || throw(ArgumentError("round log $round_log is given but the block has no release key: nothing would be released by it"))
        return nothing
    end
    for (i, a) in enumerate(arms)
        a["orbits"] == "D_B" && get(a, "release", "") != "closure" &&
            throw(ArgumentError("run block `arm[$i].release`: arm `$(a["name"])` uses the subset D_B and needs release = \"closure\" explicitly"))
    end
    r_max > 1 || throw(ArgumentError("run block `arm.release`: the panel has one round (manifest r_max = $r_max), so a release key has no meaning"))
    bases = [b.round > 0 ? merge(b, (cell = b.cell * ":r$(b.round)",)) : b for b in bases]
    if round_log === nothing
        check || throw(ArgumentError("run block `arm.release`: the arms with a release key need a round log to plan a run (round_log, or --round-log)"))
        return (bases = bases, passes = NamedTuple[], dry = true, R_plan = 0, closed = false, log = (sha256 = _EMPTY_SHA, entries = NamedTuple[]),
                r_max = r_max, sub = _subsets_at(bases, 0, false, 0))
    end
    log = read_round_log(round_log)
    st = round_state(log.entries, r_max)
    sub = _subsets_at(bases, st.R_plan, st.closed, Int(get(pnl.man, "db_size", 0)))
    return (bases = bases, passes = _round_passes(rb, bases, st.R_plan, st.closed, Set(sub.dbids)), dry = false, R_plan = st.R_plan,
            closed = st.closed, log = log, r_max = r_max, sub = sub)
end

# The digest of the ref arm for the rounds 1 to `r` (E40): SHA-256 of the sorted chunk hashes of its chunks with a round in the name, each followed by a newline.
function _ref_digest(chunks, rb, r)
    return bytes2hex(sha256(join(sort!([c.hash for c in _ref_chunks(chunks, rb, r)]) .* "\n")))
end

# The chunks of the ref arms (each_round) with a round up to `r` in the cell name.
function _ref_chunks(chunks, rb, r)
    ref = Set(a["name"] for a in rb.raw["arm"] if _release(a) == "each_round")
    return [c for c in chunks if any(in(ref), c.arms) && occursin(r":r\d+$", c.cell) && _cell_round(c.cell) <= r]
end

# H2: every log entry n says rounds 1 to n are done: a ref chunk of such a round must be in the store, complete, with the plan's chunk hash
# (root attributes only, E43). A refusal names the first entry that is not backed; nothing is written.
function _check_store(plan::RunPlan, store, round_log, hooks)
    rd = plan.rounds
    (rd === nothing || rd.log_entries == 0) && return
    for c in sort!(_ref_chunks(plan.chunks, plan.block, rd.log_entries); by = c -> (_cell_round(c.cell), c.ordinal))
        f = chunk_path(store, c.ordinal)
        a = Dict{String, Any}()
        if isfile(f)
            _read_hook(hooks, :attrs, f)
            a = try read_chunk_attrs(f) catch; Dict{String, Any}() end
        end
        why = !isfile(f) ? "is missing in the store" : get(a, "complete", false) !== true ? "is incomplete" :
              get(a, "chunk_hash", "") != c.hash ? "has a chunk hash different from the plan's" : continue
        _log_err(round_log, _cell_round(c.cell), "the store file $(basename(f)) (round $(_cell_round(c.cell)) chunk of the ref arm) $why")
    end
end

"Check the log against the plan's chunks (R26 item 1) and return the `rounds` record of the plan."
function _round_finish(rs, rb, pnl, chk, chunks, round_log)
    st = rs.dry ? (R_plan = 0, closed = false) :
         check_round_log(rs.log.entries; path = round_log, r_max = rs.r_max, panel = pnl.sha256, prereg = chk.header_sha256,
                         digest = r -> _ref_digest(chunks, rb, r))
    return (R_plan = st.R_plan, closed = st.closed, log_sha256 = rs.log.sha256, log_entries = length(rs.log.entries), D_size = length(rs.sub.DR),
            DB_size = length(rs.sub.dbids), DB_allocation = rs.sub.alloc, release = Dict{String, String}(a["name"] => _release(a) for a in rb.raw["arm"]))
end

# Plan text lines of a plan with rounds: the rounds, the log, and the items and chunks per release phase and round (R25 item 4).
function _rounds_text(rd, chunks, rb)
    io = IOBuffer()
    println(io, "R_plan=", rd.R_plan, "\nclosed=", rd.closed, "\nlog_sha256=", rd.log_sha256, "\nlog_entries=", rd.log_entries)
    rd.R_plan == 0 && println(io, "no round released (dry plan)")
    ref = Set(a["name"] for a in rb.raw["arm"] if _release(a) == "each_round")
    for phase in ("each_round", "closure"), r in 1:rd.R_plan
        cs = [c for c in chunks if (any(in(ref), c.arms) ? "each_round" : "closure") == phase && _cell_round(c.cell) == r]
        isempty(cs) || println(io, "phase ", phase, " round ", r, " items=", sum(length(c.items) for c in cs), " chunks=", length(cs))
    end
    return String(take!(io))
end

"Manifest record of a plan with rounds (R27)."
_rounds_manifest(rd) = Dict{String, Any}("R_plan" => rd.R_plan, "closed" => rd.closed, "log_entries" => rd.log_entries, "log_sha256" => rd.log_sha256,
    "D_size" => rd.D_size, "DB_size" => rd.DB_size, "DB_allocation" => rd.DB_allocation, "release" => rd.release)

# What a run may ask for (R26 item 2): a closure arm needs a closed log, a round beyond R_plan is refused. `arms` and `up_to` select the chunks to run.
function _check_request(plan::RunPlan, arms, up_to)
    names = [a["name"] for a in plan.block.raw["arm"]]
    arms === nothing || isempty(setdiff(String.(arms), names)) || throw(ArgumentError("arms names no arm of the plan: $(join(setdiff(String.(arms), names), ", "))"))
    rd = plan.rounds
    rd === nothing && return
    for a in plan.block.raw["arm"]
        (arms !== nothing && a["name"] in String.(arms)) && _release(a) == "closure" && !rd.closed &&
            throw(ArgumentError("arm `$(a["name"])` is a closure arm and the round log is not closed: it is not released"))
    end
    up_to === nothing || up_to <= rd.R_plan || throw(ArgumentError("a run for round $up_to is refused: the round log releases the rounds 1 to $(rd.R_plan)"))
end

function _selected(plan::RunPlan, arms, up_to)
    (arms === nothing && up_to === nothing) && return nothing
    want = arms === nothing ? nothing : Set(String.(arms))
    return Set(c.ordinal for c in plan.chunks if (want === nothing || any(in(want), c.arms)) && (up_to === nothing || _cell_round(c.cell) <= up_to))
end
