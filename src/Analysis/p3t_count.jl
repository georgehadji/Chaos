# Study p3-abstention (Amendment D95): the command `count-round` (E59, E62 c): after round r it verifies the log, counts T_s and W2 over the ref datasets of the rounds 1 to r, and appends one
# entry. It reads the root attributes of every chunk file and the datasets of the ref arm only (I39); it writes the appended bytes of one entry and nothing else (I38).

# The ref state of the rounds `lo` to `hi` merged into `st` (the datasets of those rounds are read here and only here), and the rows rebuilt for the rounds up to `hi`.
function _p3t_extend_state(st, c, lo, hi, on_read)
    items = st === nothing ? Dict{Tuple{Int, String}, Any}() : st.items
    N, d1 = c.p2.N_1, c.p2.directions[1]
    for ch in c.chunks
        (ch.ref && lo <= ch.round <= hi && ch.status == "ok") || continue
        on_read === nothing || on_read(:datasets, ch.path)
        for it in _p2_chunk_items(read_chunk(ch.path), ch.ordinal)
            (it.id in ("L1", "L5") && it.precision == "Float64" && it.horizon == N && it.tau == 1 && it.direction == d1) || continue      # E64 a: selected by content
            haskey(items, (it.orbit, it.id)) && _fail("item stored twice: orbit $(it.orbit), spec $(it.id) (arm ref, Float64, horizon $N, direction $d1), chunk file $(ch.file)")
            items[(it.orbit, it.id)] = it
        end
    end
    rows = [_p3t_ref_row(o, get(items, (Int(o["id"]), "L1"), nothing), get(items, (Int(o["id"]), "L5"), nothing), c.q, c.p2)
            for o in c.panel_orbits if Int(o["round"]) <= hi && Float64(o["K"]) > 0]
    return (rows = rows, items = items)
end

# The entry of round r from the counts (E62 a, b).
function _p3t_new_entry(c, r, counts, prev)
    dec = p3t_decide(counts["t_s"], counts["w2"], c.targets.t_s, c.targets.w2, r, Int(c.q.r_max))
    return (round = r, decision = dec, counts = counts, panel = c.panel_hash, prereg = c.prereg_hash, digest = p3t_ref_digest(c.chunks, r), prev = prev)
end

function _p3t_result(action, e, text, c)
    cnt = NamedTuple{Tuple(Symbol.(P3T_COUNT_NAMES))}(Tuple(e.counts[k] for k in P3T_COUNT_NAMES))
    return (status = :ok, action = action, round = e.round, entry = text, counts = cnt, decision = e.decision, targets = c.targets)
end

# A closed log (or one that has no round to add): the log is verified and nothing is appended.
function _p3t_verify_only(c, text, k, on_read)
    st = k == 0 ? (rows = NamedTuple[], items = Dict{Tuple{Int, String}, Any}()) : _p3t_extend_state(nothing, c, 1, k, on_read)
    v = p3t_verify_log(text, c, st)
    return v, st
end

"""
    count_round(; prereg, params, test_params, store, panel, round_log, p2_manifest, p2_prereg, p2_params, calibration_record, threads = 1, on_read = nothing) -> NamedTuple

The count of the stopping rule after the next round (E62): verifies the log (E63), counts, and appends one entry. Returns `(status, action, round, entry, counts, decision, targets)` with `action`
`:appended`, `:nothing` (no round to add) or `:closed`. An `ArgumentError` names the cause of every refusal; the log then has the bytes it had. `on_read(kind, path)` is called for every read
of a chunk file, `kind` in `:attrs` or `:datasets`.
"""
function count_round(; prereg, params, test_params, store, panel, round_log, p2_manifest, p2_prereg, p2_params, calibration_record, threads::Integer = 1, on_read = nothing)
    p3t_check_not_data(round_log, "the round log")
    c = p3t_context(; prereg, params, test_params, store, panel, p2_manifest, p2_prereg, p2_params, calibration_record, on_read, command = "count-round")
    text = isfile(round_log) ? read(round_log, String) : ""
    isvalid(text) || _fail("the round log $round_log is not valid UTF-8")
    entries = p3t_split_log(text)
    k = length(entries)
    closed = k > 0 && occursin(r"(?m)^decision=close\r?$", entries[end])
    if closed
        v, _ = _p3t_verify_only(c, text, k, on_read)
        v.closed || _fail("round log entry $k: the last entry does not say close")
        p3t_recheck(c.raw)
        return _p3t_result(:closed, v.entries[end], v.last_text, c)
    end
    return _p3t_open(c, text, k, round_log, on_read)
end

# E62 c (ii): no ref chunk of a round above `r` exists in the store; `what` completes the message.
function _p3t_no_round_above(c, r, what)
    above = [x.round for x in c.chunks if x.ref && x.round > r]
    isempty(above) || _fail("a ref chunk of round $(minimum(above)) exists in the store $what")
end

# The log is open: the entries are verified (so that an altered entry is named), the closure arms must not exist, the next round is r = k + 1.
function _p3t_open(c, text, k, round_log, on_read)
    v, st = _p3t_verify_only(c, text, k, on_read)
    v.closed && _fail("round log entry $k says close but the log was read as open")
    clo = [x.file for x in c.chunks if !x.ref]
    isempty(clo) || _fail("a chunk file of a closure arm exists before closure ($(first(clo)) and $(length(clo) - 1) more): the closure arms run after closure of the round log")
    r = k + 1
    r <= Int(c.q.r_max) || _fail("round $r exceeds r_max = $(Int(c.q.r_max))")
    _p3t_no_round_above(c, r, "while the next round to count is round $r: the guard of the runner was bypassed")
    last_e = isempty(v.entries) ? nothing : v.entries[end]
    if !any(x -> x.ref && x.round == r, c.chunks)
        p3t_recheck(c.raw)
        return last_e === nothing ? (status = :ok, action = :nothing, round = 0, entry = "", counts = nothing, decision = "", targets = c.targets) : _p3t_result(:nothing, last_e, v.last_text, c)
    end
    pr = p3t_round_problems(c.chunks, c.panel_orbits, r)
    isempty(pr) || _fail("round $r is not finished: " * join(pr, "; ") * " (round log of the next round not appended)")
    st = _p3t_extend_state(st, c, r, r, on_read)
    e = _p3t_new_entry(c, r, p3t_counts(st.rows, c.q, c.targets, r), isempty(v.entries) ? P3T_EMPTY_SHA : _p3t_sha(v.last_text))
    body = p3t_render_entry(e)
    p3t_recheck(c.raw)
    _p3t_append(round_log, text, body)
    return _p3t_result(:appended, e, body, c)
end

# The log becomes `old * body` by a temporary file in its folder and a rename, after a re-read that finds the log unchanged; on any error the temporary file is removed, the log has its bytes.
function _p3t_append(log, old, body)
    mkpath(dirname(abspath(log)))
    isfile(log) && (iswritable(log) || _fail("the round log $log is not writable"))      # the rename would replace a read-only log
    tmp = log * ".tmp"
    try
        write(tmp, old * body)
        (isfile(log) ? read(log, String) : "") == old || _fail("input file $log changed during the run")
        mv(tmp, log; force = true)
    catch
        rm(tmp; force = true)
        rethrow()
    end
end
