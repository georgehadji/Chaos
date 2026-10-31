# The round log (Amendment D94: E40, R26): an append-only text file written by the study, read here and never written (I37).
# An entry is the lines `round=<r>`, `decision=<continue|close>`, one `count.<name>=<integer>` per count, `panel_sha256`, `prereg_sha256`,
# `ref_digest` and `prev_sha256` (the SHA-256 of the text of the previous entry; of the empty text for the first), each ending in a newline.

const _LOG_KEYS = ("round", "decision", "panel_sha256", "prereg_sha256", "ref_digest", "prev_sha256")
const _EMPTY_SHA = bytes2hex(sha256(""))

_log_err(path, n, msg) = throw(ArgumentError("round log $path: entry $n: $msg"))

# Entry texts: a new entry starts at every line `round=`; text before the first one is not an entry.
function _split_entries(text, path)
    ents = String[]
    for line in eachline(IOBuffer(text); keep = true)
        if startswith(line, "round=")
            push!(ents, line)
        elseif isempty(ents)
            _log_err(path, 1, "a line that is not an entry: $(repr(first(strip(line), 40)))")
        else
            ents[end] *= line
        end
    end
    return ents
end

function _parse_entry(e, n, path)
    kv, counts = Dict{String, String}(), Dict{String, Int}()
    for line in eachline(IOBuffer(e))
        m = match(r"^([A-Za-z0-9_.]+)=(.*)$", line)
        m === nothing && _log_err(path, n, "a line that is not key=value: $(repr(first(strip(line), 40)))")
        k, v = String(m.captures[1]::SubString{String}), String(m.captures[2]::SubString{String})
        haskey(kv, k) && _log_err(path, n, "the key `$k` is repeated")
        if startswith(k, "count.")
            c = tryparse(Int, v)
            c === nothing && _log_err(path, n, "`$k` must be an integer, got $(repr(v))")
            counts[k[7:end]] = c
        else
            k in _LOG_KEYS || _log_err(path, n, "unknown key `$k`")
        end
        kv[k] = v
    end
    for k in _LOG_KEYS
        haskey(kv, k) || _log_err(path, n, "missing key `$k`")
    end
    r = tryparse(Int, kv["round"])
    r === nothing && _log_err(path, n, "`round` must be an integer, got $(repr(kv["round"]))")
    return (text = e, round = r, decision = kv["decision"], counts = counts, panel = kv["panel_sha256"], prereg = kv["prereg_sha256"],
            digest = kv["ref_digest"], prev = kv["prev_sha256"])
end

"""
    read_round_log(path) -> (sha256, entries)

Read the round log at `path`: the SHA-256 of its bytes and the parsed entries. A missing file or a line that is not an entry is an
`ArgumentError` naming the file (and the entry). The file is never written.
"""
function read_round_log(path::AbstractString)
    isfile(path) || throw(ArgumentError("round log $path does not exist"))
    bytes = read(path)
    text = String(copy(bytes))
    isvalid(text) || throw(ArgumentError("round log $path: the file is not valid UTF-8"))
    ents = _split_entries(text, path)
    return (sha256 = bytes2hex(sha256(bytes)), entries = [_parse_entry(e, n, path) for (n, e) in enumerate(ents)])
end

"""
Rounds the log releases (E39) before it is checked: 1 for an empty log, the last round plus 1 after `continue`, the last round after
`close`, at most `r_max`. Returns `(R_plan, closed)`.
"""
function round_state(entries, r_max::Integer)
    isempty(entries) && return (R_plan = 1, closed = false)
    last = entries[end]
    closed = last.decision == "close"
    return (R_plan = clamp(closed ? last.round : last.round + 1, 1, r_max), closed = closed)
end

"""
    check_round_log(entries; path, r_max, panel, prereg, digest) -> (R_plan, closed)

The checks of R26 item 1 for every entry: the rounds 1, 2, ... without a gap or repeat, a decision that is `continue` or `close`, the chain hash,
the panel hash, the preregistration hash, the digest of the ref arm (`digest(r)` gives the digest of the rounds 1 to `r`), no entry after
`close`, no `continue` at `r_max`. A refusal names the entry; nothing is written.
"""
function check_round_log(entries; path, r_max::Integer, panel, prereg, digest)
    prev, closed = _EMPTY_SHA, false
    for (n, e) in enumerate(entries)
        closed && _log_err(path, n, "an entry after a close entry")
        e.round == n || _log_err(path, n, "round $(e.round) where round $n is expected (a gap or a repeated round)")
        e.decision in ("continue", "close") || _log_err(path, n, "the decision $(repr(e.decision)) is neither continue nor close")
        e.prev == prev || _log_err(path, n, "the chain hash is not the SHA-256 of the previous entry")
        e.panel == panel || _log_err(path, n, "the panel hash is not the SHA-256 of the panel file")
        e.prereg == prereg || _log_err(path, n, "the preregistration hash is not the header hash of the preregistration")
        e.digest == digest(n) || _log_err(path, n, "the ref arm digest differs from the digest of the plan's chunk hashes of rounds 1 to $n")
        (e.decision == "continue" && n >= r_max) && _log_err(path, n, "`continue` at round $n, the last round (r_max = $r_max)")
        closed = e.decision == "close"
        prev = bytes2hex(sha256(e.text))
    end
    return round_state(entries, r_max)
end
