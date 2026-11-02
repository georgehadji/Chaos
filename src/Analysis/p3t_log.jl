# Study p3-abstention (Amendment D95): the round log of the stopping rule (E62 a: the entry format, E63: the verification V1 of D94). An entry is exactly 13 lines in a fixed order; the
# writer and the verifier agree by byte equality of the entry with its canonical re-rendering (I40). Counts are recomputed from the ref datasets at every call; no cache is kept.

const P3T_COUNT_NAMES = ("drift", "excluded", "orbits", "t_s", "target_t_s", "target_w2", "w2")       # ascending byte order
const P3T_EMPTY_SHA = "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855"
const P3T_HASH_LINES = ("panel_sha256", "prereg_sha256", "ref_digest", "prev_sha256")

"""
    p3t_render_entry(e) -> String

The canonical text of a round log entry: `e = (round, decision, counts, panel, prereg, digest, prev)` with `counts` a dictionary by the seven count names; 13 lines, each ending in LF.
"""
function p3t_render_entry(e)
    lines = vcat(["round=$(e.round)", "decision=$(e.decision)"], ["count.$k=$(e.counts[k])" for k in P3T_COUNT_NAMES],
                 ["panel_sha256=$(e.panel)", "prereg_sha256=$(e.prereg)", "ref_digest=$(e.digest)", "prev_sha256=$(e.prev)"])
    return join(lines, "\n") * "\n"
end

# The entries of a log text: each starts at a line `round=`; a text that does not start with such a line is one entry (entry 1, refused by the parser).
function p3t_split_log(text::AbstractString)
    isempty(text) && return String[]
    st = [m.offset for m in eachmatch(r"(?m)^round=", text)]
    (isempty(st) || st[1] != 1) && return [text]
    push!(st, ncodeunits(text) + 1)
    return String[text[st[i]:prevind(text, st[i + 1])] for i in 1:(length(st) - 1)]
end

# The values of an entry (E62 a); an `ArgumentError` naming the entry and the line for a missing, unreadable or non-canonical line.
function _p3t_parse_entry(n, s)
    kv = Dict{String, String}()
    for l in split(s, '\n')
        i = findfirst('=', l)
        i === nothing || (kv[String(l[1:(i - 1)])] = String(l[(i + 1):end]))
    end
    num(k) = (v = tryparse(Int, get(kv, k, "")); v === nothing ? _fail("round log entry $n: line $k is missing or not an integer") : v)
    for k in P3T_HASH_LINES
        occursin(r"^[0-9a-f]{64}$", get(kv, k, "")) || _fail("round log entry $n: line $k is missing or not 64 lowercase hex digits")
    end
    dec = get(kv, "decision", "")
    dec in ("continue", "close") || _fail("round log entry $n: line decision is \"$dec\", not continue or close")
    e = (round = num("round"), decision = dec, counts = Dict(k => num("count.$k") for k in P3T_COUNT_NAMES), panel = kv["panel_sha256"], prereg = kv["prereg_sha256"],
         digest = kv["ref_digest"], prev = kv["prev_sha256"])
    want = p3t_render_entry(e)
    if s != want
        a, b = split(s, '\n'), split(want, '\n')
        i = something(findfirst(j -> j > length(a) || j > length(b) || a[j] != b[j], 1:max(length(a), length(b))), 1)
        _fail("round log entry $n is not in canonical form: line $i is \"$(i <= length(a) ? a[i] : "(missing)")\" where \"$(i <= length(b) ? b[i] : "(end)")\" is expected")
    end
    return e
end

_p3t_sha(s) = bytes2hex(sha256(codeunits(s)))

# Entry n against the rules (b) to (g) of E63; `prev` is the hash of the previous entry's bytes, `closed` whether an entry before it said close.
function _p3t_check_entry(n, e, prev, closed, c, st)
    q = c.q
    closed && _fail("round log entry $n follows an entry that says close")
    e.round == n || _fail("round log entry $n: line round says $(e.round), the rounds must be 1, 2, ... without gap or repeat (expected $n)")
    (n == Int(q.r_max) && e.decision == "continue") && _fail("round log entry $n: line decision says continue at r_max = $(Int(q.r_max))")
    e.prev == prev || _fail("round log entry $n: line prev_sha256 is not the SHA-256 of the bytes of the previous entry")
    e.panel == c.panel_hash || _fail("round log entry $n: line panel_sha256 differs from the SHA-256 of the panel file")
    e.prereg == c.prereg_hash || _fail("round log entry $n: line prereg_sha256 differs from the header hash of the preregistration")
    for r in 1:n
        pr = p3t_round_problems(c.chunks, c.panel_orbits, r)
        isempty(pr) || _fail("round log entry $n: the ref chunks of round $r are not complete: " * join(pr, "; "))
    end
    e.digest == p3t_ref_digest(c.chunks, n) || _fail("round log entry $n: line ref_digest differs from the digest of the ref chunks of the rounds up to $n")
    want = p3t_counts(st.rows, q, c.targets, n)
    for k in P3T_COUNT_NAMES
        e.counts[k] == want[k] || _fail("round log entry $n: line count.$k is $(e.counts[k]), the recomputation from the ref datasets gives $(want[k])")
    end
    want_dec = p3t_decide(e.counts["t_s"], e.counts["w2"], e.counts["target_t_s"], e.counts["target_w2"], n, Int(q.r_max))
    e.decision == want_dec || _fail("round log entry $n: line decision is $(e.decision), the rule of the counts gives $want_dec")
end

"""
    p3t_verify_log(text, c, st) -> (entries, closed)

The verification V1 of E63 for the log `text`: every entry byte-identical to its canonical re-rendering, the rounds, the chain of `prev_sha256`, the hashes of the panel and the preregistration,
the digest of the ref chunks and the seven counts recomputed from the ref state `st`, and the decision. An `ArgumentError` names the entry and the line.
"""
function p3t_verify_log(text, c, st)
    prev, closed, out = P3T_EMPTY_SHA, false, NamedTuple[]
    for (n, s) in enumerate(p3t_split_log(text))
        e = _p3t_parse_entry(n, s)
        _p3t_check_entry(n, e, prev, closed, c, st)
        push!(out, e)
        prev, closed = _p3t_sha(s), e.decision == "close"
    end
    return (entries = out, closed = closed, last_text = isempty(out) ? "" : p3t_split_log(text)[end])
end
