# Preregistration hash (E18) and the machine-readable panel block (Amendment A1, item 7).

const _HASH_ROW = r"(?m)^\| SHA-256 of this file[^\n]*$"
const _HASH_CELL = r"(?m)^(\| SHA-256 of this file[^|\n]*\|)[^|\n]*\|"

_lf(text::AbstractString) = replace(text, "\r\n" => "\n")

"""
    prereg_hash(text) -> String

SHA-256 (64 hexadecimal digits) of `text` after replacing CRLF by LF and the value cell of the header row
"SHA-256 of this file" by `SELF` (E18).
"""
prereg_hash(text::AbstractString) = bytes2hex(sha256(codeunits(replace(_lf(text), _HASH_CELL => s"\1 SELF |"))))

"""
    check_prereg(path) -> (header_sha256, computed_sha256, amendments)

Compare the hash in the header row with the recomputed one. An `ArgumentError` for a missing or repeated header row, a value
that is not 64 lowercase hexadecimal digits, or any mismatch.
"""
function check_prereg(path::AbstractString)
    text = _lf(read(path, String))
    rows = collect(eachmatch(_HASH_ROW, text))
    length(rows) == 1 || throw(ArgumentError("$path: expected exactly one header row \"SHA-256 of this file\", found $(length(rows))"))
    m = match(r"^\| SHA-256 of this file[^|]*\|[ \t]*([0-9a-f]{64})(?![0-9a-f])", rows[1].match)
    m === nothing && throw(ArgumentError("$path: the header hash is not 64 hexadecimal digits"))
    header, computed = string(m.captures[1]), prereg_hash(text)
    header == computed || throw(ArgumentError("$path: header SHA-256 $header differs from the computed $computed"))
    log = match(r"(?m)^\| Amendment log \|.*$", text)
    amendments = log === nothing ? 0 : count(_ -> true, eachmatch(r"\bA\d+ \(", log.match))
    return (header_sha256 = header, computed_sha256 = computed, amendments = amendments)
end

"""
    read_prereg(path) -> (header_sha256, computed_sha256, amendments, inputs)

[`check_prereg`](@ref) plus the validated [`panel_inputs`](@ref) of the first ```toml block's `[panel]` table.
"""
function read_prereg(path::AbstractString)
    chk = check_prereg(path)
    text = _lf(read(path, String))
    m = match(r"```toml\n(.*?)\n```"s, text)
    m === nothing && throw(ArgumentError("$path: no ```toml panel block"))
    parsed = try
        TOML.parse(string(m.captures[1]))
    catch e
        throw(ArgumentError("$path: the panel block is not valid TOML ($(sprint(showerror, e)))"))
    end
    haskey(parsed, "panel") || throw(ArgumentError("$path: the TOML block has no [panel] table"))
    return (; chk..., inputs = panel_inputs(parsed["panel"]))
end
