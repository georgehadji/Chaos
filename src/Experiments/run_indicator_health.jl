# Health of the indicator records (Amendment D86, R18, E31): a non-finite value fails the chunk unless the flag of its row says the quantity is undefined;
# flags are data, never failures; `complete` false and the frequency status `failed` fail it. The module draws no conclusion from any other flag (I25).

_bit_set(flags, bits, row, name) = (i = findfirst(==(name), bits); i !== nothing && (flags[row] >> (i - 1)) & 1 == 1)

# A NaN in `col` of `row` is allowed by the flags of the row (E31); an infinity never.
function _ind_nan_ok(bits, flags, row, col)
    "degenerate" in bits && return _bit_set(flags, bits, row, "degenerate")
    "undefined window" in bits && return _bit_set(flags, bits, row, "undefined window")
    zero = _bit_set(flags, bits, row, "zero")
    col == "reldig" && return zero || _bit_set(flags, bits, row, "reldig_undefined")
    return zero && col in ("absdig", "dig")
end

"Name of the first required-finite quantity of an indicator table that is not finite (E31), empty if none."
function ind_nonfinite(id::AbstractString, t::NamedTuple)
    base = "indicators/$id/"
    t.series === nothing || all(isfinite, t.series) || return base * "series"
    t.scalars === nothing || all(isfinite, t.scalars) || return base * "scalars"
    t.final === nothing && return ""
    for row in axes(t.final, 1), (j, col) in enumerate(t.fcolumns)
        v = t.final[row, j]
        isfinite(v) && continue
        (isnan(v) && _ind_nan_ok(t.bits, t.flags, row, col)) || return base * "final column $col (row $row, $(t.rows[row]))"
    end
    return ""
end

"Defects of an indicator record other than a non-finite value: `complete` false, a failed frequency window (R18 item 1 (c), (d))."
function ind_defects(o::OrbitRecord)
    out = String[]
    for (id, t) in o.indicators
        j = t.scolumns === nothing ? nothing : findfirst(==("complete"), t.scolumns)
        j === nothing || t.scalars[j] == 1.0 || push!(out, "indicators/$id: complete = false at the run end")
        t.flags === nothing || !("failed" in t.bits) || !any(r -> _bit_set(t.flags, t.bits, r, "failed"), eachindex(t.flags)) ||
            push!(out, "indicators/$id: the frequency status is failed")
    end
    return out
end
