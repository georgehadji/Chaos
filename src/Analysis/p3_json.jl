# Study p3-abstention (Amendment D93): the calibration record of E58 as JSON text. UTF-8, LF line ends, keys in a fixed order (a NamedTuple is an object in its field order, a
# dictionary an object with sorted keys), counts as integers, rationals as the text "a/b" in lowest terms, reals with 12 significant digits, no clock and no absolute path.

_p3_jstr(s::AbstractString) = "\"" * replace(String(s), "\\" => "\\\\", "\"" => "\\\"", "\n" => "\\n", "\t" => "\\t", "\r" => "\\r") * "\""

_p3_scalar(x::Nothing) = "null"
_p3_scalar(x::Bool) = string(x)
_p3_scalar(x::Integer) = string(x)
_p3_scalar(x::Rational) = "\"$(numerator(x))/$(denominator(x))\""
_p3_scalar(x::AbstractFloat) = isfinite(x) ? string(round(x; sigdigits = 12)) : "null"
_p3_scalar(x::Union{AbstractString, Symbol}) = _p3_jstr(string(x))

_p3_isscalar(x) = x isa Union{Nothing, Bool, Integer, Rational, AbstractFloat, AbstractString, Symbol}

# A key of a dictionary as text: a tuple of keys is joined with `|`.
_p3_jkey(k::Tuple) = join(string.(k), "|")
_p3_jkey(k) = string(k)

# The text of `x` indented by `ind` levels; short lists of scalars stay on one line.
function _p3_json(io::IO, x, ind::Int = 0)
    pad(n) = "  "^n
    if _p3_isscalar(x)
        print(io, _p3_scalar(x))
    elseif x isa Union{AbstractVector, Tuple}
        if all(_p3_isscalar, x) && length(x) <= 24
            print(io, "[", join((_p3_scalar(v) for v in x), ", "), "]")
        else
            isempty(x) && return print(io, "[]")
            println(io, "[")
            for (i, v) in enumerate(x)
                print(io, pad(ind + 1)); _p3_json(io, v, ind + 1); println(io, i < length(x) ? "," : "")
            end
            print(io, pad(ind), "]")
        end
    else
        pairs_ = x isa NamedTuple ? collect(pairs(x)) : sort!(collect(pairs(x)); by = kv -> _p3_jkey(kv[1]))
        isempty(pairs_) && return print(io, "{}")
        println(io, "{")
        for (i, (k, v)) in enumerate(pairs_)
            print(io, pad(ind + 1), _p3_jstr(_p3_jkey(k)), ": "); _p3_json(io, v, ind + 1); println(io, i < length(pairs_) ? "," : "")
        end
        print(io, pad(ind), "}")
    end
    return nothing
end

"The JSON text of an ordered record, ending with a line feed."
p3_json(x) = (io = IOBuffer(); _p3_json(io, x); println(io); String(take!(io)))
