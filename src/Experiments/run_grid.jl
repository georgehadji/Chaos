# Recording grid (R2 item 1), deviation frame (R2 item 6), exact closed-form evaluation and the code version (R1 item 3).

# A real as an exact rational of its shortest decimal text ("0.1" is 1 // 10, never the binary64 value).
function _dec_rational(f::Real)
    m = match(r"^(-?)(\d*)\.?(\d*)(?:[eE]([-+]?\d+))?$", string(f))
    m === nothing && throw(ArgumentError("cannot read $f as a decimal"))
    whole, frac, ex = something(m[2], ""), something(m[3], ""), something(m[4], "0")
    digits = whole * frac
    r = parse(BigInt, isempty(digits) ? "0" : digits) // big(10)^length(frac)
    e = parse(Int, ex)
    r = e >= 0 ? r * big(10)^e : r // big(10)^(-e)
    return m[1] == "-" ? -r : r
end

_half_up(v::BigFloat) = Int(floor(v + BigFloat(1) / 2))

"""
    record_grid(N, tau; points_per_decade, first_steps, extra_fractions) -> Vector{Int}

Steps at which a work item of horizon `N` is recorded: every step `1:first_steps`, `round(10^(j / points_per_decade))` (256-bit
evaluation, half up), `round(f N)` for each extra fraction and `N`; each target is moved to the handoff `tau ceil(n / tau)`; sorted
and without duplicates.
"""
function record_grid(N::Integer, tau::Integer; points_per_decade::Integer, first_steps::Integer, extra_fractions)
    targets = Set{Int}(1:min(first_steps, N))
    setprecision(BigFloat, 256) do
        j = 0
        while (v = BigFloat(10)^(BigFloat(j) / points_per_decade)) <= N
            push!(targets, _half_up(v))
            j += 1
        end
    end
    for f in extra_fractions
        push!(targets, floor(Int, _dec_rational(f) * N + 1 // 2))
    end
    push!(targets, Int(N))
    return sort!(unique!([Int(tau) * cld(t, Int(tau)) for t in targets]))
end

"SHA-256 (64 hexadecimal digits) of the canonical text of a recording grid."
grid_hash(grid::AbstractVector{<:Integer}) = bytes2hex(sha256("chaoslab-grid-1\n" * join(grid, ",")))

# ---------------------------------------------------------------- exact evaluation of closed-form texts
const _NUM_LITERAL = r"(?<![\w.])(\d+(?:\.\d+)?(?:[eE][-+]?\d+)?)"

# Evaluate a parsed expression; a numeric literal arrives as `_n("text")` and is parsed as the exact decimal.
function _ev_exact(x)
    x === :pi && return BigFloat(pi)
    if x isa Expr && x.head === :call && x.args[1] === :_n
        return parse(BigFloat, x.args[2])
    elseif x isa Expr && x.head === :call && haskey(_FUNCS, x.args[1])
        return _FUNCS[x.args[1]](map(_ev_exact, x.args[2:end])...)
    elseif x isa Expr && x.head === :block
        return _ev_exact(last(filter(a -> !(a isa LineNumberNode), x.args)))
    end
    throw(ArgumentError("unsupported expression $x"))
end

"Value of a decimal or closed-form text (`\"0.3 + pi/5\"`) at the ambient BigFloat precision, decimals read exactly."
function _eval_exact(s::AbstractString)
    t = replace(replace(String(s), r"(\d)pi" => s"\1*pi"), _NUM_LITERAL => s"_n(\"\1\")")
    e = try
        Meta.parse(t)
    catch
        throw(ArgumentError("cannot parse \"$s\""))
    end
    return _ev_exact(e)
end

# Decimal text of a closed form, kept as written when it already is a decimal.
function _decimal_text(s::AbstractString)
    tryparse(BigFloat, s) === nothing || return String(s)
    return setprecision(() -> string(_eval_exact(s)), BigFloat, 256)
end

_to_type(::Type{T}, x::BigFloat) where {T} = T(x)
_to_type(::Type{BigFloat}, x::BigFloat) = BigFloat(x; precision = precision(BigFloat))

"""
    deviation_frame(T, direction_text, k) -> SMatrix{2,k,T}

First column `(cos a, sin a)` with `a` evaluated at 256 bits from the closed-form text, second column (`k = 2`) its
counterclockwise quarter turn `(-W[2,1], W[1,1])`.
"""
function deviation_frame(::Type{T}, text::AbstractString, k::Integer) where {T}
    k in (1, 2) || throw(ArgumentError("deviation_frame: k must be 1 or 2, got $k"))
    c, s = setprecision(BigFloat, max(256, precision(BigFloat))) do
        a = _eval_exact(text)
        (cos(a), sin(a))
    end
    c, s = _to_type(T, c), _to_type(T, s)
    return k == 1 ? SMatrix{2, 1, T}(c, s) : SMatrix{2, 2, T}(c, s, -s, c)
end

# ---------------------------------------------------------------- code version
function _lock_file(root)
    for f in (joinpath(root, "Manifest-v$(VERSION.major).$(VERSION.minor).toml"), joinpath(root, "Manifest.toml"))
        isfile(f) && return f
    end
    return ""
end

_file_digest(f) = sha256(codeunits(replace(String(read(f)), "\r\n" => "\n")))

"""
    code_version(root = <package root>) -> String

Content identifier (hex) of the files under `root/src` and of the dependency lock (`Manifest-v<major>.<minor>.toml`, else
`Manifest.toml`); no other file enters.
"""
function code_version(root::AbstractString = _root())
    io = IOBuffer()
    files = sort!([relpath(joinpath(d, f), root) for (d, _, fs) in walkdir(joinpath(root, "src")) for f in fs]; by = p -> replace(p, '\\' => '/'))
    for rel in files
        println(io, replace(rel, '\\' => '/')," ", bytes2hex(_file_digest(joinpath(root, rel))))
    end
    lock = _lock_file(root)
    println(io, "lock ", isempty(lock) ? "none" : bytes2hex(_file_digest(lock)))
    return bytes2hex(sha256(take!(io)))
end

# SHA-256 of the lock in use (the manifest records it as the package set).
packages_hash(root::AbstractString = _root()) = (l = _lock_file(root); isempty(l) ? "none" : bytes2hex(_file_digest(l)))
