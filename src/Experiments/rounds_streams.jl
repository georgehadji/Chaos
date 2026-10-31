# Streams and draws of Amendment D94: screening points (E32), round streams (E33), pools and the used set (E34, E35), subset D (E37)
# and subset D_B (E38). Pure functions of (seed, text); nothing here depends on threads or devices.

_seed64(text) = foldl((acc, b) -> (acc << 8) | UInt64(b), sha256(text)[1:8]; init = UInt64(0))

"Seed of the screening stream of K index `i` (E32): first 8 bytes, big-endian, of SHA-256 over `chaoslab-screen-1|<seed>|<i>`."
screen_stream_seed(seed::Integer, i::Integer) = _seed64("chaoslab-screen-1|$(seed)|$(i)")

"Seed of the stream of round cell (`r`, `i`, `class`) (E33): the text of E15 with `|<r>` appended."
round_stream_seed(seed::Integer, i::Integer, class, r::Integer) = _seed64("chaoslab-panel-1|$(seed)|$(i)|$(class)|$(r)")

"Seed of the stream of the subset D of round cell (`r`, `i`, `class`) (E37)."
d_stream_seed(seed::Integer, r::Integer, i::Integer, class) = _seed64("chaoslab-d-1|$(seed)|$(r)|$(i)|$(class)")

"Seed of the stream of stratum `i` of subset D_B for `R` rounds run (E38)."
db_stream_seed(seed::Integer, R::Integer, i::Integer) = _seed64("chaoslab-db-1|$(seed)|$(R)|$(i)")

"""
    random_screen_points(seed, i, M) -> Vector{SVector{2,Float64}}

The first `M` screening points of K index `i` (E32): `(x, p) = (RN(2 pi u), RN(2 pi u'))` for the uniforms of the stream in order, `RN` the
binary64 nearest to the 256-bit product. Sequential on one stream, so independent of threads and chunking.
"""
function random_screen_points(seed::Integer, i::Integer, M::Integer)
    rng = StableRNG(screen_stream_seed(seed, i))
    return setprecision(BigFloat, 256) do
        tp = 2 * BigFloat(pi)
        [SVector(Float64(tp * BigFloat(rand(rng))), Float64(tp * BigFloat(rand(rng)))) for _ in 1:M]
    end
end

"`m = floor(f n + 1/2)` orbits of a round cell with `n` orbits drawn and share `f` (E37), in exact integer arithmetic (I36)."
d_size(f::Rational, n::Integer) = Int(div(2 * big(numerator(f)) * n + big(denominator(f)), 2 * big(denominator(f))))

"Sorted draw indices (of `1:n`) of subset D in round cell (`r`, `i`, `class`) with share `f` (E37)."
d_indices(seed::Integer, r::Integer, i::Integer, class, n::Integer, f::Rational) =
    sort!(draw_candidates(StableRNG(d_stream_seed(seed, r, i, class)), 1:n, d_size(f, n)))

"""
    db_allocation(n, b) -> Vector{Int}

Largest-remainder allocation of `b` to the strata of sizes `n` (E38): `floor(b n_i / sum(n))` plus one for each of the largest remainders,
ties to the smaller index; `n` itself when `sum(n) <= b`. Exact integer arithmetic.
"""
function db_allocation(n::AbstractVector{<:Integer}, b::Integer)
    T = sum(n)
    T <= b && return collect(Int, n)
    q = [Int(div(big(b) * ni, T)) for ni in n]
    rem = [mod(big(b) * ni, T) for ni in n]
    out = copy(q)
    for i in first(sortperm(eachindex(n); by = i -> (-rem[i], i)), b - sum(q))
        out[i] += 1
    end
    return out
end

"""
    db_subset(seed, R, DR, b) -> Vector{Int}

Subset D_B (E38) as sorted ids: `DR` holds `(id, k_index)` of the orbits of D in the rounds run, `b` is `db_size`, `R` the number of rounds run.
"""
function db_subset(seed::Integer, R::Integer, DR, b::Integer)
    isempty(DR) && return Int[]
    nK = maximum(x.k_index for x in DR)
    al = db_allocation([count(x -> x.k_index == i, DR) for i in 1:nK], b)
    out = Int[]
    for i in 1:nK
        al[i] == 0 && continue
        ids = sort!([x.id for x in DR if x.k_index == i])
        append!(out, draw_candidates(StableRNG(db_stream_seed(seed, R, i)), ids, al[i]))
    end
    return sort!(out)
end

"Status and reason of a round cell with `n` orbits drawn of size `q` (E35)."
function _round_status(cls, n, q)
    q == 0 && return "not drawn", ""
    n == q && return "populated", ""
    n > 0 && return "short", "fewer qualifying points than the quota"
    return "empty", cls == "R" ? "no elliptic periodic point" : "no qualifying point"
end

const _POOL_BIT = Dict("T" => _BIT_T, "C" => _BIT_C, "S" => _BIT_S)

"""
    draw_pool_rounds(seed, i, bits; quota, r_max, active = 1:r_max)

The draw of the pool classes S, C, T of K index `i` by round (E34, E35) from the class bit vector `bits` (one `UInt8` per point id; T `0x01`,
C `0x02`, S `0x04`). `quota` maps a class to its size, `active` lists the rounds in which the cell is drawn. A point drawn leaves every
pool of the K (the used set). Returns one `(round, class, pool, q, drawn, status, reason)` per round and class, in the order of the draw;
`pool` is the candidate count at the moment of the draw.
"""
function draw_pool_rounds(seed::Integer, i::Integer, bits::AbstractVector{UInt8}; quota, r_max::Integer, active = 1:r_max)
    used = falses(length(bits))
    cells = NamedTuple[]
    for r in 1:r_max, c in ("S", "C", "T")
        q = r in active ? get(quota, c, 0) : 0
        bit = _POOL_BIT[c]
        cand = [j for j in eachindex(bits) if bits[j] & bit != 0 && !used[j]]
        drawn = q > 0 ? draw_candidates(StableRNG(round_stream_seed(seed, i, c, r)), cand, q) : Int[]
        foreach(j -> used[j] = true, drawn)
        status, reason = _round_status(c, length(drawn), q)
        push!(cells, (round = r, class = c, pool = length(cand), q = q, drawn = drawn, status = status, reason = reason))
    end
    return cells
end
