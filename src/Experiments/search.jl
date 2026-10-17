# The periodic set of one K: exact fixed points plus a Newton search from every screening cell centre for each period (spec
# Algorithm 3), the rational lines at K = 0 (E7).

"""
The periodic set at one K: `orbits` (sorted by period and representative), `lines` the `(a, b)` pairs of the K = 0 lines
`p = 2 pi a / b` (empty for K > 0), `stats` the Newton counts, `min_separation` the smallest distance between points of
distinct orbits, `ill` the ill-conditioned rejections per period (count and smallest `|2 - t|`). `points` and `elliptic` hold the Float64 points of all orbits and of the elliptic and parabolic ones.
"""
struct PeriodicSet
    orbits::Vector{PeriodicOrbit}
    lines::Vector{Tuple{Int, Int}}
    stats::NamedTuple{(:tried, :accepted, :singular, :diverged, :cap, :ill_conditioned), NTuple{6, Int}}
    ill::Vector{NamedTuple{(:q, :count, :min_den), Tuple{Int, Int, Float64}}}
    min_separation::Float64
    points::Vector{SVector{2, Float64}}
    elliptic::Vector{SVector{2, Float64}}
end

_f64_points(orbits, keep) = [SVector(Float64(pt[1]), Float64(pt[2])) for o in orbits if keep(o) for pt in o.points]

function PeriodicSet(orbits, lines, stats, ill, sep)
    return PeriodicSet(orbits, lines, stats, ill, sep, _f64_points(orbits, _ -> true),
                       _f64_points(orbits, o -> o.class !== :hyperbolic))
end

const _NO_ILL = NamedTuple{(:q, :count, :min_den), Tuple{Int, Int, Float64}}[]
const _NO_STATS = (tried = 0, accepted = 0, singular = 0, diverged = 0, cap = 0, ill_conditioned = 0)

# Newton solves of period q from every centre, in parallel; results by cell id.
function _search_period(Kb, centres, q::Integer, rho, nmax::Integer, delta_dup, nt::Integer)
    res = Vector{Tuple{Int, BigFloat, BigFloat, Float64}}(undef, length(centres))
    _par(() -> nothing, length(centres), nt) do _, i
        res[i] = _newton_seed(Kb, BigFloat(centres[i][1]), BigFloat(centres[i][2]), q, rho, nmax, delta_dup)
    end
    return res
end

function _count_codes(res, st)
    c = zeros(Int, 5)
    for r in res
        c[r[1] + 1] += 1
    end
    return (tried = st.tried + length(res), accepted = st.accepted + c[1], singular = st.singular + c[2],
            diverged = st.diverged + c[3], cap = st.cap + c[4], ill_conditioned = st.ill_conditioned + c[5])
end

"""
    periodic_set(K; q_max, N_g, rho, newton_max, delta_dup, bits = 256, ntasks = nothing)

Periodic orbits of period up to `q_max` of the standard map at `K`, found by Newton's method at `bits` bits from the `N_g^2`
screening cell centres, plus the exact fixed points (E5). At `K = 0` there is no search: `lines` lists the reduced fractions
`a / b`, `b <= q_max`, of E7. The result does not depend on `ntasks`. A Newton solution is accepted only if its position error bound is at most `delta_dup / 10`
(Amendment D51); the others are counted as ill-conditioned. The separation rule of E9 is applied by the caller on `min_separation`.
"""
function periodic_set(K::Real; q_max::Integer, N_g::Integer, rho, newton_max::Integer, delta_dup, bits::Integer = 256,
                      ntasks = nothing)
    if K == 0
        lines = [(a, b) for b in 1:q_max for a in 0:(b - 1) if gcd(a, b) == 1]
        return PeriodicSet(PeriodicOrbit[], lines, _NO_STATS, _NO_ILL, Inf)
    end
    centres, nt = grid_centres(N_g), _ntasks(ntasks)
    return setprecision(BigFloat, bits) do
        Kb = BigFloat(K)
        zero_b = zero(BigFloat)
        sols = [(q = 1, z = (zero_b, zero_b), id = -1), (q = 1, z = (BigFloat(pi), zero_b), id = -1)]
        stats, ill = _NO_STATS, similar(_NO_ILL, 0)
        for q in 1:q_max
            res = _search_period(Kb, centres, q, rho, newton_max, delta_dup, nt)
            stats = _count_codes(res, stats)
            bad = [r[4] for r in res if r[1] == 4]
            isempty(bad) || push!(ill, (q = q, count = length(bad), min_den = minimum(bad)))
            append!(sols, (q = q, z = (r[2], r[3]), id = i - 1) for (i, r) in enumerate(res) if r[1] == 0)
        end
        orbits = collect_orbits(K, sols; rho = rho, delta_dup = delta_dup, bits = bits, check = false)
        PeriodicSet(orbits, Tuple{Int, Int}[], stats, ill, _min_separation(orbits))
    end
end

"""
    periodic_distance(z, pset; elliptic_only = false)

Torus distance (E8) from `z` to the periodic set: to the lines `p = 2 pi a / b` at K = 0, else to the nearest orbit point. With
`elliptic_only` parabolic points count as elliptic and hyperbolic ones are left out (E11). `Inf` for an empty set.
"""
function periodic_distance(z, ps::PeriodicSet; elliptic_only::Bool = false)
    if !isempty(ps.lines)
        return minimum(torus_delta(z[2] - _TWO_PI * a / b) for (a, b) in ps.lines)
    end
    best = Inf
    for pt in (elliptic_only ? ps.elliptic : ps.points)
        best = min(best, torus_distance(z, pt))
    end
    return best
end
