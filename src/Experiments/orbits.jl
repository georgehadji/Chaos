# Orbit records, primitive period, canonical representative and deduplication (E9, E20).

"""
A periodic orbit of primitive period `q`: windings `(m, j)` of the stored representative `points[1]` (E2, E3), the points wrapped
to [0, 2 pi) at working precision, the trace of the monodromy matrix, its class (E4), the exponent per map step, the return
residual of the representative and the grid cell id of the solution kept (-1 for an exact fixed point).
"""
struct PeriodicOrbit
    q::Int
    m::Int
    j::Int
    points::Vector{SVector{2, BigFloat}}
    trace::BigFloat
    class::Symbol
    exponent::BigFloat
    residual::Float64
    grid_id::Int
end

# x mod 2 pi in [0, 2 pi) at the ambient precision.
function _wrap_big(x)
    tp = 2 * oftype(x, pi)
    y = x - tp * floor(x / tp)
    y < 0 && (y += tp)
    return y >= tp ? zero(y) : y
end

_step_wrapped(K, x, p) = (q = p + K * sin(x); (_wrap_big(x + q), _wrap_big(q)))

# Sort key of the canonical representative (E20): the wrapped point rounded to binary64; a coordinate that rounds up to 2 pi
# counts as 0, the torus point it is.
function _key(pt)
    f(v) = (y = Float64(v); y >= _TWO_PI ? 0.0 : y)
    return (f(pt[1]), f(pt[2]))
end

# The q points of the orbit of (x, p), wrapped, the canonical representative first (the others in iteration order).
function _orbit_points(K, x, p, q::Integer)
    pts = Vector{SVector{2, BigFloat}}(undef, q)
    x, p = _wrap_big(x), _wrap_big(p)
    for t in 1:q
        pts[t] = SVector(x, p)
        x, p = _step_wrapped(K, x, p)
    end
    return circshift(pts, 1 - argmin(map(_key, pts)))
end

# Smallest divisor of `q` at which the point returns within `rho` (its primitive period, E9); 0 if none does.
function _primitive_period(K, x, p, q::Integer, rho)
    for d in 1:q
        q % d == 0 && accept_residual(_residual(K, x, p, d), rho) && return d
    end
    return 0
end

function _make_orbit(K, pts, q::Integer, rho, id::Integer)
    x, p = pts[1][1], pts[1][2]
    fx, fp, m, j, a, _, _, d = _return_map(K, x, p, q)
    tr = a + d
    cls = classify_trace(tr, rho)
    ex = cls === :hyperbolic ? orbit_exponent(tr, q) : zero(tr)
    return PeriodicOrbit(q, m, j, pts, tr, cls, ex, Float64(max(abs(fx), abs(fp))), id)
end

# Smallest torus distance between points of different orbits (Float64, sweep over x with the torus copy at x + 2 pi); Inf for
# fewer than two orbits.
function _min_separation(orbits)
    entries = Tuple{Float64, Float64, Int}[]
    for (k, o) in enumerate(orbits), pt in o.points
        x, p = Float64(pt[1]), Float64(pt[2])
        push!(entries, (x, p, k), (x + _TWO_PI, p, k))
    end
    sort!(entries; by = first)
    best = Inf
    for i in eachindex(entries)
        xi, pi_, ki = entries[i]
        for jj in (i + 1):length(entries)
            xj, pj, kj = entries[jj]
            xj - xi >= best && break
            kj == ki && continue
            best = min(best, torus_distance((xi, pi_), (xj, pj)))
        end
    end
    return best
end

function _check_separation(sep, delta_dup)
    sep >= 10 * delta_dup || throw(ArgumentError(
        "delta_dup = $delta_dup: distinct orbit points are only $sep apart, below 10 delta_dup (ambiguous orbit identity)"))
    return sep
end

"""
    collect_orbits(K, sols; rho, delta_dup, bits = 256, check = true)

Deduplicate solutions `(q, z, id)` of the return equation at `K` into periodic orbits (E9): each solution is reduced to its primitive
period, solutions are taken in ascending `(period, id)` order (id -1 is an exact fixed point and comes first), and one lying
within `delta_dup` of a point of an orbit already kept is dropped. The orbits are returned sorted by `(q, representative)`.
With `check`, an `ArgumentError` naming `delta_dup` when distinct orbit points are closer than `10 delta_dup`.
"""
function collect_orbits(K, sols; rho, delta_dup, bits::Integer = 256, check::Bool = true)
    return setprecision(BigFloat, bits) do
        Kb = BigFloat(K)
        red = NamedTuple{(:q, :id, :x, :p), Tuple{Int, Int, BigFloat, BigFloat}}[]
        for s in sols
            x, p = BigFloat(s.z[1]), BigFloat(s.z[2])
            q = _primitive_period(Kb, x, p, s.q, rho)
            q > 0 && push!(red, (q = q, id = s.id, x = _wrap_big(x), p = _wrap_big(p)))
        end
        sort!(red; by = r -> (r.q, r.id))
        orbits = PeriodicOrbit[]
        known = Dict{Int, Vector{SVector{2, Float64}}}()
        for r in red
            list = get!(known, r.q, SVector{2, Float64}[])
            any(k -> torus_distance(k, (Float64(r.x), Float64(r.p))) <= delta_dup, list) && continue
            o = _make_orbit(Kb, _orbit_points(Kb, r.x, r.p, r.q), r.q, rho, r.id)
            push!(orbits, o)
            append!(list, (SVector(Float64(pt[1]), Float64(pt[2])) for pt in o.points))
        end
        sort!(orbits; by = o -> (o.q, _key(o.points[1])))
        check && _check_separation(_min_separation(orbits), delta_dup)
        orbits
    end
end
