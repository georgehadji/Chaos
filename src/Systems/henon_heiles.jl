# Henon-Heiles (catalog F3): the system parts and the energy-surface starts. Spec: workspaces/module-build/stages/01-spec/output/henon-heiles-spec.md
# (Equations 4, 5, 7, 8, Algorithm 2 and 4). Every function is generic in the element type; the section line is q1 = 0 with coordinates (q2, p2).

@inline _hh_u(q2) = q2^2 / 2 - q2^3 / 3            # Equation 7: V(0, q2)

"""
    critical_points(flow::HenonHeilesFlow, T) -> NTuple{4, SVector{2, T}}

The four critical points of V (Equation 4): the minimum `O = (0, 0)` and the saddles `S0 = (0, 1)`, `S1 = (sqrt(3)/2, -1/2)`, `S2 = (-sqrt(3)/2, -1/2)`, correctly rounded in `T`.
"""
function critical_points(::HenonHeilesFlow, ::Type{T}) where {T}
    s = sqrt(T(3)) / 2; h = T(1) / 2
    return (SVector{2, T}(0, 0), SVector{2, T}(0, 1), SVector{2, T}(s, -h), SVector{2, T}(-s, -h))
end

"`escape_energy(flow::HenonHeilesFlow, T)`: the energy `1/6` of the three saddles, above which orbits escape (Equation 6), in `T`."
escape_energy(::HenonHeilesFlow, ::Type{T}) where {T} = T(1) / 6

"""
    d3_elements(T) -> NTuple{6, SMatrix{2, 2, T, 4}}

The symmetry group D3 of V (Equation 5) as orthogonal matrices acting on q and p alike: identity, the rotations by 2 pi/3 and 4 pi/3, the reflection `diag(-1, 1)`,
then the two rotations times the reflection. Entries are correctly rounded.
"""
function d3_elements(::Type{T}) where {T}
    c = -T(1) / 2; s = sqrt(T(3)) / 2
    R1 = SMatrix{2, 2, T, 4}(c, s, -s, c); R2 = SMatrix{2, 2, T, 4}(c, -s, s, c); F = SMatrix{2, 2, T, 4}(-1, 0, 0, 1)
    return (one(SMatrix{2, 2, T, 4}), R1, R2, F, R1 * F, R2 * F)
end

# root of u(q2) = E in [lo, hi] by bisection to the resolution of T; f(lo) and f(hi) have opposite signs
function _hh_bisect(E, lo::T, hi::T) where {T}
    sgn = _hh_u(lo) - E > 0
    for _ in 1:4096
        mid = (lo + hi) / 2
        (mid == lo || mid == hi) && break
        ((_hh_u(mid) - E > 0) == sgn) ? (lo = mid) : (hi = mid)
    end
    return (lo + hi) / 2
end

"""
    section_interval(flow::HenonHeilesFlow, E) -> (a, b)

The roots `a < 0 < b` of `u(q2) = q2^2/2 - q2^3/3 = E` that bound the origin component of the section line `q1 = 0` (Equation 7); `(-1/2, 1)` at `E = 1/6`.
`ArgumentError` unless `0 < E <= 1/6`.
"""
function section_interval(::HenonHeilesFlow, E::T) where {T <: AbstractFloat}
    (isfinite(E) && 0 < E <= T(1) / 6) || throw(ArgumentError("section_interval needs 0 < E <= 1/6, got $E"))
    E == T(1) / 6 && return (-T(1) / 2, T(1))
    return (_hh_bisect(E, -T(1), zero(T)), _hh_bisect(E, zero(T), T(1)))
end
section_interval(flow::HenonHeilesFlow, E::Real) = section_interval(flow, float(E))

"""
    energy_start(flow::HenonHeilesFlow, E, q2, p2) -> (z, status)

The start `z = (0, q2, p1, p2)` on the energy surface `H = E` with `p1 = sqrt(2 (E - p2^2/2 - u(q2))) > 0` (Equation 8), `status = :ok`. `(nothing, :outside)` if
`p2^2/2 + u(q2) >= E` (the boundary `p1 = 0` included, Amendment D117) or `abs(q2) > 1` (the start must lie in the origin component for `E <= 1/6` and inside `r <= 1` above it); `(nothing, :nonfinite)` for a non-finite
argument. Never a NaN state.
"""
function energy_start(::HenonHeilesFlow, E, q2, p2)
    T = float(promote_type(typeof(E), typeof(q2), typeof(p2)))
    (isfinite(E) && isfinite(q2) && isfinite(p2)) || return (nothing, :nonfinite)
    r = 2 * (T(E) - T(p2)^2 / 2 - _hh_u(T(q2)))
    (r > 0 && abs(q2) <= 1) || return (nothing, :outside)                 # D117 (L1): p1 = 0 is not an upward crossing
    return (SVector{4, T}(0, T(q2), sqrt(r), T(p2)), :ok)
end

"""
    random_energy_start(flow::HenonHeilesFlow, E, q2window, p2window, rng; maxtries = 10_000) -> (z, status)

Rejection draw (Algorithm 2): `(q2, p2)` uniform in the two windows (tuples `(lo, hi)`), accepted by [`energy_start`](@ref). `(z, :ok)`, or `(nothing, :no_admissible_point)` after
`maxtries` rejections, or `(nothing, :nonfinite)` for a non-finite `E` or window. The uniform draws are binary64 (`rand(rng)`) converted to the element type.
"""
function random_energy_start(flow::HenonHeilesFlow, E, q2w::Tuple, p2w::Tuple, rng; maxtries::Integer = 10_000)
    T = float(promote_type(typeof(E), typeof.(q2w)..., typeof.(p2w)...))
    (isfinite(E) && all(isfinite, q2w) && all(isfinite, p2w)) || return (nothing, :nonfinite)
    for _ in 1:maxtries
        q2 = T(q2w[1]) + (T(q2w[2]) - T(q2w[1])) * T(rand(rng)); p2 = T(p2w[1]) + (T(p2w[2]) - T(p2w[1])) * T(rand(rng))
        z, st = energy_start(flow, T(E), q2, p2)
        st === :ok && return (z, st)
    end
    return (nothing, :no_admissible_point)
end

"""
    energy_grid(flow::HenonHeilesFlow, E, n) -> Vector{SVector{4, T}}

The portrait starts (Algorithm 4): cell centres of the `n x n` grid over `[a, b] x [-sqrt(2E), sqrt(2E)]` ([`section_interval`](@ref)), kept where inside the section domain,
`p1` from Equation 8. Deterministic; `0 < E <= 1/6`.
"""
function energy_grid(flow::HenonHeilesFlow, E::T, n::Integer) where {T <: AbstractFloat}
    a, b = section_interval(flow, E); pm = sqrt(2E); zs = SVector{4, T}[]
    for i in 1:n, j in 1:n
        z, st = energy_start(flow, E, a + (i - T(1) / 2) * (b - a) / n, -pm + (j - T(1) / 2) * 2pm / n)
        st === :ok && push!(zs, z)
    end
    return zs
end
