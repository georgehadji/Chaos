# Escape of a flow orbit above the escape energy (henon-heiles spec Equation 14, Algorithm 5): the first passage of the circle r = R_esc, located by a sub-step, and the exit
# sector by the polar angle. Uses the root solver of section.jl.

@inline _radius(w, c) = sqrt(w[1]^2 + w[2]^2)
@inline _rad_gap(w, R) = _radius(w, R) - R
@inline _rad_rate(w, c) = (w[1] * w[3] + w[2] * w[4]) / _radius(w, c)

"""
    escape_sector(flow::HenonHeilesFlow, q) -> Int

The exit `0`, `1` or `2` whose centre angle (`pi/2`, `7 pi/6`, `11 pi/6`) is nearest to the polar angle of `q` (Equation 14). The rotation by `2 pi/3` maps `k` to `k + 1` (mod 3), the
reflection `q1 -> -q1` swaps 1 and 2 (also at `q1 = +-0`: a 1-2 tie goes to 1 when the sign bit of `q1` is set, Amendment D117). Allocation free.
"""
function escape_sector(::HenonHeilesFlow, q)
    T = float(eltype(q)); h = T(1) / 2; s = sqrt(T(3)) / 2
    d0 = q[2]; d1 = -s * q[1] - h * q[2]; d2 = s * q[1] - h * q[2]      # q . u_k for the unit exit vectors u_k: the nearest centre angle is the largest dot product
    return (d0 >= d1 && d0 >= d2) ? 0 : ((d1 > d2 || (d1 == d2 && signbit(q[1]))) ? 1 : 2)    # D117 (L2): a 1-2 tie (q1 = +-0) by the sign bit, so the reflection swaps it
end

"""
    locate_escape(m, flow, z, tau, R) -> (zc, s)

For a state `z` with `r <= R` whose step of size `tau` has `r > R`: the root `s` in `(0, tau]` of `r(Phi_s z) = R` by the method's own step, and `zc = advance(m, flow, z, s)`
with `abs(r(zc) - R) <= 4 eps R`. Allocation free for `isbits` states.
"""
locate_escape(m, flow, z::SVector{4, T}, tau, R) where {T} = _solve_substep(m, flow, z, tau, _rad_gap, _rad_rate, T(R), 4 * effective_eps(T) * T(R))

"""
    escape_orbit(m, flow, z0, tau, n; R_esc = 10) -> (t, sector, status, z)

Integrates `n` steps and declares escape at the first step with `r > R_esc`, located by [`locate_escape`](@ref): `status = :escaped`, `t` the passage time, `z` the state at
`r = R_esc` and `sector` of [`escape_sector`](@ref). `:bound`: no escape (`t` NaN, `sector` -1, `z` the final state). `:nonfinite`: a non-finite state or overflow (`z` the last finite
state). Never throws.
"""
function escape_orbit(m, flow, z0::SVector{4, T}, tau, n::Integer; R_esc = 10) where {T}
    R = T(R_esc); z = z0
    all(isfinite, z) || return (t = T(NaN), sector = -1, status = :nonfinite, z = z)
    _radius(z, R) > R && return (t = zero(T), sector = escape_sector(flow, SVector(z[1], z[2])), status = :escaped, z = z)
    for k in 0:(n - 1)
        zn = advance(m, flow, z, tau)
        all(isfinite, zn) || return (t = T(NaN), sector = -1, status = :nonfinite, z = z)
        if _radius(zn, R) > R
            zc, s = locate_escape(m, flow, z, tau, R)
            return (t = T(k) * tau + s, sector = escape_sector(flow, SVector(zc[1], zc[2])), status = :escaped, z = zc)
        end
        z = zn
    end
    return (t = T(NaN), sector = -1, status = :bound, z = z)
end
