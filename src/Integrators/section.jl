# Surface of section of a flow orbit (henon-heiles spec Equations 10 to 12, Algorithm 3): the plane q1 = 0 crossed upward (p1 > 0), each crossing located by a sub-step of the
# method's own one-step map, with the optional 2 x 2 Jacobian of the return map in (q2, p2). A root solver shared with the escape tools (escape.jl).

# Root s in (0, tau] of g(Phi_s z, c) = 0 with g(z) < 0 <= g(Phi_tau z): Newton steps with derivative dg along the flow, bisection when a step leaves the bracket.
# Returns (Phi_s z, s) with |g| <= tol, or the upper end of a collapsed bracket. Allocation free for isbits states.
function _solve_substep(m, flow, z, tau, g::G, dg::D, c, tol) where {G, D}
    lo = zero(tau); hi = tau; s = tau / 2
    w = advance(m, flow, z, s)
    for _ in 1:200
        gv = g(w, c)
        abs(gv) <= tol && return (w, s)
        gv < 0 ? (lo = s) : (hi = s)
        sn = s - gv / dg(w, c)
        mid = (lo + hi) / 2
        (mid == lo || mid == hi) && return (advance(m, flow, z, hi), hi)
        s = (lo < sn < hi) ? sn : mid
        w = advance(m, flow, z, s)
    end
    return (w, s)
end

@inline _q1(w, _) = w[1]
@inline _dq1(w, _) = w[3]

"""
    locate_crossing(m, flow, z, tau) -> (zc, s)

For a state `z` with `q1 < 0` whose step of size `tau` has `q1 >= 0`: the root `s` in `(0, tau]` of `q1(Phi_s z) = 0` found with the method's own step (Newton with a bisection
safeguard), and `zc = advance(m, flow, z, s)` with `abs(zc[1]) <= 4 eps`. Allocation free for `isbits` states.
"""
locate_crossing(m, flow, z::SVector{4, T}, tau) where {T} = _solve_substep(m, flow, z, tau, _q1, _dq1, zero(T), 4 * effective_eps(T))

# 2 x 2 return-map Jacobian in (q2, p2) of the leg from the section point zp to zc (Equation 12): DP = (I - f e1' / p1) Phi restricted to the energy surface (dp1 from dH = 0)
function _leg_jacobian(flow, Phi::SMatrix{4, 4, T}, zp, zc) where {T}
    g2 = force(flow, SVector(zp[1], zp[2]))[2]; f = vector_field(flow, zc)
    w1 = Phi * SVector{4, T}(0, 1, -g2 / zp[3], 0); w2 = Phi * SVector{4, T}(0, 0, -zp[4] / zp[3], 1)
    w1 = w1 - f * (w1[1] / f[1]); w2 = w2 - f * (w2[1] / f[1])
    return SMatrix{2, 2, T, 4}(w1[2], w1[4], w2[2], w2[4])
end

# the next state, and the tangent frame advanced across the same step when a Jacobian is wanted
@inline _frame_step(m, flow, z, Phi, tau, ::Val{false}) = Phi
@inline _frame_step(m, flow, z, Phi, tau, ::Val{true}) = last(tangent_advance(m, flow, z, Phi, tau))

"""
    section_map(m, flow, z0, tau, n; jacobian = false, max_crossings = typemax(Int))

The upward crossings of the plane `q1 = 0` (`p1 > 0`) of the orbit of `z0` over `n` steps of `tau`. Returns a `NamedTuple` `(t, z, J, grazing, status)`: crossing times, full crossing
states (`q1` about 0), the Jacobians of the legs in `(q2, p2)` (only with `jacobian`; `J[i]` is the leg from section point `i - 1` to `i`, `J[1]` the identity when the start is on
the section and a placeholder identity otherwise), `grazing` the count of bracketed crossings with `abs(p1) < 8 eps` (flagged, not recorded), and `status` in `:ok`, `:nonfinite`,
`:buffer_full`. A start with `q1 = 0`, `p1 > 0` is the first point at `t = 0`. Never throws on non-finite values.
"""
function section_map(m, flow, z0::SVector{4, T}, tau, n::Integer; jacobian::Bool = false, max_crossings::Integer = typemax(Int)) where {T}
    ts = T[]; zs = SVector{4, T}[]; Js = SMatrix{2, 2, T, 4}[]
    all(isfinite, z0) || return (t = ts, z = zs, J = Js, grazing = 0, status = :nonfinite)
    if z0[1] == 0 && z0[3] > 0
        push!(ts, zero(T)); push!(zs, z0); jacobian && push!(Js, one(SMatrix{2, 2, T, 4}))
    end
    z = z0; Phi = one(SMatrix{4, 4, T}); graz = 0; status = :ok; jv = Val(jacobian)
    length(ts) >= max_crossings && return (t = ts, z = zs, J = Js, grazing = 0, status = :buffer_full)
    for k in 0:(n - 1)
        zn = advance(m, flow, z, tau)
        all(isfinite, zn) || (status = :nonfinite; break)
        Phin = _frame_step(m, flow, z, Phi, tau, jv)
        if z[1] < 0 <= zn[1] && zn[3] > 0
            zc, s = locate_crossing(m, flow, z, tau)
            if abs(zc[3]) < 8 * eps(T)
                graz += 1
            else
                jacobian && push!(Js, isempty(zs) ? one(SMatrix{2, 2, T, 4}) : _leg_jacobian(flow, _frame_step(m, flow, z, Phi, s, jv), zs[end], zc))
                push!(ts, T(k) * tau + s); push!(zs, zc)
                jacobian && (Phin = _frame_step(m, flow, zc, one(SMatrix{4, 4, T}), tau - s, jv))
                length(ts) >= max_crossings && (status = :buffer_full; break)
            end
        end
        z = zn; Phi = Phin
    end
    return (t = ts, z = zs, J = Js, grazing = graz, status = status)
end
