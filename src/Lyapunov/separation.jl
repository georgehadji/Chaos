# Separation hooks of the two-trajectory estimator (spec lyap-two-trajectory T2, T3, Algorithm 1 and 3).

"""
    separation(sys, ref, shadow) -> SVector{d}

Separation vector `shadow - ref` on the `d` tangent coordinates. The default is the plain difference (lifted states);
`WrappedStandardMap` reduces each coordinate to its minimal image in `(-pi, pi]` (spec T3).
"""
separation(sys, ref, shadow) = shadow - ref

"""
    displaced(sys, ref, delta) -> state

Shadow state at `ref + delta`. The default is `ref + delta`; `WrappedStandardMap` wraps both coordinates and adds the
reduction counts to those of `ref`, so the lifted separation is exactly the minimal-image one (spec T2, T3).
"""
displaced(sys, ref, delta) = ref + delta

# minimal image into (-pi, pi]
@inline function _minimal_image(u, tp)
    h = tp / 2
    return u > h ? u - tp : (u <= -h ? u + tp : u)
end

function separation(::WrappedStandardMap, ref::SVector{4}, shadow::SVector{4})
    tp = 2 * Systems._pi(ref[1])
    return SVector(_minimal_image(shadow[1] - ref[1], tp), _minimal_image(shadow[2] - ref[2], tp))
end

function displaced(::WrappedStandardMap, ref::SVector{4}, delta::SVector{2})
    tp = 2 * Systems._pi(ref[1])
    x, a = Systems._reduce(ref[1] + delta[1], tp)
    p, b = Systems._reduce(ref[2] + delta[2], tp)
    return SVector(x, p, ref[3] + a, ref[4] + b)
end

state_mode(::Any) = "lifted"
state_mode(::WrappedStandardMap) = "wrapped"
