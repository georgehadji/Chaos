# Screening grid (E10) and torus distance (E8).

const _TWO_PI = 2 * pi

"""
    torus_delta(s) -> |s - 2 pi round(s / 2 pi)|, the distance of `s` to the nearest multiple of 2 pi (E8).
"""
torus_delta(s) = (tp = 2 * oftype(s, pi); abs(s - tp * round(s / tp)))

"""
    torus_distance(u, v)

Euclidean distance of two points of the torus in the minimum image (E8); `u`, `v` are 2-indexable.
"""
function torus_distance(u, v)
    a = torus_delta(u[1] - v[1])
    b = torus_delta(u[2] - v[2])
    return sqrt(a * a + b * b)
end

"""
    grid_centres(N_g) -> Vector{SVector{2,Float64}}

Centres of the `N_g^2` screening cells (E10); index `id + 1` holds cell `id = j N_g + i`. Each coordinate is the binary64
nearest to `2 pi (2 i + 1) / (2 N_g)`, a single correct rounding from a 256-bit evaluation.
"""
function grid_centres(N_g::Integer)
    N_g >= 1 || throw(ArgumentError("grid_centres: N_g must be positive"))
    col = setprecision(BigFloat, 256) do
        [Float64(2 * BigFloat(pi) * (2i + 1) / (2N_g)) for i in 0:(N_g - 1)]
    end
    return [SVector(col[i], col[j]) for j in 1:N_g for i in 1:N_g]
end
