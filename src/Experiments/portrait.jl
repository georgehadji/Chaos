# Phase portrait content (spec Algorithm 12): data for the Figures layer, derived from the panel file only.

const _SYMBOL = Dict("hyperbolic" => :cross, "elliptic" => :circle, "parabolic" => :diamond)

function _start_marker(o)
    w = wrap(SVector(Float64(o["x0"]), Float64(o["p0"])))
    return (id = o["id"], class = o["class"], x = w[1], p = w[2])
end

function _iterates(o, n_plot)
    sys = StandardMap(Float64(o["K"]))
    z = SVector(Float64(o["x0"]), Float64(o["p0"]))
    xy = Vector{SVector{2, Float64}}(undef, n_plot)
    for t in 1:n_plot
        xy[t] = wrap(z)
        z = forward(sys, z)
    end
    return (id = o["id"], class = o["class"], xy = xy)
end

function _periodic_marker(q, pt)
    w = wrap(SVector(Float64(parse(BigFloat, pt[1])), Float64(parse(BigFloat, pt[2]))))
    return (id = q["id"], symbol = _SYMBOL[q["class"]], x = w[1], p = w[2])
end

function _notes(cells)
    return ["class $(c["class"]): $(c["status"]) ($(c["reason"]))" for c in cells if c["status"] in ("empty", "short")]
end

"""
    portrait_content(panel_or_path; n_plot)

Content of the phase portraits, one NamedTuple per K (spec Algorithm 12): `starts` (id, class, wrapped start), `iterates` (id, class,
the first `n_plot` wrapped iterates, the first being the start), `periodic` (one marker with a symbol per periodic orbit point:
`:cross` hyperbolic, `:circle` elliptic, `:diamond` parabolic) and `notes` on empty and short cells. No screening data.
"""
function portrait_content(panel::AbstractDict; n_plot::Integer)
    man, orbs = panel["manifest"], panel["orbit"]
    return setprecision(BigFloat, 256) do
        map(enumerate(man["K_binary64"])) do (i, K)
            here = [o for o in orbs if get(o, "k_index", 0) == i]
            per = [q for q in man["periodic"] if q["k_index"] == i]
            cells = [c for c in man["cells"] if c["k_index"] == i]
            (k_index = i, K = K, starts = map(_start_marker, here), iterates = [_iterates(o, n_plot) for o in here],
             periodic = [_periodic_marker(q, pt) for q in per for pt in q["points"]], notes = _notes(cells))
        end
    end
end

portrait_content(path::AbstractString; n_plot::Integer) = portrait_content(TOML.parsefile(path); n_plot = n_plot)
