# Study p2-indicators (Amendment D89): the bootstrap of E34. Streams are the streams of E12 (`stream_seed`, `bootstrap_rng`) keyed by the label
# `p2|<item>|<set>|<set rule>|<verdict rule>|<variant>`; the draws of a replicate depend on the label and the replicate only (I20).

"The stream label of E34."
p2_label(item, set, set_rule, rule, variant) = "p2|$item|$set|$set_rule|$rule|$variant"

"The `m x B` matrix of indices: column r holds the m draws (uniform with replacement) of replicate r, in the order of the stream of `label` (E34)."
function p2_resample(label::AbstractString, m::Integer, B::Integer; seed::Integer)
    rng = bootstrap_rng(label, seed)
    return Int[rand(rng, 1:m) for _ in 1:m, _ in 1:B]
end

# The scalar statistics that carry an interval, by hypothesis.
function _scalars(h::Symbol, s)
    h === :H1 && return Dict{Symbol, Any}(:U => s.U)
    h === :H2 && return Dict{Symbol, Any}(:Q => s.Q)
    h === :H3 && return Dict{Symbol, Any}(Symbol(:R_, X) => v for (X, v) in s.R)
    h === :H4 && return Dict{Symbol, Any}(:dphi => s.dphi)
    h === :H5 && return Dict{Symbol, Any}(:F_med => s.F_med)
    throw(ArgumentError("hypothesis $h has no bootstrap interval"))
end

"""
    p2_bootstrap(h, sets, labels, p; I = I_7, seed = p.bootstrap_seed) -> (interval, dropped, defined)

Percentile bootstrap (E34) of the scalar statistics of hypothesis `h` (`:H1` to `:H5`): `p.bootstrap_resamples` resamples of the orbits with replacement, one
stream per set (`labels` is a text, or two texts for H3: W then S, resampled independently), all statistics of a label sharing the draws of the replicate.
A replicate on which a statistic is undefined is dropped and counted; above 5 percent dropped the interval is `nothing`. `interval` maps `:U`, `:Q`, `:R_<X>`,
`:dphi`, `:F_med` to `(lo, hi)`.
"""
function p2_bootstrap(h::Symbol, sets::Tuple, labels, p; I = P2_I7, seed::Integer = p.bootstrap_seed)
    B = p.bootstrap_resamples
    ls = labels isa AbstractString ? (labels,) : Tuple(labels)
    ms = Tuple(length(s) for s in sets)
    any(==(0), ms) && return (interval = nothing, dropped = B, defined = false)
    rngs = Tuple(bootstrap_rng(l, seed) for l in ls)
    kp = _kp(p)
    reps = Dict{Symbol, Vector{Float64}}()
    dropped = 0
    for _ in 1:B
        rs = Tuple(sets[j][Int[rand(rngs[j], 1:ms[j]) for _ in 1:ms[j]]] for j in eachindex(sets))
        sc = _scalars(h, p2_stat(h, rs, p; I = I, kp = kp))
        if any(isnothing, values(sc))
            dropped += 1
            continue
        end
        foreach(((k, v),) -> push!(get!(reps, k, Float64[]), Float64(v)), sc)
    end
    a = (1 - p.bootstrap_level) / 2
    iv = Dict{Symbol, NamedTuple}(k => (lo = _p2_percentile(sort(v), a), hi = _p2_percentile(sort(v), 1 - a)) for (k, v) in reps)
    ok = dropped <= 0.05 * B
    return (interval = ok ? iv : nothing, dropped = dropped, defined = ok)
end
