# The numbers behind each figure F1 to F7 (E16, E19): tables with a header and rows, rendered into the tables file and drawn by the extension.

_pmap(f, xs, threads::Integer) = begin
    out = Vector{Any}(undef, length(xs))
    if threads > 1 && Threads.nthreads() > 1
        Threads.@threads for i in eachindex(xs)
            out[i] = f(xs[i])
        end
    else
        for i in eachindex(xs)
            out[i] = f(xs[i])
        end
    end
    return out
end

_quantile(v, q) = _quantile_sorted(sort(Float64.(v)), q)

# Groups of non-control panel orbits by (class, k_index), in sorted order, as indices into the panel list.
function _cell_groups(orbits)
    g = Dict{Tuple{String, Int}, Vector{Int}}()
    for (i, o) in enumerate(orbits)
        o["control"] || push!(get!(g, (String(o["class"]), Int(o["k_index"])), Int[]), i)
    end
    return sort!(collect(g); by = kv -> (kv[1][2], kv[1][1]))
end

# F1: all methods and the Delta_V band for the first drawn orbit of each cell.
function _f1(orbits, recs, idx, p)
    rows = Vector{Any}[]
    for ((cls, k), ix) in _cell_groups(orbits)
        i = ix[argmin([orbits[j]["draw_index"] for j in ix])]
        r = recs[i]
        l2 = _primary(idx, r.id, "L2", p.N_1, p).lam
        for (s, n) in enumerate(r.steps)
            v = [r.series[m][s] for m in V_METHODS]
            push!(rows, Any[cls, k, r.id, n, v[1], l2[s], v[2], v[3], v[4], v[5], v[6], v[7], v[8], maximum(v) - minimum(v)])
        end
    end
    return (columns = ["class", "k_index", "orbit", "n", "L1", "L2", "L3", "L4", "L4b", "L5", "L5b", "L6", "L9", "DeltaV"], rows = rows)
end

# F2: R(n) median and bootstrap band per cell.
function _f2(orbits, recs, sets, p, seed, threads)
    dropped = Set(d.id for d in sets.dropped)
    keep(i) = !(Int(orbits[i]["id"]) in dropped) && recs[i].lam_ref > 0 && all(isfinite, recs[i].R)
    groups = [(key, [i for i in ix if keep(i)]) for (key, ix) in _cell_groups(orbits)]
    wset = Set(sets.W)
    push!(groups, (("W", 0), [i for i in eachindex(orbits) if Int(orbits[i]["id"]) in wset && keep(i)]))      # pooled over the W orbits of H2
    groups = [(k, ix) for (k, ix) in groups if !isempty(ix)]
    bands = _pmap(g -> bootstrap_band(permutedims(reduce(hcat, [recs[i].R for i in g[2]])), "F2/$(g[1][1])/$(g[1][2])"; seed = seed,
                                       resamples = p.bootstrap_resamples, level = p.bootstrap_level), groups, threads)
    rows = Vector{Any}[]
    for ((key, ix), b) in zip(groups, bands), (s, n) in enumerate(recs[ix[1]].steps)
        push!(rows, Any[key[1], key[2], n, length(ix), b.median[s], b.lo[s], b.hi[s]])
    end
    return (columns = ["class", "k_index", "n", "orbits", "median", "lo", "hi"], rows = rows)
end

# F3: convergence time of every kept orbit by cell (the empirical distribution functions are drawn from these).
function _f3(orbits, recs, sets)
    dropped = Set(d.id for d in sets.dropped)
    rows = Vector{Any}[]
    for ((cls, k), ix) in _cell_groups(orbits), i in sort(ix; by = j -> recs[j].conv.n)
        recs[i].id in dropped || push!(rows, Any[cls, k, recs[i].id, recs[i].conv.n, recs[i].conv.defined])
    end
    return (columns = ["class", "k_index", "orbit", "n_conv", "defined"], rows = rows)
end

# F4: share of orbits forced to unresolved by each test, among the orbits the test could be applied to.
function _f4(recs)
    rows = Vector{Any}[]
    for t in (:ladder, :direction, :estimator, :diagnostics, :excluded)
        app = [r for r in recs if t in r.applicable || t === :excluded]
        k = count(r -> t in r.label.forced, app)
        push!(rows, Any[string(t), length(app), k, isempty(app) ? NaN : k / length(app)])
    end
    raw = count(r -> r.label.raw.instability == UNRESOLVED || r.label.raw.structure == UNRESOLVED, recs)
    fin = count(r -> r.label.instability == UNRESOLVED || r.label.structure == UNRESOLVED, recs)
    push!(rows, Any["unresolved before the tests", length(recs), raw, raw / length(recs)], Any["unresolved after the tests", length(recs), fin, fin / length(recs)])
    return (columns = ["test", "applicable", "forced", "share"], rows = rows)
end

# F5: quantiles over the orbits of the stored running maxima of each diagnostic, with the threshold in force.
function _f5(recs, idx, p)
    rows = Vector{Any}[]
    for m in V_METHODS
        its = [_primary(idx, r.id, m, p.N_1, p) for r in recs]
        for nm in sort!(collect(keys(first(its).diag)))
            M = reduce(hcat, [it.diag[nm] for it in its])
            thr = _quantile([it.thresholds[nm] for it in its], 0.5)
            for (s, n) in enumerate(first(its).t)
                push!(rows, Any[m, nm, n, _quantile(M[s, :], 0.5), _quantile(M[s, :], 0.9), maximum(M[s, :]), thr])
            end
        end
    end
    return (columns = ["method", "diagnostic", "n", "q50", "q90", "max", "threshold"], rows = rows)
end

# F6: thread-nanoseconds per orbit-step by method and precision, and relative to L5 of the same precision.
function _f6(items)
    g = Dict{Tuple{String, String}, Vector{Float64}}()
    for it in items
        length(it.cost) >= 4 && it.cost[1] > 0 && push!(get!(g, (it.method, it.precision), Float64[]), 1e9 * it.cost[4] / it.cost[1])
    end
    mean_of(k) = sum(g[k]) / length(g[k])
    rows = [Any[m, pr, length(v), mean_of((m, pr)), haskey(g, ("L5", pr)) ? mean_of((m, pr)) / mean_of(("L5", pr)) : NaN] for ((m, pr), v) in sort!(collect(g))]
    return (columns = ["method", "precision", "items", "ns_per_orbit_step", "relative_to_L5"], rows = rows)
end

function _factor(it, p)
    it.tau != 1 && return ("tau", string(it.tau))
    it.eps !== nothing && it.eps != p.h1c_eps && return ("eps", string(it.eps))
    it.direction != "pi/4" && return ("direction", it.direction)
    it.mode == "wrapped" && it.method != "L2" && return ("wrapped", "wrapped")
    return nothing
end

# F7: effect of tau, eps, direction and coordinate mode on lambda(n_H2) and on n_conv, paired with the primary item, over the kept S_D orbits.
function _f7(items, idx, recs, sets, p)
    byid = Dict(r.id => r for r in recs)
    skip = Set(d.id for d in sets.dropped)
    g = Dict{Tuple{String, String, String}, Vector{Tuple{Float64, Float64}}}()
    br = Dict{Tuple{String, String, String}, Int}()
    detail = Vector{Any}[]
    for it in items
        it.horizon == p.N_1 && it.precision == "Float64" && haskey(byid, it.orbit) && byid[it.orbit].in_SD && !(it.orbit in skip) || continue
        f = _factor(it, p)
        f === nothing && continue
        key = (f[1], f[2], it.method)
        if "breakdown_tau" in it.flags
            br[key] = get(br, key, 0) + 1
            push!(detail, Any[it.orbit, it.method, f[1], f[2], get(it.attrs, "breakdown_step", "-"), get(it.attrs, "breakdown_kind", "-")])
            continue
        end
        r, prim = byid[it.orbit], _primary(idx, it.orbit, it.method, p.N_1, p)
        c1, c0 = convergence_time(it.t, it.lam, r.lam_ref, p.h3_band), convergence_time(prim.t, prim.lam, r.lam_ref, p.h3_band)
        push!(get!(g, key, Tuple{Float64, Float64}[]), (_val(it, p.n_H2) - _val(prim, p.n_H2), log10(c1.n / c0.n)))
    end
    for k in keys(br)
        haskey(g, k) || (g[k] = Tuple{Float64, Float64}[])
    end
    med(v, j) = isempty(v) ? NaN : _quantile([x[j] for x in v], 0.5)
    rows = [Any[k[1], k[2], k[3], length(v), med(v, 1), med(v, 2), get(br, k, 0)] for (k, v) in sort!(collect(g))]
    extra = isempty(detail) ? String[] : vcat(["", "Items with an orthonormalisation breakdown (a result of the sweep, A5 item 2):", ""],
                                              md_table(["orbit", "method", "factor", "level", "breakdown_step", "breakdown_kind"], sort!(detail; by = r -> (r[3], string(r[4]), r[2], r[1]))))
    return (columns = ["factor", "level", "method", "orbits", "median_dlambda_nH2", "median_log10_nconv_ratio", "breakdown_tau"], rows = rows, extra = extra)
end

"""
    figure_data(orbits, recs, items, idx, sets, p, seed, threads) -> Dict

The seven tables behind F1 to F7, each `(columns, rows)`.
"""
function figure_data(orbits, recs, items, idx, sets, p, seed, threads)
    return Dict{String, Any}("F1" => _f1(orbits, recs, idx, p), "F2" => _f2(orbits, recs, sets, p, seed, threads), "F3" => _f3(orbits, recs, sets),
                             "F4" => _f4(recs), "F5" => _f5(recs, idx, p), "F6" => _f6(items), "F7" => _f7(items, idx, recs, sets, p))
end
