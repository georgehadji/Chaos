# Tables of the panel file: periodic orbits, controls, orbit list with ids, subsets and the cell records (spec Algorithm 3, 7, 8, 11).

const FORMAT_VERSION = 1
const STREAM_RULE = "first 8 bytes, big-endian, of SHA-256 over the UTF-8 text \"chaoslab-panel-1|<seed>|<K index>|<class>\"; StableRNGs.StableRNG"

_dec(x) = string(BigFloat(x))

"Rows of the periodic table over all K (ids 1, 2, ... in K order) and the id of every orbit per K."
function _periodic_rows(psets)
    rows = Dict{String, Any}[]
    pids = [Int[] for _ in psets]
    for (i, ps) in enumerate(psets)
        for (a, b) in ps.lines              # K = 0: the line p = 2 pi a / b, as the orbit of (0, p) with windings (a, 0) (E7)
            p = 2 * BigFloat(pi) * a / b
            pts = [[string(_wrap_big(t * p)), string(_wrap_big(p))] for t in 0:(b - 1)]
            push!(rows, Dict{String, Any}("id" => length(rows) + 1, "k_index" => i, "q" => b, "m" => a, "j" => 0, "class" => "parabolic",
                                           "trace" => 2.0, "exponent" => 0.0, "residual" => 0.0, "points" => pts))
        end
        for o in ps.orbits
        push!(rows, Dict{String, Any}("id" => length(rows) + 1, "k_index" => i, "q" => o.q, "m" => o.m, "j" => o.j,
                                       "class" => String(o.class), "trace" => Float64(o.trace), "exponent" => Float64(o.exponent),
                                       "residual" => o.residual, "points" => [[string(pt[1]), string(pt[2])] for pt in o.points]))
            push!(pids[i], length(rows))
        end
    end
    return rows, pids
end

function _control_dict(inp, i, o, kind)
    pt = o.points[1]
    d = Dict{String, Any}("class" => "control", "control" => true, "control_kind" => kind, "system" => "standard", "k_index" => i,
                          "K_text" => _ktext(inp.raw["K"][i]), "K" => inp.K[i], "run_seed" => string(inp.seed),
                          "in_SD" => true, "in_SL" => false, "period" => o.q, "wind_m" => o.m, "wind_j" => o.j,
                          "stability" => String(o.class), "exponent" => Float64(o.exponent))
    return _with_start!(d, pt[1], pt[2])
end

_startkey(d) = (d["x0"], d["p0"])

"""
Controls (E19): the exact fixed points by K index (0, 0) before (pi, 0), then the Newton orbits by (period, K index,
representative); the first `n_ctrl_max` of that list, then the Y5 control.
"""
function _controls(inp, psets)
    fixed, newton = Dict{String, Any}[], Dict{String, Any}[]
    for (i, ps) in enumerate(psets), o in ps.orbits
        isfix = o.grid_id == -1 && o.q == 1
        push!(isfix ? fixed : newton, _control_dict(inp, i, o, isfix ? "fixed" : "newton"))
    end
    sort!(fixed; by = d -> (d["k_index"], _startkey(d)))
    sort!(newton; by = d -> (d["period"], d["k_index"], _startkey(d)))
    list = first(vcat(fixed, newton), inp.n_ctrl_max)
    y = y5_control(inp.mu)
    push!(list, Dict{String, Any}("class" => "control", "control" => true, "control_kind" => "y5", "system" => "linear", "k_index" => 0,
                                  "run_seed" => string(inp.seed), "in_SD" => true, "in_SL" => false, "exponent" => y.exponent,
                                  "mu" => inp.mu, "x0" => 0.0, "p0" => 0.0, "x0_dec" => "0.0", "p0_dec" => "0.0"))
    return list
end

"All orbits in the order of E20 with ids, per-orbit seeds and subset flags set; returns the orbit dictionaries and the subset record."
function _number_orbits!(inp, cells, controls)
    orbs = Dict{String, Any}[]
    members = NamedTuple{(:id, :cell, :rank, :control), Tuple{Int, Int, Int, Bool}}[]
    for (c, cell) in enumerate(cells), d in cell.orbits
        push!(orbs, d)
        push!(members, (id = length(orbs), cell = c, rank = d["draw_index"], control = false))
    end
    for d in controls
        push!(orbs, d)
        push!(members, (id = length(orbs), cell = 0, rank = 0, control = true))
    end
    for (id, d) in enumerate(orbs)
        d["id"] = id
        d["orbit_seed"] = string(orbit_seed(inp.seed, id))
    end
    sub = select_subsets(members; s_D = inp.s_D, cap_D = inp.cap_D, s_L = inp.s_L, cap_L = inp.cap_L)
    for id in sub.SD
        orbs[id]["in_SD"] = true
    end
    for id in sub.SL
        orbs[id]["in_SL"] = true
    end
    return orbs, sub
end

"Cell records (K index and class order) with counts, status and reason."
function _cell_rows(inp, all_cells)
    rows = Dict{String, Any}[]
    for (i, cells) in enumerate(all_cells), c in cells
        undefined = hasproperty(c, :undefined)
        status, reason = undefined ? ("undefined", "class undefined at K = 0") : _status(c.class, length(c.orbits), inp.n_c)
        push!(rows, Dict{String, Any}("k_index" => i, "class" => c.class, "qualifying" => c.qualifying, "drawn" => length(c.orbits),
                                       "status" => status, "reason" => reason))
    end
    return rows
end

_sum_stats(psets) = Dict{String, Any}(String(k) => sum(getproperty(ps.stats, k) for ps in psets) for k in keys(_NO_STATS))

function _counts(psets)
    rows = Dict{String, Any}[]
    for (i, ps) in enumerate(psets), q in sort!(unique(o.q for o in ps.orbits)), c in (:elliptic, :parabolic, :hyperbolic)
        n = count(o -> o.q == q && o.class === c, ps.orbits)
        n > 0 && push!(rows, Dict{String, Any}("k_index" => i, "q" => q, "class" => String(c), "orbits" => n))
    end
    return rows
end

function _inputs_table(inp)
    t = Dict{String, Any}()
    for (k, v) in inp.raw
        t[k] = Dict{String, Any}("value" => v, "source" => "preregistration")
    end
    for (k, v) in (("tau", inp.tau), ("seed", inp.seed), ("bigfloat_bits", inp.bigfloat_bits))
        haskey(t, k) || (t[k] = Dict{String, Any}("value" => v, "source" => inp.sources[k]))
    end
    return t
end


"""
Ill-conditioned rejections of the search (Amendment D51) per K index and period: count and the smallest `|2 - t|` among them.
"""
function _ill_rows(psets)
    return [Dict{String, Any}("k_index" => i, "q" => r.q, "count" => r.count, "min_abs_2_minus_t" => r.min_den)
            for (i, ps) in enumerate(psets) for r in ps.ill]
end
