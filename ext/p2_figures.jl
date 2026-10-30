# D91: the figures F1 to F3 of study p2-indicators (E38 to E41) and the helpers shared by F1 to F6 (this file, p2_figures_b.jl).
# Included by ChaoslabCairoMakieExt.jl; imports nothing beyond Chaoslab and CairoMakie (I29). Every text is 16 units on a page of at most 1000 units drawn with
# pt_per_unit = 1, so the smallest text scaled to 504 pt is at least 8 pt (E38 e). A value that cannot be drawn raises an ErrorException that names figure and row.

const P2_I8 = ["LE", "SALI", "FLI", "MEGNO", "E1", "E2", "WBA", "RLI"]
const P2_SET_ORDER = ["W2", "S2", "T_s", "regular"]   # the name of the regular set in the block of F3
# chaotic/C orange, regular/R dark blue, undecided/U light grey: three steps in lightness (grey-scale print), orange/blue is safe for the common colour-vision deficiencies
const P2_STATE_COLOR = Dict(:C => to_color("#E69F00"), :R => to_color("#1B4F9C"), :U => to_color("#E2E2E2"))
# indicators: eight colours away from the class colours, plus one marker shape each (the shape carries the identity when the colours merge in grey)
const P2_IND_COLOR = Dict(zip(vcat(P2_I8, "STR"), to_color.(["#882255", "#117733", "#CC6677", "#999933", "#44AA99", "#AA4499", "#88CCEE", "#8C510A", "#000000"])))
const P2_IND_MARKER = Dict(zip(vcat(P2_I8, "STR"), [:circle, :rect, :diamond, :utriangle, :dtriangle, :hexagon, :star5, :cross, :xcross]))
const P2_THEME = Theme(fontsize = 16, Axis = (titlesize = 16, xlabelsize = 16, ylabelsize = 16, xticklabelsize = 16, yticklabelsize = 16),
                       Colorbar = (labelsize = 16, ticklabelsize = 16), Legend = (labelsize = 16, titlesize = 16, framevisible = false))

p2fail(fig, r, msg) = throw(ErrorException("$fig row $r: $msg"))

function p2col(t, name, fig)
    j = findfirst(==(name), t.columns)
    j === nothing && throw(ErrorException("$fig: the block has no column $name"))
    return j
end

# One number of a block, checked for drawing: finite, positive on a logarithmic axis, inside the fixed range of its axis (E38 b: never altered, never clipped).
function p2num(fig, r, name, v; lo = -Inf, hi = Inf, pos = false)
    v isa Real && isfinite(v) || p2fail(fig, r, "$name value $v cannot be drawn")
    pos && v <= 0 && p2fail(fig, r, "$name value $v cannot be drawn on a logarithmic axis")
    lo - 1e-9 <= v <= hi + 1e-9 || p2fail(fig, r, "$name value $v is outside the axis range [$lo, $hi]")
    return Float64(v)
end

# The values of the named columns of every row, each row checked by `check(r, row, c)`; c maps a column name to its index.
function p2rows(t, fig, names, check)
    c = Dict(n => p2col(t, n, fig) for n in names)
    return [check(r, row, c) for (r, row) in enumerate(t.rows)]
end

function p2indicator(fig, r, x, allowed)
    x in allowed || p2fail(fig, r, "indicator $x is not one of $(join(allowed, ", "))")
    return String(x)
end

function p2param(res, key)
    p = hasproperty(res, :params) ? res.params : nothing
    p isa AbstractDict && return get(p, key, get(p, string(key), nothing))
    return p !== nothing && hasproperty(p, key) ? getproperty(p, key) : nothing
end

# A budget is called b_P or N_c only when it equals that key of the parameter table (E38 c).
function p2budget(res, v)
    v == p2param(res, :b_primary) && return "b_P = $v"
    v == p2param(res, :N_c) && return "N_c = $v"
    return "budget $v"
end

# Tick labels of the budget index k: the values b_k when the parameter table is at hand, else the index.
function p2klabels(res, ks)
    b = p2param(res, :budgets)
    return [b !== nothing && k + 1 <= length(b) ? string(b[Int(k) + 1]) : string(Int(k)) for k in ks]
end

p2order(vals, lead) = (u = unique(vals); vcat([s for s in lead if s in u], [s for s in u if !(s in lead)]))
p2fmt(v) = isinteger(v) ? string(Int(v)) : string(v)
p2range(v) = (lo = minimum(v); hi = maximum(v); lo == hi ? p2fmt(lo) : "$(p2fmt(lo)) to $(p2fmt(hi))")

# lines with a marker at every finite point (E38 d); NaN is a gap; marker size falls with the rank so that coinciding points stay visible
function p2series!(ax, x, y, name, rank; steps = false)
    kw = (color = P2_IND_COLOR[name], linestyle = name == "STR" ? :dash : :solid, linewidth = 2)
    steps ? stairs!(ax, x, y; step = :post, kw...) : lines!(ax, x, y; kw...)
    ok = isfinite.(y)
    scatter!(ax, x[ok], y[ok]; color = P2_IND_COLOR[name], marker = P2_IND_MARKER[name], markersize = 15 - rank)
end

function p2legend_entries(names)
    elems = [[LineElement(color = P2_IND_COLOR[x], linestyle = x == "STR" ? :dash : :solid, linewidth = 2),
              MarkerElement(color = P2_IND_COLOR[x], marker = P2_IND_MARKER[x], markersize = 12)] for x in names]
    return elems, collect(names)
end

# The page must sit inside the bounds of E38 e. It is drawn to a temporary file; draw_p2_figures moves the six into place only when all six are drawn (E38 b: whole or not at all).
function p2save(fig, out, id)
    w, h = size(fig.scene)
    (w <= 1000 && h <= 1.3 * w) || throw(ErrorException("$id: the page $w x $h pt is outside the bounds (width 1000, height 1.3 x width)"))
    tmp = joinpath(out, "$(Analysis.P2_SLUG)-fig-$id.tmp.pdf")
    try
        save(tmp, fig; pt_per_unit = 1)
    catch
        rm(tmp; force = true)
        rethrow()
    end
    return tmp
end

# ---------------------------------------------------------------- F1 (E39)
function p2f1_rows(t)
    return p2rows(t, "F1", ("set", "budget", "X", "Y", "N", "A", "kappa"), function (r, row, c)
        x, y = findfirst(==(row[c["X"]]), P2_I8), findfirst(==(row[c["Y"]]), P2_I8)
        (x === nothing || y === nothing || x == y) && p2fail("F1", r, "the pair $(row[c["X"]]), $(row[c["Y"]]) is not two different indicators of I_8")
        p2num("F1", r, "budget", row[c["budget"]]; pos = true)
        return (set = String(row[c["set"]]), budget = row[c["budget"]], i = min(x, y), j = max(x, y), N = p2num("F1", r, "N", row[c["N"]]; lo = 0),
                A = p2num("F1", r, "A", row[c["A"]]; lo = 0, hi = 1), kappa = p2num("F1", r, "kappa", row[c["kappa"]]; lo = -1, hi = 1))
    end)
end

function p2f1_matrix(rows, set, budget, field)
    M = fill(NaN, 8, 8)
    for r in rows
        (r.set == set && r.budget == budget) || continue
        isnan(M[r.i, r.j]) || throw(ErrorException("F1: the pair $(P2_I8[r.i]), $(P2_I8[r.j]) appears twice for $set at budget $budget"))
        M[r.i, r.j] = getfield(r, field)
    end
    return M
end

function p2f1_note(rows)
    rli = [(r.set, r.N) for r in rows if r.i == 8 || r.j == 8]
    isempty(rli) && return "*RLI exists for a subset of the orbits; this block holds no pair with RLI"
    return "*RLI exists for a subset of the orbits, N = " * join(["$(p2fmt(n)) in $s" for (s, n) in unique(rli)], ", ")
end

function p2f1(res, t, out)
    rows = p2f1_rows(t)
    sets = p2order([r.set for r in rows], ["W2", "S2"])
    cols = [(s, b) for s in sets for b in sort(unique([r.budget for r in rows if r.set == s]))]
    ps = min(185, (1000 - 240) / length(cols) - 12)
    names = [x == "RLI" ? "RLI*" : x for x in P2_I8]
    fig = Figure(figure_padding = 12)
    for (c, (s, b)) in enumerate(cols)
        ns = [r.N for r in rows if r.set == s && r.budget == b && r.j != 8]
        ns = isempty(ns) ? [r.N for r in rows if r.set == s && r.budget == b] : ns
        for (ri, field, cm, lim) in ((1, :A, :viridis, (0, 1)), (2, :kappa, :vik, (-1, 1)))
            ax = Axis(fig[ri, c]; width = ps, height = ps, yreversed = true, xticks = (1:8, names), yticks = (1:8, names), xticklabelrotation = pi / 2, xgridvisible = false, ygridvisible = false,
                      xlabel = ri == 2 ? "indicator" : "", ylabel = c == 1 ? "indicator" : "", title = ri == 1 ? "$s, $(p2budget(res, b))\nN = $(p2range(ns))" : "")
            hm = heatmap!(ax, 1:8, 1:8, permutedims(p2f1_matrix(rows, s, b, field)); colormap = cm, colorrange = lim, nan_color = :transparent)
            c > 1 && (ax.yticklabelsvisible = false)
            c == length(cols) && Colorbar(fig[ri, length(cols) + 1], hm; label = field == :A ? "A_XY" : "kappa_XY", height = ps)
        end
    end
    Label(fig[3, 1:length(cols)], p2f1_note(rows); halign = :left, tellwidth = false)
    resize_to_layout!(fig)
    return p2save(fig, out, "F1")
end

# ---------------------------------------------------------------- F2 (E40)
const P2_OUTCOMES = ["correct", "wrong", "undetected"]
const P2_OUTCOME_GREY = [to_color("#2B2B2B"), to_color("#8C8C8C"), to_color("#D9D9D9")]

# t_det counts by (set, indicator) as (budget, count) pairs, and the three outcome counts by (set, indicator)
function p2f2_data(t)
    det, outc = Dict{Tuple{String, String}, Vector{Tuple{Float64, Float64}}}(), Dict{Tuple{String, String}, Dict{String, Float64}}()
    seen = Set{Tuple}()
    p2rows(t, "F2", ("set", "indicator", "kind", "value", "count"), function (r, row, c)
        key = (String(row[c["set"]]), p2indicator("F2", r, row[c["indicator"]], P2_I8))
        n, kind = p2num("F2", r, "count", row[c["count"]]; lo = 0), row[c["kind"]]
        id = (key, kind, string(row[c["value"]]))
        id in seen && p2fail("F2", r, "the row $(key[1]) $(key[2]) $kind $(row[c["value"]]) appears twice")
        push!(seen, id)
        kind == "t_det" && push!(get!(det, key, Tuple{Float64, Float64}[]), (p2num("F2", r, "t_det value", row[c["value"]]; pos = true), n))
        kind == "outcome" && (get!(outc, key, Dict{String, Float64}())[string(row[c["value"]])] = n)
        kind in ("t_det", "outcome", "t_censored") || p2fail("F2", r, "kind $kind is not t_det, t_censored or outcome")
    end)
    return det, outc
end

function p2f2_share(key, det, outc)
    haskey(outc, key) && all(haskey(outc[key], o) for o in P2_OUTCOMES) || throw(ErrorException("F2: $(key[1]) $(key[2]) has no complete outcome rows"))
    N = sum(values(outc[key]))
    N > 0 || throw(ErrorException("F2: $(key[1]) $(key[2]) has no orbit (the three outcome counts sum to 0)"))
    d = sort(det[key])
    return first.(d), cumsum(last.(d)) ./ N, N
end

function p2f2(res, t, out)
    det, outc = p2f2_data(t)
    keys_ = collect(keys(det)); sets = p2order(first.(keys_), ["S2", "W2", "T_s"])
    xs = sort(unique(first.(reduce(vcat, values(det)))))
    fig = Figure(size = (1000, 900), figure_padding = 12)
    for (c, s) in enumerate(sets)
        inds = [x for x in P2_I8 if (s, x) in keys_]
        sh = Dict(x => p2f2_share((s, x), det, outc) for x in inds)
        top = Axis(fig[1, c]; title = "$s (N = $(p2fmt(maximum(last(sh[x]) for x in inds))))", xscale = log10, xticks = (xs, p2fmt.(xs)), xticklabelrotation = pi / 2,
                   xlabel = "t_det", ylabel = c == 1 ? "share of N, t_det at most the budget" : "", limits = (nothing, (-0.03, 1.03)))
        c > 1 && (top.yticklabelsvisible = false)
        for (rk, x) in enumerate(inds)
            p2series!(top, sh[x][1], sh[x][2], x, rk; steps = true)
        end
        p2f2_bars!(Axis(fig[2, c]; xlabel = "share of the N orbits", ylabel = c == 1 ? "indicator" : "", limits = ((0, 1), nothing)), inds, sh, outc, s)
    end
    Legend(fig[3, 1:length(sets)], p2legend_entries(P2_I8)..., "indicator (top row)"; orientation = :horizontal, tellwidth = false, titleposition = :left)
    Legend(fig[4, 1:length(sets)], [PolyElement(color = g) for g in P2_OUTCOME_GREY], P2_OUTCOMES, "outcome of N (bottom row)"; orientation = :horizontal, tellwidth = false, titleposition = :left)
    return p2save(fig, out, "F2")
end

function p2f2_bars!(ax, inds, sh, outc, s)
    full = maximum(last(sh[x]) for x in inds)
    pos = Float64[]; val = Float64[]; stk = Int[]
    for (rk, x) in enumerate(inds), (g, o) in enumerate(P2_OUTCOMES)
        push!(pos, rk); push!(val, outc[(s, x)][o] / last(sh[x])); push!(stk, g)
    end
    barplot!(ax, pos, val; stack = stk, direction = :x, color = P2_OUTCOME_GREY[stk])
    lab = [last(sh[x]) == full ? x : "$x (N = $(p2fmt(last(sh[x]))))" for x in inds]
    ax.yticks = (1:length(inds), lab); ax.yreversed = true
end

# ---------------------------------------------------------------- F3 (E41)
function p2f3(res, t, out)
    rows = p2rows(t, "F3", ("set", "indicator", "k", "decided", "fraction"), (r, row, c) -> (set = String(row[c["set"]]), x = p2indicator("F3", r, row[c["indicator"]], P2_I8),
                  k = p2num("F3", r, "k", row[c["k"]]; lo = 0), dec = p2num("F3", r, "decided", row[c["decided"]]; lo = 0), y = p2num("F3", r, "fraction", row[c["fraction"]]; lo = 0, hi = 1)))
    keys_ = [(r.set, r.x, r.k) for r in rows]
    length(unique(keys_)) == length(keys_) || throw(ErrorException("F3: a set, indicator and budget index appear in more than one row"))
    sets = p2order([r.set for r in rows], P2_SET_ORDER); ks = sort(unique([r.k for r in rows]))
    bud = p2param(res, :budgets)
    marks = bud === nothing ? [] : [(findfirst(==(v), bud), nm, ls) for (v, nm, ls) in ((p2param(res, :b_primary), "b_P", :dash), (p2param(res, :N_c), "N_c", :dot)) if v in bud]
    fig = Figure(size = (1000, 800), figure_padding = 12)
    for (i, s) in enumerate(sets)
        ax = p2f3_axis!(fig, i, s, rows, ks, p2klabels(res, ks))
        for (kk, _, ls) in marks
            vlines!(ax, kk - 1; color = :black, linestyle = ls, linewidth = 1.5)
        end
    end
    inds = [x for x in P2_I8 if any(r.x == x for r in rows)]
    el, _ = p2legend_entries(inds)
    nm = [p2f3_name(x, sets, rows) for x in inds]
    for (kk, n, ls) in marks
        push!(el, [LineElement(color = :black, linestyle = ls, linewidth = 1.5)]); push!(nm, "$n = $(bud[kk])")
    end
    Legend(fig[1:cld(length(sets), 2), 3], el, nm)
    return p2save(fig, out, "F3")
end

function p2f3_axis!(fig, i, s, rows, ks, labels)
    ix = findall(r -> r.set == s, rows)
    last_row = fld(i - 1, 2) + 1 == cld(length(unique(r.set for r in rows)), 2)
    note = s == "regular" ? "regular set: share of chaotic\nverdicts, a false-chaotic rate\n" : "$s: "
    ax = Axis(fig[fld(i - 1, 2) + 1, mod(i - 1, 2) + 1]; title = note * "up to $(p2fmt(maximum(rows[j].dec for j in ix))) orbits decided", xticks = (ks, labels),
              xticklabelrotation = pi / 2, xlabel = last_row ? "budget b_k" : "", ylabel = mod(i - 1, 2) == 0 ? "fraction chaotic" : "", limits = ((-0.4, maximum(ks) + 0.4), (-0.03, 1.03)))
    for (rk, x) in enumerate(P2_I8)
        y = [(j = findfirst(r -> r.set == s && r.x == x && r.k == k, rows); j === nothing ? NaN : rows[j].y) for k in ks]
        any(isfinite, y) && p2series!(ax, ks, y, x, rk)
    end
    return ax
end

# RLI exists for a subset of the orbits, so the legend gives its own decided count per set
function p2f3_name(x, sets, rows)
    x == "RLI" || return x
    own(s) = maximum([r.dec for r in rows if r.set == s && r.x == "RLI"])
    return "RLI (own decided count,\n" * join(["$s $(p2fmt(own(s)))" for s in sets if any(r.set == s && r.x == "RLI" for r in rows)], ", ") * ")"
end
