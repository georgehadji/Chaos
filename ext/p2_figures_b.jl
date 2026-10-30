# D91: the figures F4 to F6 of study p2-indicators (E42 to E44) and the driver. Included by ChaoslabCairoMakieExt.jl after p2_figures.jl (helpers, constants).

p2gap(v) = v === nothing || v == "nothing"      # an AUC the block marks as not computed
p2kval(K) = K isa Real ? Float64(K) : parse(Float64, string(K))
p2ktext(K) = K isa AbstractFloat && isfinite(K) ? string(round(K; sigdigits = 12)) : string(K)   # the digits of the block (12 significant)

# ---------------------------------------------------------------- F4 (E42)
function p2f4_rows(t)
    rows = p2rows(t, "F4", ("indicator", "k", "K", "source", "n_C", "n_R", "AUC"), function (r, row, c)
        row[c["source"]] == "pooled" || return nothing       # panel and grid rows stay in the block (P37)
        auc = p2gap(row[c["AUC"]]) ? NaN : p2num("F4", r, "AUC", row[c["AUC"]]; lo = 0, hi = 1)
        return (x = p2indicator("F4", r, row[c["indicator"]], vcat(P2_I8, "STR")), k = p2num("F4", r, "k", row[c["k"]]; lo = 0), K = row[c["K"]],
                nC = p2num("F4", r, "n_C", row[c["n_C"]]; lo = 0), nR = p2num("F4", r, "n_R", row[c["n_R"]]; lo = 0), auc = auc)
    end)
    rows = filter(!isnothing, rows)
    keys_ = [(r.x, r.k, r.K) for r in rows]
    length(unique(keys_)) == length(keys_) || throw(ErrorException("F4: an indicator, budget index and K appear in more than one pooled row"))
    return rows
end

function p2f4_panel!(ax, sel, ks)
    if all(r -> isnan(r.auc), sel)
        text!(ax, sum(ks) / length(ks), 0.5; text = "AUC not computed:\na class has fewer\norbits than the\nminimum set size", align = (:center, :center), fontsize = 16)
        return
    end
    hlines!(ax, 0.5; color = :gray40, linestyle = :dash, linewidth = 1)
    for (rk, x) in enumerate(vcat(P2_I8, "STR"))
        y = [(j = findfirst(r -> r.x == x && r.k == k, sel); j === nothing ? NaN : sel[j].auc) for k in ks]
        any(isfinite, y) && p2series!(ax, ks, y, x, rk)
    end
end

# The fourth title line: the indicators with no AUC in the panel (an AUC the block marks as not computed leaves a gap in its line)
function p2f4_missing(sel, inds)
    miss = [x for x in inds if !any(r -> r.x == x && isfinite(r.auc), sel)]
    isempty(miss) && return "AUC of every indicator"
    length(miss) == length(inds) && return "no AUC at all"
    return length(miss) > 3 ? "no AUC: $(length(miss)) of $(length(inds))" : "no AUC: " * join(miss, ", ")
end

function p2f4(res, t, out)
    rows = p2f4_rows(t)
    inds = [x for x in vcat(P2_I8, "STR") if any(r.x == x for r in rows)]
    Ks = sort(unique([r.K for r in rows]); by = p2kval); ks = sort(unique([r.k for r in rows]))
    nr = cld(length(Ks), 5)
    fig = Figure(size = (1000, 380 * nr + 110), figure_padding = 12)
    for (i, K) in enumerate(Ks)
        sel = [r for r in rows if r.K == K]
        ax = Axis(fig[fld(i - 1, 5) + 1, mod(i - 1, 5) + 1]; title = "K = $(p2ktext(K))\nn_C $(p2range([r.nC for r in sel]))\nn_R $(p2range([r.nR for r in sel]))\n$(p2f4_missing(sel, inds))", xticks = (ks, p2klabels(res, ks)),
                  xticklabelrotation = pi / 2, xlabel = "budget b_k", ylabel = mod(i - 1, 5) == 0 ? "AUC" : "", limits = ((-0.4, maximum(ks) + 0.4), (-0.03, 1.03)))
        mod(i - 1, 5) > 0 && (ax.yticklabelsvisible = false)
        p2f4_panel!(ax, sel, ks)
    end
    el, nm = p2legend_entries(inds)
    Legend(fig[nr + 1, 1:min(5, length(Ks))], vcat(el, [[LineElement(color = :gray40, linestyle = :dash, linewidth = 1)]]), vcat(nm, "chance level (AUC = 0.5)"); orientation = :horizontal, nbanks = 2, tellwidth = false)
    return p2save(fig, out, "F4")
end

# ---------------------------------------------------------------- F5 (E43)
const P2_STATE_WORD = Dict("C" => :C, "R" => :R, "U" => :U, "chaotic" => :C, "regular" => :R, "undecided" => :U)
const P2_BUDGET_COLS = (("b_P", :b_primary), ("N_c", :N_c))

function p2f5_rows(t)
    vcols = [["$x $b" for x in P2_I8] for (b, _) in P2_BUDGET_COLS]
    state(r, v) = haskey(P2_STATE_WORD, v) ? P2_STATE_WORD[v] : p2fail("F5", r, "state $v is not C, R, U, chaotic, regular or undecided")
    int(r, name, v) = (x = p2num("F5", r, name, v; lo = 0); isinteger(x) ? Int(x) : p2fail("F5", r, "$name value $v is not an integer"))
    return p2rows(t, "F5", vcat(["K index", "i", "j", "rho"], vcols...), (r, row, c) -> (kx = int(r, "K index", row[c["K index"]]), i = int(r, "i", row[c["i"]]), j = int(r, "j", row[c["j"]]),
                  rho = state(r, row[c["rho"]]), st = [[state(r, row[c[n]]) for n in vc] for vc in vcols]))
end

# One colour matrix per (group, K index, column); column 1 is rho, columns 2 to 9 the indicators. A position without a row stays white.
function p2f5_maps(rows)
    ni, nj = maximum(r.i for r in rows) + 1, maximum(r.j for r in rows) + 1
    kxs = sort(unique([r.kx for r in rows]))
    maps = Dict((g, kx, c) => fill(to_color("white"), ni, nj) for g in 1:2, kx in kxs, c in 1:9)
    seen = Set{NTuple{3, Int}}()
    for r in rows
        key = (r.kx, r.i, r.j)
        key in seen && throw(ErrorException("F5: K index $(r.kx) holds the start i = $(r.i), j = $(r.j) twice"))
        push!(seen, key)
        for g in 1:2
            maps[(g, r.kx, 1)][r.i + 1, r.j + 1] = P2_STATE_COLOR[r.rho]
            for c in 2:9
                maps[(g, r.kx, c)][r.i + 1, r.j + 1] = P2_STATE_COLOR[r.st[g][c - 1]]
            end
        end
    end
    return maps, kxs, ni, nj
end

# Page geometry in pt: one data unit of the single axis is one pt, a cell is s x s (square by construction, no interpolation).
function p2f5_geometry(nk, ni, nj)
    s = min(4.0, 830 / (9 * ni)); mw, mh, gx, gy, x0, head = ni * s, nj * s, 6.0, 4.0, 100.0, 72.0
    gh = head + nk * (mh + gy) + 16
    return (s = s, mw = mw, mh = mh, gx = gx, gy = gy, x0 = x0, head = head, gh = gh, W = ceil(x0 + 9 * mw + 8 * gx + 16), H = ceil(2 * gh + 44))
end

function p2f5_group!(ax, res, g, kxs, maps, L, ni, nj)
    top = 10 + (g - 1) * L.gh
    Y(yd) = L.H - yd
    budget = P2_BUDGET_COLS[g]
    v = p2param(res, budget[2])
    text!(ax, 12, Y(top); text = "verdict of each indicator at budget $(budget[1])$(v === nothing ? "" : " = $v")", align = (:left, :top), fontsize = 16)
    text!(ax, 12, Y(top + 24); text = "each map is the lattice of its K index: i = 0 to $(ni - 1) along x, j = 0 to $(nj - 1) along y (up); x0: 0 to 2 pi, p0: 0 to pi", align = (:left, :top), fontsize = 16)
    for (c, name) in enumerate(vcat("rho", P2_I8))
        text!(ax, L.x0 + (c - 1) * (L.mw + L.gx) + L.mw / 2, Y(top + 48); text = name, align = (:center, :top), fontsize = 16)
    end
    for (r, kx) in enumerate(kxs)
        yb = top + L.head + (r - 1) * (L.mh + L.gy)
        text!(ax, L.x0 - 8, Y(yb + L.mh / 2); text = "K index $kx", align = (:right, :center), fontsize = 16)
        for c in 1:9
            xl = L.x0 + (c - 1) * (L.mw + L.gx)
            image!(ax, xl .. (xl + L.mw), Y(yb + L.mh) .. Y(yb), maps[(g, kx, c)]; interpolate = false)
            poly!(ax, Rect2f(xl, Y(yb + L.mh), L.mw, L.mh); color = :transparent, strokecolor = :gray60, strokewidth = 0.5)
        end
    end
end

function p2f5(res, t, out)
    rows = p2f5_rows(t)
    maps, kxs, ni, nj = p2f5_maps(rows)
    L = p2f5_geometry(length(kxs), ni, nj)
    fig = Figure(size = (L.W, L.H), figure_padding = 0)
    ax = Axis(fig[1, 1]; limits = (0, L.W, 0, L.H))
    hidedecorations!(ax); hidespines!(ax)
    for g in 1:2
        p2f5_group!(ax, res, g, kxs, maps, L, ni, nj)
    end
    for (n, (st, word)) in enumerate((:C => "C chaotic", :R => "R regular", :U => "U undecided"))
        x = 12 + (n - 1) * 150
        poly!(ax, Rect2f(x, 12, 16, 16); color = P2_STATE_COLOR[st], strokecolor = :gray40, strokewidth = 0.5)
        text!(ax, x + 22, 20; text = word, align = (:left, :center), fontsize = 16)
    end
    text!(ax, 470, 20; text = "rho: reference class; the eight columns: verdicts", align = (:left, :center), fontsize = 16)
    return p2save(fig, out, "F5")
end

# ---------------------------------------------------------------- F6 (E44)
const P2_CLASS_STYLE = Dict("C" => (P2_STATE_COLOR[:C], :pentagon), "R" => (P2_STATE_COLOR[:R], :star4), "T_s" => (to_color("#5E5E5E"), :ltriangle))

function p2pow10(vals)   # powers of ten in plain digits: a superscript would be set below 8 pt at print width
    ks = floor(Int, log10(minimum(vals))):ceil(Int, log10(maximum(vals)))
    return (10.0 .^ ks, [k >= 0 ? string(10^k) : "0." * "0"^(-k - 1) * "1" for k in ks])
end

function p2f6_name(cl, n)
    length(unique(n)) == 1 && return "$cl ($(p2fmt(n[1])) orbits)"
    return "$cl ($(p2fmt(n[1])) orbits at step 1;\nthe count changes, $(p2range(n)))"
end

function p2f6(res, t, out)
    qs = ["p5", "p25", "p50", "p75", "p95"]
    rows = p2rows(t, "F6", vcat(["class", "step", "count"], qs), function (r, row, c)
        cl = String(row[c["class"]])
        haskey(P2_CLASS_STYLE, cl) || p2fail("F6", r, "class $cl is not C, R or T_s")
        return (cl = cl, step = p2num("F6", r, "step", row[c["step"]]; pos = true), n = p2num("F6", r, "count", row[c["count"]]; lo = 0),
                q = [p2num("F6", r, q, row[c[q]]; pos = true) for q in qs])
    end)
    fig = Figure(size = (1000, 560), figure_padding = 12)
    ax = Axis(fig[1, 1]; xscale = log10, yscale = log10, xlabel = "recorded step n", ylabel = "D_str", title = "stretching distance D_str against n, by class",
              xticks = p2pow10([r.step for r in rows]), yticks = p2pow10(vcat([r.q for r in rows]...)))
    el, nm = Any[], String[]
    for cl in p2order([r.cl for r in rows], ["C", "R", "T_s"])
        sel = [r for r in rows if r.cl == cl]; x = [r.step for r in sel]; q(i) = [r.q[i] for r in sel]; (col, mk) = P2_CLASS_STYLE[cl]
        band!(ax, x, q(1), q(5); color = (col, 0.14)); band!(ax, x, q(2), q(4); color = (col, 0.3))
        lines!(ax, x, q(3); color = col, linewidth = 2, linestyle = cl == "T_s" ? :dash : :solid); scatter!(ax, x, q(3); color = col, marker = mk, markersize = 6)
        push!(el, [LineElement(color = col, linewidth = 2, linestyle = cl == "T_s" ? :dash : :solid), MarkerElement(color = col, marker = mk, markersize = 10)]); push!(nm, p2f6_name(cl, [r.n for r in sel]))
    end
    push!(el, [PolyElement(color = (:gray30, 0.3))], [PolyElement(color = (:gray30, 0.14))]); push!(nm, "band p25 to p75", "band p5 to p95")
    Legend(fig[1, 2], el, nm, "line: median p50")
    Label(fig[2, 1:2], "T_s is a subset of the class C: the two share orbits"; halign = :left, tellwidth = false)
    return p2save(fig, out, "F6")
end

# All six pages are drawn to temporary files first; a failure in any figure leaves the output folder as it was (the files of an earlier run stay whole) and names the figure.
function draw_p2_figures(res, out)
    f = res.figures
    tmps = String[]
    try
        with_theme(P2_THEME) do
            for (id, draw) in (("F1", p2f1), ("F2", p2f2), ("F3", p2f3), ("F4", p2f4), ("F5", p2f5), ("F6", p2f6))
                haskey(f, id) || throw(ErrorException("$id: the result holds no figure table $id"))
                try
                    push!(tmps, draw(res, f[id], out))
                catch e
                    msg = e isa ErrorException ? e.msg : sprint(showerror, e)
                    startswith(msg, id) ? rethrow() : throw(ErrorException("$id: $msg"))
                end
            end
        end
        paths = [chopsuffix(t, ".tmp.pdf") * ".pdf" for t in tmps]
        foreach((t, p) -> mv(t, p; force = true), tmps, paths)
        return paths
    finally
        foreach(t -> rm(t; force = true), tmps)
    end
end
