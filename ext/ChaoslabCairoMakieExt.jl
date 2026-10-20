"""
    ChaoslabCairoMakieExt

The Figures layer of the analysis (decision D60): draws F1 to F7 as PDF files `<study>-fig-F<k>.pdf` from the figure tables of `Chaoslab.Analysis.analyze`.
Loaded when CairoMakie is present (the verify environment); the analysis itself never needs it.
"""
module ChaoslabCairoMakieExt

using Chaoslab, CairoMakie
using Chaoslab.Analysis: Analysis

const VM = ["L1", "L3", "L4", "L4b", "L5", "L5b", "L6", "L9"]
const COLORS = Dict("L1" => :black, "L2" => :gray, "L3" => :red, "L4" => :orange, "L4b" => :gold, "L5" => :blue, "L5b" => :cyan, "L6" => :green, "L9" => :purple)

col(t, name) = (j = findfirst(==(name), t.columns); [r[j] for r in t.rows])

# Row indices of a figure table grouped by the values of the named columns, in order of first appearance.
function groups(t, names...)
    keys_ = [Tuple(r[findfirst(==(n), t.columns)] for n in names) for r in t.rows]
    return [(k, findall(==(k), keys_)) for k in unique(keys_)]
end

function savefig(fig, out, study, id)
    path = joinpath(out, "$study-fig-$id.pdf")
    save(path, fig)
    return path
end

function grid(n; cols = 4)
    return (cols = min(n, cols), rows = cld(n, cols))
end

function f1(t, out, study)
    gs = groups(t, "class", "k_index")
    g = grid(length(gs))
    fig = Figure(size = (320 * g.cols, 260 * g.rows))
    for (i, (key, ix)) in enumerate(gs)
        ax = Axis(fig[fld(i - 1, g.cols) + 1, mod(i - 1, g.cols) + 1]; title = "class $(key[1]), K index $(key[2])", xscale = log10, xlabel = "n", ylabel = "lambda_1(n)")
        n = col(t, "n")[ix]
        lo = [minimum(t.rows[j][findfirst(==(m), t.columns)] for m in VM) for j in ix]
        hi = [maximum(t.rows[j][findfirst(==(m), t.columns)] for m in VM) for j in ix]
        band!(ax, n, lo, hi; color = (:gray, 0.3))
        for m in vcat(VM, "L2")
            lines!(ax, n, col(t, m)[ix]; color = COLORS[m], label = m)
        end
        i == 1 && axislegend(ax; position = :rt, nbanks = 2, labelsize = 8)
    end
    return savefig(fig, out, study, "F1")
end

function f2(t, out, study)
    gs = groups(t, "class", "k_index")
    g = grid(length(gs))
    fig = Figure(size = (320 * g.cols, 260 * g.rows))
    for (i, (key, ix)) in enumerate(gs)
        ax = Axis(fig[fld(i - 1, g.cols) + 1, mod(i - 1, g.cols) + 1]; title = "class $(key[1]), K index $(key[2])", xscale = log10, xlabel = "n", ylabel = "R(n)")
        n = col(t, "n")[ix]
        band!(ax, n, col(t, "lo")[ix], col(t, "hi")[ix]; color = (:blue, 0.25))
        lines!(ax, n, col(t, "median")[ix]; color = :blue)
    end
    return savefig(fig, out, study, "F2")
end

function f3(t, out, study)
    gs = groups(t, "class", "k_index")
    fig = Figure(size = (640, 420))
    ax = Axis(fig[1, 1]; xscale = log10, xlabel = "n_conv", ylabel = "fraction of orbits", title = "convergence time")
    for (key, ix) in gs
        v = sort(Float64.(col(t, "n_conv")[ix]))
        stairs!(ax, v, (1:length(v)) ./ length(v); label = "$(key[1]) K$(key[2])")
    end
    Legend(fig[1, 2], ax; labelsize = 8)
    return savefig(fig, out, study, "F3")
end

function f4(t, out, study)
    fig = Figure(size = (640, 420))
    names = string.(col(t, "test"))
    ax = Axis(fig[1, 1]; xticks = (1:length(names), names), xticklabelrotation = pi / 6, ylabel = "share forced to unresolved", title = "label stability")
    barplot!(ax, 1:length(names), Float64.(col(t, "share")))
    return savefig(fig, out, study, "F4")
end

function f5(t, out, study)
    gs = groups(t, "diagnostic")
    g = grid(length(gs); cols = 3)
    fig = Figure(size = (400 * g.cols, 300 * g.rows))
    for (i, (key, ix)) in enumerate(gs)
        ax = Axis(fig[fld(i - 1, g.cols) + 1, mod(i - 1, g.cols) + 1]; title = key[1], xscale = log10, yscale = log10, xlabel = "n", ylabel = "running maximum")
        meth = col(t, "method")
        for m in VM
            sel = [j for j in ix if meth[j] == m]
            isempty(sel) && continue
            lines!(ax, col(t, "n")[sel], max.(Float64.(col(t, "max")[sel]), 1e-20); color = COLORS[m], label = m)
        end
        thr = col(t, "threshold")[ix]
        isempty(thr) || hlines!(ax, [maximum(thr)]; color = :black, linestyle = :dash)
    end
    return savefig(fig, out, study, "F5")
end

function f6(t, out, study)
    fig = Figure(size = (640, 420))
    prec = unique(string.(col(t, "precision")))
    meth = unique(string.(col(t, "method")))
    ax = Axis(fig[1, 1]; xticks = (1:length(meth), meth), yscale = log10, ylabel = "thread ns per orbit-step", title = "cost")
    xs = [findfirst(==(string(m)), meth) for m in col(t, "method")]
    dodge = [findfirst(==(string(p)), prec) for p in col(t, "precision")]
    barplot!(ax, xs, Float64.(col(t, "ns_per_orbit_step")); dodge = dodge, color = dodge)
    return savefig(fig, out, study, "F6")
end

function f7(t, out, study)
    gs = groups(t, "factor")
    isempty(gs) && (gs = [(("none",), Int[])])
    g = grid(length(gs); cols = 2)
    fig = Figure(size = (480 * g.cols, 320 * g.rows))
    for (i, (key, ix)) in enumerate(gs)
        lab = ["$(t.rows[j][2]) $(t.rows[j][3])" for j in ix]
        ax = Axis(fig[fld(i - 1, g.cols) + 1, mod(i - 1, g.cols) + 1]; title = key[1], xticks = (1:length(ix), lab), xticklabelrotation = pi / 4, ylabel = "median difference of lambda(n_H2)")
        isempty(ix) || scatter!(ax, 1:length(ix), Float64.(col(t, "median_dlambda_nH2")[ix]))
    end
    return savefig(fig, out, study, "F7")
end

function draw_figures(res, out)
    mkpath(out)
    f = res.figures
    paths = [f1(f["F1"], out, res.study), f2(f["F2"], out, res.study), f3(f["F3"], out, res.study), f4(f["F4"], out, res.study),
             f5(f["F5"], out, res.study), f6(f["F6"], out, res.study), f7(f["F7"], out, res.study)]
    return paths
end

__init__() = (Analysis.FIGURE_DRAWER[] = draw_figures)

end
