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

# n axes: ticks at whole powers of ten (every second one), labelled 10^k
const LOGN = (xticks = (lo, hi) -> 10.0 .^ (2 * fld(floor(Int, log10(lo)), 2):2:ceil(Int, log10(hi))),
              xtickformat = vs -> [rich("10", superscript(string(round(Int, log10(v))))) for v in vs])
# F3: colour by design class, line style by K index
const CLASS_COLORS = Dict("C" => :red, "R" => :blue, "S" => :green, "T" => :orange, "X" => :purple)
const K_STYLES = Dict(1 => (:dash, :loose), 2 => :solid, 3 => :dash, 4 => :dot, 5 => :dashdot, 6 => :dashdotdot)
const PREC_COLORS = [:steelblue, :darkorange, :seagreen, :purple]   # F6: by order of first appearance of the precision label

method_legend(fig, pos; extra = Pair[]) = Legend(fig[pos...], vcat([LineElement(color = COLORS[m]) for m in vcat(VM, "L2")], first.(extra)),
                                                 vcat(vcat(VM, "L2"), last.(extra)), "method"; labelsize = 12, framevisible = false)

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
        ax = Axis(fig[fld(i - 1, g.cols) + 1, mod(i - 1, g.cols) + 1]; title = "class $(key[1]), K index $(key[2])", xscale = log10, LOGN..., xlabel = "n", ylabel = "lambda_1(n)")
        n = col(t, "n")[ix]
        lo = [minimum(t.rows[j][findfirst(==(m), t.columns)] for m in VM) for j in ix]
        hi = [maximum(t.rows[j][findfirst(==(m), t.columns)] for m in VM) for j in ix]
        band!(ax, n, lo, hi; color = (:gray, 0.3))
        for m in vcat(VM, "L2")
            lines!(ax, n, col(t, m)[ix]; color = COLORS[m], label = m)
        end
    end
    method_legend(fig, (1:g.rows, g.cols + 1); extra = [PolyElement(color = (:gray, 0.3)) => "Delta_V band (L2 excluded)"])
    return savefig(fig, out, study, "F1")
end

function f2(t, out, study)
    gs = groups(t, "class", "k_index")
    g = grid(length(gs))
    fig = Figure(size = (320 * g.cols, 260 * g.rows))
    for (i, (key, ix)) in enumerate(gs)
        title = key[1] == "W" ? "set W of H2, pooled" : "class $(key[1]), K index $(key[2])"
        ax = Axis(fig[fld(i - 1, g.cols) + 1, mod(i - 1, g.cols) + 1]; title = title, xscale = log10, LOGN..., xlabel = "n", ylabel = "R(n)")
        n = col(t, "n")[ix]
        band!(ax, n, col(t, "lo")[ix], col(t, "hi")[ix]; color = (:blue, 0.25))
        lines!(ax, n, col(t, "median")[ix]; color = :blue)
    end
    Legend(fig[1:g.rows, g.cols + 1], [LineElement(color = :blue), PolyElement(color = (:blue, 0.25))],
           ["median over the orbits", "95 percent bootstrap band"]; labelsize = 12, framevisible = false)
    return savefig(fig, out, study, "F2")
end

function f3(t, out, study)
    gs = groups(t, "class", "k_index")
    fig = Figure(size = (820, 480))
    ax = Axis(fig[1, 1]; xscale = log10, LOGN..., xlabel = "n_conv", ylabel = "fraction of orbits", title = "convergence time")
    for (key, ix) in gs
        v = sort(Float64.(col(t, "n_conv")[ix]))
        stairs!(ax, v, (1:length(v)) ./ length(v); color = CLASS_COLORS[key[1]], linestyle = K_STYLES[key[2]])
    end
    cls = sort(unique(first.(first.(gs)))); ks = sort(unique(last.(first.(gs))))
    Legend(fig[1, 2], [[LineElement(color = CLASS_COLORS[c]) for c in cls], [LineElement(color = :black, linestyle = K_STYLES[k]) for k in ks]],
           [cls, ["K index $k" for k in ks]], ["design class", "K index"]; labelsize = 12, framevisible = false)
    return savefig(fig, out, study, "F3")
end

function f4(t, out, study)
    fig = Figure(size = (640, 420))
    names = string.(col(t, "test"))
    ax = Axis(fig[1, 1]; xticks = (1:length(names), names), xticklabelrotation = pi / 6, ylabel = "share of orbits", title = "label stability")
    summary = startswith.(names, "unresolved")   # the last two rows: unresolved share of all orbits, not forced by one test
    barplot!(ax, 1:length(names), Float64.(col(t, "share")); color = [s ? :gray : :steelblue for s in summary])
    Legend(fig[2, 1], [PolyElement(color = :steelblue), PolyElement(color = :gray)],
           ["forced to unresolved by the test (of the orbits it applies to)", "unresolved, all orbits"]; orientation = :horizontal, labelsize = 11, framevisible = false)
    return savefig(fig, out, study, "F4")
end

function f5(t, out, study)
    gs = groups(t, "diagnostic")
    g = grid(length(gs); cols = 3)
    fig = Figure(size = (400 * g.cols, 300 * g.rows))
    for (i, (key, ix)) in enumerate(gs)
        ax = Axis(fig[fld(i - 1, g.cols) + 1, mod(i - 1, g.cols) + 1]; title = key[1], xscale = log10, LOGN..., yscale = log10, xlabel = "n",
                  ylabel = "largest running maximum over the orbits")
        meth = col(t, "method")
        for m in VM
            sel = [j for j in ix if meth[j] == m]
            isempty(sel) && continue
            lines!(ax, col(t, "n")[sel], max.(Float64.(col(t, "max")[sel]), 1e-20); color = COLORS[m])
        end
        # one dashed line per distinct threshold: black when several methods share it, else in the method's colour
        thr = Dict{Float64, Vector{String}}()
        for j in ix
            push!(get!(thr, round(Float64(col(t, "threshold")[j]); sigdigits = 6), String[]), string(meth[j]))   # equal up to rounding
        end
        for (v, ms) in thr
            hlines!(ax, [v]; color = length(unique(ms)) > 1 ? :black : COLORS[first(ms)], linestyle = :dash, linewidth = 1)
        end
    end
    method_legend(fig, (1:g.rows, g.cols + 1); extra = [LineElement(color = :black, linestyle = :dash) => "threshold, median over the orbits
(black: shared by several methods)"])
    return savefig(fig, out, study, "F5")
end

function f6(t, out, study)
    fig = Figure(size = (640, 420))
    prec = unique(string.(col(t, "precision")))
    meth = unique(string.(col(t, "method")))
    ax = Axis(fig[1, 1]; xticks = (1:length(meth), meth), yscale = log10, ylabel = "thread ns per orbit-step", title = "cost")
    xs = [findfirst(==(string(m)), meth) for m in col(t, "method")]
    dodge = [findfirst(==(string(p)), prec) for p in col(t, "precision")]
    barplot!(ax, xs, Float64.(col(t, "ns_per_orbit_step")); dodge = dodge, color = [PREC_COLORS[d] for d in dodge])
    Legend(fig[1, 2], [PolyElement(color = PREC_COLORS[i]) for i in eachindex(prec)], prec, "precision"; framevisible = false)
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
