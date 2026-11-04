# Study p3-abstention, stage 02 (checkpoint): phase portraits of the test panel from `portrait_content`, plotted as the p1 portraits (verify/experiments.jl, section PORTRAITS).
# Run from the repo root (verify environment, plotting stays out of the package):  julia --project=verify verify/p3-abstention-portraits.jl <panel file> <png>
using Printf, Chaoslab
import CairoMakie as CM
const EX = Chaoslab.Experiments
const TWO_PI = 2 * pi
const CLASS_ORDER = ["T", "R", "C", "S", "X"]
fmt(x::AbstractFloat) = @sprintf("%.4g", x)
fmt(x) = string(x)
panel, png = ARGS
occursin("/data/", replace(abspath(png), "\\" => "/")) && error("refused: $png is inside data/")
pc = EX.portrait_content(panel; n_plot = 1000)
colors = Dict("T" => "#0072B2", "R" => "#E69F00", "C" => "#D55E00", "S" => "#009E73", "X" => "#CC79A7")
fig = CM.Figure(size = (2100, 1500), fontsize = 15)
for (n, c) in enumerate(pc)
    r, col = fldmod1(n, 3)
    ax = CM.Axis(fig[r, col]; title = "k index $(c.k_index), K = $(fmt(c.K)), $(length(c.starts)) orbits",
                 subtitle = (isempty(c.notes) ? "all cells populated" : join(c.notes, "; ")) * "; " * string(count(q -> q.symbol == :cross, c.periodic)) * " hyperbolic periodic points not drawn",
                 subtitlesize = 11, aspect = CM.DataAspect(), limits = (0, TWO_PI, 0, TWO_PI), xlabel = "x", ylabel = "p",
                 xticks = ([0, pi, TWO_PI], ["0", "π", "2π"]), yticks = ([0, pi, TWO_PI], ["0", "π", "2π"]))
    for cls in CLASS_ORDER
        its = [it for it in c.iterates if it.class == cls]
        isempty(its) && continue
        xs = Float64[x[1] for it in its for x in it.xy]; ps = Float64[x[2] for it in its for x in it.xy]
        CM.scatter!(ax, xs, ps; markersize = 1.6, color = (colors[cls], 0.55))
    end
    for cls in CLASS_ORDER
        st = [s for s in c.starts if s.class == cls]
        isempty(st) || CM.scatter!(ax, [s.x for s in st], [s.p for s in st]; marker = :star5, markersize = 9, color = colors[cls], strokecolor = :black, strokewidth = 0.6)
    end
    for (sym, mk) in ((:circle, :circle), (:diamond, :diamond))
        pp = [q for q in c.periodic if q.symbol == sym]
        isempty(pp) && continue
        CM.scatter!(ax, [q.x for q in pp], [q.p for q in pp]; marker = mk, markersize = sym == :circle ? 6 : 7, color = :transparent, strokecolor = (:black, 0.8), strokewidth = 0.8)
    end
end
els = [CM.MarkerElement(color = colors[k], marker = :circle, markersize = 10) for k in CLASS_ORDER]
labs = ["T torus", "R resonant", "C chaotic sea", "S sticky", "X near separatrix"]
push!(els, CM.MarkerElement(color = :black, marker = :star5, markersize = 11)); push!(labs, "orbit start (star)")
push!(els, CM.MarkerElement(color = :transparent, strokecolor = :black, strokewidth = 1, marker = :circle, markersize = 10)); push!(labs, "periodic point, elliptic")
push!(els, CM.MarkerElement(color = :transparent, strokecolor = :black, strokewidth = 1, marker = :diamond, markersize = 10)); push!(labs, "parabolic")
CM.Legend(fig[3, 1:3], els, labs; orientation = :horizontal, nbanks = 2, tellheight = true, framevisible = false)
CM.Label(fig[0, 1:3], "p3-abstention test panel: first 1000 iterates (wrapped) of every panel orbit, all rounds, $(sum(length(c.starts) for c in pc)) orbits", fontsize = 18)
CM.save(png, fig; px_per_unit = 1)
println("portraits written: ", png, "; starts per K ", join((length(c.starts) for c in pc), ", "), "; periodic markers per K ", join((length(c.periodic) for c in pc), ", "))
