# Study p3-abstention, stage 04: figures F1 to F3 of the Analysis plan, drawn from the data sections the analysis wrote (the p3 path draws no figure, reading P110).
# F1 from "F1 data" of the calibration tables, F2 and F3 from "F2 data" and "F3 data" of the test tables. Nothing is recomputed.
# Run from the repo root (verify environment, plotting stays out of the package):  julia --project=verify verify/p3-abstention-figures.jl
using Printf
import CairoMakie as CM
const OUT = joinpath(dirname(@__DIR__), "workspaces", "study", "stages", "04-analyze", "output")
const CAL = joinpath(OUT, "p3-abstention-calibration-tables.md")
const TEST = joinpath(OUT, "p3-abstention-tables.md")

# rows of the first markdown table under the heading `## <name>`, as vectors of cell strings (header and separator dropped)
function section_rows(path, name)
    lines = readlines(path)
    i = findfirst(==("## " * name), lines)
    i === nothing && error("$path: no section \"$name\"")
    rows = Vector{Vector{String}}()
    for l in lines[i+1:end]
        startswith(l, "## ") && break
        startswith(l, "|") || (isempty(rows) ? continue : break)
        push!(rows, strip.(split(strip(l, '|'), '|')))
    end
    return rows[3:end]
end
frac(s) = (p = split(s, '/'); length(p) == 1 ? parse(Float64, s) : parse(Float64, p[1]) / parse(Float64, p[2]))

f1 = section_rows(CAL, "F1 data")
@assert length(f1) == 382 "F1: expected |F| = 382 rows (Amendment A1), got $(length(f1))"
cr = frac.(getindex.(f1, 3)); w2 = frac.(getindex.(f1, 4)); ts = frac.(getindex.(f1, 5))
tol = getindex.(f1, 7) .== "true"; star = getindex.(f1, 8) .== "true"
@assert count(tol) == 20 && count(star) == 1 "F1: expected |Theta_tol| = 20 and one theta*"
fig = CM.Figure(size = (1000, 430), fontsize = 13)
for (n, (y, lab)) in enumerate(((w2, "FR(W2; θ)"), (ts, "FR(T_s; θ)")))
    ax = CM.Axis(fig[1, n]; xlabel = "C_R(θ, b_P) (thread-ns)", ylabel = lab, title = "F1: calibration landscape over F (|F| = 382)")
    CM.scatter!(ax, cr[.!tol], y[.!tol]; markersize = 5, color = (:gray40, 0.5), label = "θ in F")
    CM.scatter!(ax, cr[tol .& .!star], y[tol .& .!star]; markersize = 7, color = "#0072B2", label = "θ in Θ_tol")
    CM.scatter!(ax, cr[star], y[star]; marker = :star5, markersize = 16, color = "#D55E00", strokecolor = :black, strokewidth = 0.6, label = "θ*")
    n == 2 && CM.axislegend(ax; position = :rt, framevisible = false)
end
CM.save(joinpath(OUT, "p3-abstention-fig-F1.pdf"), fig)

f2 = section_rows(TEST, "F2 data")
@assert length(f2) == 10 "F2: expected 3 rule points and 7 single indicators, got $(length(f2))"
fig = CM.Figure(size = (620, 480), fontsize = 13)
ax = CM.Axis(fig[1, 1]; xlabel = "decided share on E_test", ylabel = "error e among the decided", title = "F2: risk against coverage (matched budgets)")
CM.xlims!(ax, nothing, 1.025)  # room for the labels of the points at decided share 1
for (kind, col, mk) in (("rule", "#D55E00", :star5), ("single", "#0072B2", :circle))
    r = [x for x in f2 if x[1] == kind]
    xs = parse.(Float64, getindex.(r, 3)); ys = parse.(Float64, getindex.(r, 4))
    CM.scatter!(ax, xs, ys; marker = mk, markersize = kind == "rule" ? 15 : 10, color = col, strokecolor = :black, strokewidth = 0.5, label = kind == "rule" ? "rule family R_a" : "single indicator")
    CM.text!(ax, xs, ys; text = getindex.(r, 2), fontsize = 10, offset = (6, 2))
end
CM.axislegend(ax; position = :lt, framevisible = false)
CM.save(joinpath(OUT, "p3-abstention-fig-F2.pdf"), fig)

f3 = section_rows(TEST, "F3 data")
fig = CM.Figure(size = (520, 500), fontsize = 13)
ax = CM.Axis(fig[1, 1]; xlabel = "false positive rate", ylabel = "true positive rate (label failure f)", title = "F3: ROC of s_I and s_P on D", aspect = 1, limits = (0, 1, 0, 1))
CM.lines!(ax, [0, 1], [0, 1]; color = :gray70, linestyle = :dash)
for (curve, col) in (("s_I", "#D55E00"), ("s_P", "#0072B2"))
    r = [x for x in f3 if x[1] == curve]
    auc = parse(Float64, only(x[3] for x in r if x[2] == "AUC"))
    pts = sort([(parse(Float64, x[3]), parse(Float64, x[4])) for x in r if x[2] != "AUC"])
    xs = vcat(0.0, first.(pts)); ys = vcat(0.0, last.(pts))
    CM.lines!(ax, xs, ys; color = col, linewidth = 2, label = @sprintf("%s (AUC %.4f)", curve, auc))
    CM.scatter!(ax, xs, ys; color = col, markersize = 7)
end
CM.axislegend(ax; position = :rb, framevisible = false)
CM.save(joinpath(OUT, "p3-abstention-fig-F3.pdf"), fig)
println("figures written: ", join(("p3-abstention-fig-F$k.pdf" for k in 1:3), ", "), " in ", OUT)
