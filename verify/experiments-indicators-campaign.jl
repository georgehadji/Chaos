# Child process of verify/experiments-indicators.jl: one uninterrupted reduced campaign (D86y) in its own process, so that the wall time and the peak
# host memory (peak working set, Sys.maxrss) belong to the campaign alone. Usage:
#   julia --project=verify -t 23 verify/experiments-indicators-campaign.jl <store dir> <manifest> <tag> [capacity]
# The one injected chunk failure (a NaN in a series column of the second orbit of the first tangent chunk, on attempts 1 and 2) is the same as in the main script.
include(joinpath(@__DIR__, "experiments-indicators-prep.jl"))
store, manifest, tag = ARGS[1], ARGS[2], ARGS[3]
cap = length(ARGS) >= 4 ? parse(Int, ARGS[4]) : nothing
p = (prereg = joinpath(SCR, "prereg-reduced.md"), panel = joinpath(SCR, "panel-reduced.toml"))
isfile(p.prereg) && isfile(p.panel) || prepare()
plan = EX.plan_run(p.prereg, p.panel)
ord_fail = first(c.ordinal for c in plan.chunks if c.run_type.method == "tg" && c.run_type.precision == "Float64" && c.run_type.coordinate_mode == "wrapped")
rm(store; recursive = true, force = true); rm(manifest; force = true)
rss0 = Sys.maxrss()
t = @elapsed r = EX.run_campaign(p.prereg, p.panel; store = store, manifest = manifest, threads = NT, thread_override = NT != 23, allow_dirty = true,
                                 capacity = cap, hooks = EX.RunHooks(mutate = nan_fault(ord_fail, "tg"; attempts = (1, 2))))
println("CAMPAIGN $tag status=$(r.status) chunks=$(length(plan.chunks)) ran=$(length(r.ran)) failed=$(length(r.failed)) threads=$NT wall_s=$(round(t; digits = 2)) " *
        "maxrss_before_bytes=$rss0 maxrss_bytes=$(Sys.maxrss()) injected_ordinal=$ord_fail capacity=$(cap === nothing ? "memory" : cap) patch=$(get(ENV, "EXI_PATCH", "1"))")
