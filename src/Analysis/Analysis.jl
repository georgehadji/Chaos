"""
    Analysis

Study analysis (decision D60, `workspaces/module-build/stages/01-spec/output/analysis-spec.md`): reads a finished study run (manifest, run
store, panel file, preregistration) and computes the preregistered tables, the decision of each hypothesis, the labels of D36 and D37 and the
numbers behind each figure. Reads only; owns no scientific threshold (they come from a parameter table). Plotting is a package extension.
"""
module Analysis

using ..Experiments: Experiments, read_chunk, read_json, check_prereg, read_run_block, default_store, git_state, package_version
using SHA
using StableRNGs
using TOML

export analyze, read_params, read_store

include("params.jl")
include("stats.jl")
include("h1.jl")
include("sets.jl")
include("h4.jl")
include("labels.jl")
include("store.jl")
include("checks.jl")
include("orbits.jl")
include("figdata.jl")
include("pipeline.jl")
include("tables.jl")
include("p2_rules.jl")
include("p2_agree.jl")
include("p2_stat.jl")
include("p2_boot.jl")
include("p2_sets.jl")
include("p2_store.jl")
include("p2_params.jl")
include("p2_check.jl")
include("p2_orbits.jl")
include("p2_hyp.jl")
include("p2_report.jl")
include("p2_items.jl")
include("p2_hitems.jl")
include("p2_ritems.jl")
include("p2_render.jl")
include("p2_pipeline.jl")
include("p2_dispatch.jl")
include("cli.jl")

end
