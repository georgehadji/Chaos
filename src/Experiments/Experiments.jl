"""
    Experiments

Study panel builder (experiments, decision D49): reads the machine-readable panel block of a frozen preregistration,
screens a grid per K, finds the periodic orbits, draws the panel, writes the panel file with its manifest and returns the
content of the phase portraits. Second scope (decision D53): the study runner and the run store (run plan, chunked and
resumable executor with health checks, one HDF5 file per chunk, manifest, pilot) and the entry point of the command line.
"""
module Experiments

using ..Foundation
using ..Systems
using ..Tangent
using ..Lyapunov
using ..Diagnostics
using StaticArrays
using LinearAlgebra
using StableRNGs
using SHA
using TOML

export panel_inputs, read_prereg, check_prereg, prereg_hash, build_panel, portrait_content
export error_bound, well_conditioned, periodic_set, collect_orbits, classify_trace, orbit_exponent, newton_correction, accept_residual, periodic_distance
export grid_centres, torus_delta, torus_distance
export screen_buffer, screen_point, screen_point!, screen_bytes
export block_median, stuck_length, pred_T, pred_C, pred_S
export cell_stream_seed, cell_rng, draw_candidates, place_X, place_R, place_R0, orbit_seed, y5_control
export select_subsets, git_state, package_version
export run_block, read_run_block, block_text, block_hash, plan_run, plan_text, chunk_text, item_text, item_hash, record_grid, grid_hash
export deviation_frame, launch_capacity, launch_groups, resolve_run_backend, backend_table, cpu_only_method, code_version
export run_campaign, RunHooks, read_chunk, write_chunk, chunk_layout, read_json, run_pilot, rerun_check, exit_status, cli
export default_store, hot_loop_bytes

include("parallel.jl")
include("geometry.jl")
include("inputs.jl")
include("prereg.jl")
include("periodic.jl")
include("newton.jl")
include("orbits.jl")
include("search.jl")
include("screening.jl")
include("predicates.jl")
include("draw.jl")
include("subsets.jl")
include("provenance.jl")
include("cells.jl")
include("panel_tables.jl")
include("panel.jl")
include("portrait.jl")
include("run_block.jl")
include("run_grid.jl")
include("run_plan.jl")
include("run_backend.jl")
include("run_record.jl")
include("run_store.jl")
include("run_health.jl")
include("run_campaign.jl")
include("run_manifest.jl")
include("run_resume.jl")
include("run_groups.jl")
include("run_session.jl")
include("run_pilot.jl")
include("run_cli.jl")

end
