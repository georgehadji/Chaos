"""
    Experiments

Study panel builder (experiments, decision D49): reads the machine-readable panel block of a frozen preregistration,
screens a grid per K, finds the periodic orbits, draws the panel, writes the panel file with its manifest and returns the
content of the phase portraits. The runner, the run store and the command-line entry point are not part of this version.
"""
module Experiments

using ..Foundation
using ..Systems
using ..Tangent
using ..Lyapunov
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

end
