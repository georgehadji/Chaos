# experiments suite (stage 02, RED). Spec: workspaces/module-build/stages/01-spec/output/experiments-spec.md
# API under test: submodule Chaoslab.Experiments (names listed in workspaces/module-build/stages/02-tests/output/experiments-test-plan.md).
# Nothing here calls Experiments at load time: a missing module errors the testsets that need it and nothing else.
include(joinpath(@__DIR__, "ex_helpers.jl"))
include(joinpath(@__DIR__, "ex_inputs.jl"))
include(joinpath(@__DIR__, "ex_geometry.jl"))
include(joinpath(@__DIR__, "ex_periodic.jl"))
include(joinpath(@__DIR__, "ex_screening.jl"))
include(joinpath(@__DIR__, "ex_predicates.jl"))
include(joinpath(@__DIR__, "ex_draw.jl"))
include(joinpath(@__DIR__, "ex_subsets.jl"))
include(joinpath(@__DIR__, "ex_panel.jl"))
include(joinpath(@__DIR__, "ex_determinism.jl"))
# Amendment D53 (study runner and run store): fixtures and reference code in ex_run_helpers.jl, then one file per group of criteria.
include(joinpath(@__DIR__, "ex_run_helpers.jl"))
include(joinpath(@__DIR__, "ex_run_inputs.jl"))
include(joinpath(@__DIR__, "ex_plan.jl"))
include(joinpath(@__DIR__, "ex_chunks.jl"))
include(joinpath(@__DIR__, "ex_store.jl"))
include(joinpath(@__DIR__, "ex_resume.jl"))
include(joinpath(@__DIR__, "ex_health.jl"))
include(joinpath(@__DIR__, "ex_manifest.jl"))
include(joinpath(@__DIR__, "ex_pilot.jl"))
include(joinpath(@__DIR__, "ex_device.jl"))
include(joinpath(@__DIR__, "ex_record.jl"))
include(joinpath(@__DIR__, "ex_review.jl"))   # Amendment D53 review fixes
include(joinpath(@__DIR__, "ex_a5.jl"))   # Amendment D56 and preregistration A5
include(joinpath(@__DIR__, "ex_a6.jl"))   # Amendment D57 and preregistration A6
include(joinpath(@__DIR__, "ex_prior.jl"))   # Amendment D58: --prior-record
include(joinpath(@__DIR__, "ex_backup.jl"))   # Amendment D59: manifest backups
# Amendment D86 (indicators in the runner and the store): helpers, then one file per group of criteria
include(joinpath(@__DIR__, "ex_d86_helpers.jl"))
include(joinpath(@__DIR__, "ex_d86_grammar.jl"))   # D86a to D86d
include(joinpath(@__DIR__, "ex_d86_chains.jl"))    # D86e to D86h
include(joinpath(@__DIR__, "ex_d86_windows.jl"))   # D86i to D86m
include(joinpath(@__DIR__, "ex_d86_run.jl"))       # D86n to D86r
include(joinpath(@__DIR__, "ex_d86_store.jl"))     # D86s to D86v
include(joinpath(@__DIR__, "ex_d86_guards.jl"))    # D86w, D86x
include(joinpath(@__DIR__, "ex_d86_review.jl"))   # stage 04 review fixes
# Amendment D94 (panel builder in rounds and the study runner with the round log): helpers, then one file per group of criteria
include(joinpath(@__DIR__, "ex_d94_helpers.jl"))
include(joinpath(@__DIR__, "ex_d94_grammar.jl"))   # D94a, D94b, D94k, D94u
include(joinpath(@__DIR__, "ex_d94_streams.jl"))   # D94c, D94f, D94h, D94i (function level)
include(joinpath(@__DIR__, "ex_d94_panel.jl"))     # D94c to D94j, D94p, D94t (builder level)
include(joinpath(@__DIR__, "ex_d94_plan.jl"))      # D94i, D94l, D94m
include(joinpath(@__DIR__, "ex_d94_run.jl"))       # D94n, D94o
# Amendment D96 (indicator specs held only by closure arms)
include(joinpath(@__DIR__, "ex_d96.jl"))   # D96a, D96b
# Amendment D125 (higher-dimensional maps and flows): helpers, reference, golden records, then the criteria D125a to D125u
include(joinpath(@__DIR__, "ex_d125_helpers.jl"))
include(joinpath(@__DIR__, "ex_d125_ref.jl"))
include(joinpath(@__DIR__, "ex_d125_golden_fn.jl"))
include(joinpath(@__DIR__, "ex_d125_a.jl"))        # D125a (golden), D125u
include(joinpath(@__DIR__, "ex_d125_grammar.jl"))  # D125b, D125i
include(joinpath(@__DIR__, "ex_d125_panel.jl"))    # D125d to D125h
include(joinpath(@__DIR__, "ex_d125_items.jl"))    # D125c, D125j to D125m
include(joinpath(@__DIR__, "ex_d125_est.jl"))      # D125n, D125o
include(joinpath(@__DIR__, "ex_d125_store.jl"))    # D125p to D125t
# Amendment D126 (finite-time forms of the E58 pairing threshold and of the control tolerance; lower-layer, passes before stage 03)
include(joinpath(@__DIR__, "ex_d126_bounds.jl"))  # D126a, D126b
