# analysis suite (stage 02, RED). Spec: workspaces/module-build/stages/01-spec/output/analysis-spec.md
# API under test: submodule Chaoslab.Analysis (names listed in workspaces/module-build/stages/02-tests/output/analysis-test-plan.md).
# Nothing here calls Analysis at load time: a missing module errors the testsets that need it and nothing else.
include(joinpath(@__DIR__, "an_helpers.jl"))
include(joinpath(@__DIR__, "an_fixture.jl"))      # fixture self-test: passes before the module exists
include(joinpath(@__DIR__, "an_params.jl"))
include(joinpath(@__DIR__, "an_h1.jl"))
include(joinpath(@__DIR__, "an_sets.jl"))
include(joinpath(@__DIR__, "an_bootstrap.jl"))
include(joinpath(@__DIR__, "an_h4.jl"))
include(joinpath(@__DIR__, "an_labels.jl"))
include(joinpath(@__DIR__, "an_h5.jl"))
include(joinpath(@__DIR__, "an_pipeline.jl"))
include(joinpath(@__DIR__, "an_cli.jl"))
include(joinpath(@__DIR__, "an_breakdown.jl"))     # regression of stage 04: breakdown_tau items on the real store
include(joinpath(@__DIR__, "an_review.jl"))        # regressions of the stage 04 review

# D89 stage 02 (Amendment D89, study p2-indicators). The module has no D89 functionality yet: the files below are RED, except the fixture self-test and the guards of the
# P1 path (an_p2_dispatch.jl, first testset). Helpers first, then one file per theme.
include(joinpath(@__DIR__, "an_p2_ref.jl"))          # independent reference implementation of the decision layer
include(joinpath(@__DIR__, "an_p2_helpers.jl"))      # parameter table, preregistration, run blocks
include(joinpath(@__DIR__, "an_p2_study.jl"))        # designed orbits and the writer of the chunk records
include(joinpath(@__DIR__, "an_p2_build.jl"))        # assembly of the synthetic study, editing helpers
include(joinpath(@__DIR__, "an_p2_refstore.jl"))     # expected tables computed from the designed study
include(joinpath(@__DIR__, "an_p2_fixture.jl"))      # fixture self-test: passes before the module exists
include(joinpath(@__DIR__, "an_p2_dispatch.jl"))     # D89a, D89t
include(joinpath(@__DIR__, "an_p2_params.jl"))       # D89b
include(joinpath(@__DIR__, "an_p2_verdicts.jl"))     # D89e, D89f, D89g, D89y
include(joinpath(@__DIR__, "an_p2_detection.jl"))    # D89h
include(joinpath(@__DIR__, "an_p2_agreement.jl"))    # D89i, D89j
include(joinpath(@__DIR__, "an_p2_sets.jl"))         # D89k, D89l
include(joinpath(@__DIR__, "an_p2_stats.jl"))        # D89m, D89n (statistics and decisions)
include(joinpath(@__DIR__, "an_p2_bootstrap.jl"))    # D89o
include(joinpath(@__DIR__, "an_p2_refusals.jl"))     # D89c, D89d, D89r after refusals
include(joinpath(@__DIR__, "an_p2_pipeline.jl"))     # sets, verdicts, detection, decisions, extension, R2, gate, AUC end to end
include(joinpath(@__DIR__, "an_p2_variant.jl"))      # drift orbit, sticky orbit 5, gate violation and second versions, NOT ROBUST
include(joinpath(@__DIR__, "an_p2_reported.jl"))     # T11, T15, T16, T17
include(joinpath(@__DIR__, "an_p2_output.jl"))       # D89p, D89q, D89r, D89s, D89u
include(joinpath(@__DIR__, "an_p2_review.jl"))       # review findings of stage 04 (K value pooling, infinite percentiles, stripped attributes, text against value, ...)
include(joinpath(@__DIR__, "an_p2_figures.jl"))      # D91 (module side of the figure drawers: a fake FIGURE_DRAWER, status, command lines, layering); the pictures are verify/analysis-p2-figures.jl

# D93 stage 02 (Amendment D93, study p3-abstention). The module has no D93 functionality yet: the files below are RED, except the fixture self-test and the guards of the P1 and P2 paths
# (an_p3_dispatch.jl). Helpers first, then one file per theme.
include(joinpath(@__DIR__, "an_p3_ref.jl"))          # independent reference implementation of the calibration layer (exact rationals, generators, the HB simulation)
include(joinpath(@__DIR__, "an_p3_helpers.jl"))      # parameter table of 71 keys, preregistration text, designed study, P1 reference with N_1 = 80000, prescribed tables
include(joinpath(@__DIR__, "an_p3_fixture.jl"))      # fixture self-test: passes before the module exists
include(joinpath(@__DIR__, "an_p3_rule.jl"))         # D93g, D93h, D93i
include(joinpath(@__DIR__, "an_p3_sets.jl"))         # D93d (the rule), D93e
include(joinpath(@__DIR__, "an_p3_select.jl"))       # D93j (O18 a to e, the order of C6)
include(joinpath(@__DIR__, "an_p3_pp.jl"))           # D93k, D93l, D93m
include(joinpath(@__DIR__, "an_p3_params.jl"))       # D93b
include(joinpath(@__DIR__, "an_p3_dispatch.jl"))     # D93a, D93q
include(joinpath(@__DIR__, "an_p3_refusals.jl"))     # D93c, D93p after refusals
include(joinpath(@__DIR__, "an_p3_pipeline.jl"))     # D93d to D93f, D93j to D93l, D93n to D93p end to end
include(joinpath(@__DIR__, "an_p3_gate.jl"))         # D93c P68: a gate orbit outside Cal with a missing verdict item (stage 03 rework, owner decision 3)

# D95 stage 02 (Amendment D95, study p3-abstention, test side). The module has no D95 functionality yet: the files below are RED, except the fixture self-test (an_p3t_fixture.jl), the
# McNemar and stream guards of an_p3t_stats.jl and the layering guards of an_p3t_dispatch.jl. Reference and helpers first, then one file per theme.
include(joinpath(@__DIR__, "an_p3t_ref.jl"))         # independent reference of the test side (round log entry, count rule, type 7 percentile, AUC, ROC, strata, HA, HB, HC, power at achieved sizes)
include(joinpath(@__DIR__, "an_p3t_helpers.jl"))     # toy campaign: 31 orbits, 3 rounds, 50 chunk files; panel with rounds; test-side table, record, preregistration with Amendment A1; instances
include(joinpath(@__DIR__, "an_p3t_fixture.jl"))     # fixture self-test: passes before the module exists
include(joinpath(@__DIR__, "an_p3t_stats.jl"))       # D95h D95i D95k on prescribed tables, D95m, O21 (a) to (g)
include(joinpath(@__DIR__, "an_p3t_params.jl"))      # D95b
include(joinpath(@__DIR__, "an_p3t_log.jl"))         # D95d, D95e
include(joinpath(@__DIR__, "an_p3t_refusals.jl"))    # D95c
include(joinpath(@__DIR__, "an_p3t_dispatch.jl"))    # D95a, D95p
include(joinpath(@__DIR__, "an_p3t_pipeline.jl"))    # D95f to D95l, D95n, D95o end to end
