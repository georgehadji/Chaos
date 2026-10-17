# Experiments

Study inputs built from a frozen preregistration. First scope (decision D49): the orbit panel of a study, read by study/02 and then study/03. The runner and run store follow in a later amendment of the same slug, before study/03. Imports Foundation, Systems, Tangent and Lyapunov; nothing imports it. CPU only: class assignment must be bit-identical (spec, Compute).

Spec: `workspaces/module-build/stages/01-spec/output/experiments-spec.md` (amendments D50 withdrawn, D51, D52). Report: `stages/04-verify/output/experiments-verification.md`. Oracle script: `verify/experiments.jl`.

## Files

| File | Holds |
|------|-------|
| `Experiments.jl` | module, exports, include order |
| `inputs.jl` | `panel_inputs`: closed-form strings at BigFloat precision, validation before any computation |
| `prereg.jl` | `read_prereg`, `check_prereg`, `prereg_hash`: the SELF hash rule and the machine-readable panel block |
| `geometry.jl` | `grid_centres`, `torus_delta`, `torus_distance` |
| `screening.jl` | `screen_point!`, `screen_point`, `screen_buffer`, `screen_bytes`: allocation-free L5 screening with a compensated log sum |
| `predicates.jl` | `block_median`, `stuck_length`, `pred_T`, `pred_C`, `pred_S` |
| `periodic.jl` | residual, Newton step, `classify_trace`, `orbit_exponent`, `accept_residual` |
| `newton.jl` | one Newton solve with the D51 conditioning test: `error_bound`, `well_conditioned` |
| `orbits.jl` | orbit records, primitive period, canonical representative, deduplication and the separation refusal (E9, D52) |
| `search.jl` | `periodic_set`, `collect_orbits`: fixed points plus the Newton search per K and period |
| `draw.jl` | `cell_stream_seed`, `cell_rng`, `draw_candidates`, `place_X`, `place_R`, `place_R0`, `orbit_seed`, `y5_control` |
| `cells.jl` | screening of one K and the draws of its cells |
| `subsets.jl` | `select_subsets`: nested S_D and S_L with deterministic caps |
| `panel_tables.jl` | tables of the panel file: periodic orbits, controls, orbits, subsets, cells |
| `panel.jl` | `build_panel`: preregistration in, panel file and provenance sidecar out |
| `portrait.jl` | `portrait_content`: portrait data for the plotting script in `verify/` (plotting stays out of the package) |
| `parallel.jl` | work distribution written by item index, so output never depends on the task count |
| `provenance.jl` | `git_state`, `package_version` |

## Determinism

Same inputs, same machine and Julia version: byte-identical panel file at any thread count (D52). The manifest records the Julia version and the StableRNGs version; host, time and thread data go only to the sidecar.
