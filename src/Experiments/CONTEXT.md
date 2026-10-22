# Experiments

Study inputs built from a frozen preregistration, and the study runner. First scope (decision D49): the orbit panel of a study, read by study/02 and then study/03. Second scope (decision D53): the run plan, the chunked and resumable executor with health checks, the run store (one HDF5 file per chunk), the manifest, the pilot and the command-line entry point `bin/chaoslab.jl` (panel, run, pilot). Imports Foundation, Systems, Tangent, Lyapunov and Diagnostics; nothing imports it. The panel builder is CPU only (class assignment must be bit-identical); the recorder of the runner runs on the CPU, and on the device for the eight kernel-path methods (L1, L2, L3, L4, L4b, L5, L5b, L9; L6 is CPU only).

Spec: `workspaces/module-build/stages/01-spec/output/experiments-spec.md` (amendments D50 withdrawn, D51, D52, D53, D54, D55, D56, D57, D58, D59). Store layout: `schema/run-store.md` (version 2). Run block: preregistration Amendment A4, merged with the blocks of A5 and A6 in file order. Report: `stages/04-verify/output/experiments-verification.md`. Oracle script: `verify/experiments.jl`.

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
| `run_block.jl` | `run_block`, `read_run_block`, `block_text`, `block_hash`: the A4 `[run]` block, validated before any computation (D53a) |
| `run_grid.jl` | `record_grid`, `grid_hash`, `deviation_frame`, `code_version`; exact closed-form evaluation (decimals read as decimals) |
| `run_plan.jl` | `plan_run`, `plan_text`, `chunk_text`, `item_text`, `item_hash`: work items, chunks, hashes, step counts (R1, I11) |
| `run_backend.jl` | `resolve_run_backend`, `backend_table`, `cpu_only_method`, `launch_capacity`, `launch_groups` (R3, R4) |
| `run_device.jl` | `device_records`, `compute_group`: the isbits recorder on the KernelAbstractions kernels (CPU backend in tests, CUDA on the device), step-chunked launches, grid records in preallocated per-orbit arrays, fallback per group |
| `run_record.jl` | `OrbitRecord`, the non-raising residual monitor, the per-orbit recording loop on the shared tangent body, `hot_loop_bytes` (R2) |
| `run_store.jl` | `write_chunk` (atomic), `read_chunk`, `chunk_layout`, `read_json` (R6, R10) |
| `run_health.jl` | health checks of a chunk attempt: non-finite values, exceedance fraction, drift, setup check, exact-exponent controls (R5) |
| `run_campaign.jl` | `RunHooks`, session start checks (git tree, thread count), `exit_status` |
| `run_manifest.jl` | the `Session` state and the manifest (R7) |
| `run_resume.jl` | store scan, stale files, temporary files (R8) |
| `run_groups.jl` | launch-group loop, rerun on CPU, stop and pause rules (R3, R5) |
| `run_session.jl` | `run_campaign`, `default_store` |
| `run_pilot.jl` | `run_pilot`, `rerun_check` (R9, O12) |
| `run_cli.jl` | `cli`: argument parsing of `bin/chaoslab.jl` |

## Determinism

Same inputs, same machine and Julia version: byte-identical panel file at any thread count (D52). The manifest records the Julia version and the StableRNGs version; host, time and thread data go only to the sidecar.

Runner: the stored datasets of a work item depend on its inputs and the code version only (I12): bit-identical for any thread count, launch capacity, chunk grouping and run order on the CPU. `auto` for Float64 follows the per-method record on the recording chain (D54: `[<device>."Float64".recorder.<method>]` in `data/gpu-benchmark.toml`, written by `verify/experiments-device.jl`); a missing entry, another thread count or `gpu_faster` false give the CPU. A device chunk records `backend` gpu and agrees with the CPU loop within backend_agreement (bit for bit on the KernelAbstractions CPU backend); a failed device chunk is rerun on the CPU. A resumed campaign equals an uninterrupted one in every dataset and attribute except times, host fields and session lists (I14). No random numbers are used by the runner; pilot starts come from a stream keyed by the pilot seed label.

## Obligations from the Indicators layer

For the module run that wires indicators into studies (M4, P2):

- Run `check_frame` on every supplied initial frame before the stream (D67); the stream does not.
- Build indicator chains under the one-replacer rule (D71, D72): exactly one reducer replaces the columns and it is last, so `FLI(replace = false)` before `SALI`, and `FLI()` only when alone (likewise `MEGNO(replace = false)` and `MEGNO()`). Record the chain, its order and each reducer's options in the run record; the reducers cannot (fli-ofli review R1, R2).
- RLI (rli, D75 to D77): build each pair with `pair_start` on the base map and run it on `PairMap` with the block frame `(xi_0, 0)`, `(0, xi_0)` of the pair id's seeded `xi_0`, after `check_frame(RLI(), W0)`; record Dx, Dx_r, the mode, the direction rule (u = xi_0), the seed, the frame source and the chain (the stream record cannot hold them, D76 item 6). On the pair stream only RLI and MEGNO are meaningful; SALI, GALI and the FLI basis maximum are not.
