# Experiments

Study inputs built from a frozen preregistration, and the study runner. First scope (decision D49): the orbit panel of a study, read by study/02 and then study/03. Second scope (decision D53): the run plan, the chunked and resumable executor with health checks, the run store (one HDF5 file per chunk), the manifest, the pilot and the command-line entry point `bin/chaoslab.jl` (panel, run, pilot). Third scope (decision D86): the six kinds of chaos indicator (tangent chains, pairs, stretching numbers, windowed frequency, phase advance, weighted Birkhoff) recorded by the same runner, store, manifest and pilot. Fourth scope (decision D94): the panel builder in rounds (random screening points, draws without replacement across rounds, subsets D and D_B) and the run block in rounds (release keys, round log and its guard). Imports Foundation, Systems, Tangent, Lyapunov, Diagnostics and Indicators; nothing imports it. The panel builder is CPU only (class assignment must be bit-identical); the recorder of the runner runs on the CPU, and on the device for the eight kernel-path methods (L1, L2, L3, L4, L4b, L5, L5b, L9; L6 is CPU only) and for the indicator kinds tangent, pair, phase_advance and birkhoff (stretching and frequency are CPU only).

Spec: `workspaces/module-build/stages/01-spec/output/experiments-spec.md` (amendments D50 withdrawn, D51, D52, D53, D54, D55, D56, D57, D58, D59, D86, D87, D88, D94). Store layout: `schema/run-store.md` (version 2). Run block: preregistration Amendment A4, merged with the blocks of A5 and A6 in file order. Reports: `stages/04-verify/output/experiments-verification.md` (first and second scope), `experiments-d86-verification.md` (third scope), `experiments-d94-verification.md` (fourth scope). Oracle scripts: `verify/experiments.jl`, `verify/experiments-d94*.jl` and `experiments-d94-oracle.py` (O23).

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
| `run_indicator.jl` | `[[run.indicator]]` and the arm key `indicators`: validation in two passes, `spec_text`, `spec_hash`, `ind_run_length` (R11, E27) |
| `run_indicator_chain.jl` | chains under the one-replacer rule, observables, reducer tuples, `orbit_bytes`, `group_capacity`, plan-time checks (E23 to E28, R13) |
| `run_indicator_record.jl` | the CPU recording loops of the six kinds, `IndBuf`, the result tables (E24, E25), frame checks (R14) |
| `run_indicator_device.jl` | the series and window kernels, `ind_device_records` (R16) |
| `run_indicator_store.jl` | `indicators/<id>/` datasets, provenance groups, pooled counts, root attributes (R15) |
| `run_indicator_health.jl` | finiteness table E31, `complete` and `failed` defects (R18) |
| `run_indicator_report.jl` | manifest table of the specs, flag counts, per-spec cost lines (R19) |
| `inputs_rounds.jl` | the D94 panel keys (R22): source, `screen_points`, `r_max`, `gate_rounds`, `quota`, `d_fraction`, `db_size`, `newton_grid_side`, `controls` |
| `rounds_streams.jl` | `random_screen_points`, `draw_pool_rounds`, `db_allocation`, `db_subset`: streams and draws of E32 to E38 |
| `rounds_panel.jl` | round cells, subset D, draw checks (R24), pool sidecar, manifest keys of the amendment |
| `round_log.jl` | `read_round_log`, `round_state`, `check_round_log`: the log is read, never written (R26, I37) |
| `run_rounds.jl` | plan with rounds: released rounds and arms (E39), cells `k<i>:<class>:r<round>` (E41), D_B, ref digest, store check of the ref chunks |

## Determinism

Same inputs, same machine and Julia version: byte-identical panel file at any thread count (D52). The manifest records the Julia version and the StableRNGs version; host, time and thread data go only to the sidecar.

Runner: the stored datasets of a work item depend on its inputs and the code version only (I12): bit-identical for any thread count, launch capacity, chunk grouping and run order on the CPU. `auto` for Float64 follows the per-method record on the recording chain (D54: `[<device>."Float64".recorder.<method>]` in `data/gpu-benchmark.toml`, written by `verify/experiments-device.jl`); a missing entry, another thread count or `gpu_faster` false give the CPU. A device chunk records `backend` gpu and agrees with the CPU loop within backend_agreement (bit for bit on the KernelAbstractions CPU backend); a failed device chunk is rerun on the CPU. A resumed campaign equals an uninterrupted one in every dataset and attribute except times, host fields and session lists (I14). No random numbers are used by the runner; pilot starts come from a stream keyed by the pilot seed label.

## Indicators in the runner (D86)

Spec: `experiments-spec.md`, "Amendment D86" (R11 to R21, E21 to E31), with D87 (eleven accepted deviations) and D88 (pooled counts carry their range, flags, E31). Store layout: `schema/run-store.md`, "Added by Amendment D86". What the Indicators layer asked of its caller (chain forms, frame checks, pair setup, the stretching workspace, the series feed of E1) is implemented here:

- Chains (`run_indicator.jl`, `run_indicator_chain.jl`): the module sets the forms (every member but the last non-replacing, SALI only last, D71, D72) and records `chain`, `members`, `spec` and a `provenance/<member>` group per member in every chunk file and the manifest; `check_frame` runs on every frame before the first chunk (a refusal writes nothing), the block frame of a pair per pair (a failure fails its chunk).
- Pair (`run_indicator_record.jl`): `pair_start` on the base map, block frame `(xi_0, 0)`, `(0, xi_0)`, `xi_0` seeded (`initial_deviations(run_seed, pair id)`) or the direction of the arm; `Dx`, `Dx_r`, `mode`, `xi0` are orbit attributes; a separation outside `[Dx/2, 2 Dx]` fails the chunk. The wrapped state and `Dx = 1e-12` stay the study's choice; the pair on Y5 overflows near step 770 (verification F2), so the preregistration keeps its horizon below about 700 or leaves Y5 out.
- Stretching (`run_indicator_record.jl`, `run_indicator_store.jl`): one `UInt32` workspace per orbit, zeroed at step 0; D at every grid step; scalars, `range`, counts for the named subset, pooled exact counts per chunk (`/pooled/<id>`, only for orbits that share one range).
- Frequency (`run_indicator_record.jl`): the orbit pass copies the observable into a per-orbit series buffer (32 W_max bytes in Float64), the pure `frequency_map` runs on it on the CPU, only the table and flags are stored; this is naff-fma criterion 16 (D86j). Phase advance and Birkhoff are read once at the run end; window kinds step one handoff at a time with the deviation column reset (it is never read).
- Memory (`run_indicator_chain.jl`): `orbit_bytes` (M + R + S + B) from the reducers, `group_capacity`; the device batch reads the real chain bytes. Flags are data, never read by a health rule (`run_indicator_health.jl`); the module stores no threshold, class or detection time.
- Backend (`run_backend.jl`, `run_indicator_device.jl`): `stretching` and `frequency` are CPU only and record the sixth fallback reason (a reducer in the chain that is CPU-only); `tangent`, `pair`, `phase_advance` and `birkhoff` follow the per-spec record `[<device>."Float64".recorder.<spec id>]` with its `spec_sha256` (`_config/architecture.md`, "Compute"); Double64 falls back to the CPU.

Limits that remain: the device path was measured on the local RTX 5060 only. Per-spec cost records vary between takes (tangent 1.4 to 3.1 times slower than 23 CPU threads, phase advance 0.6 to 1.1), so the study takes its own record and `auto` for phase advance is not stable. The independent reader of the files is HDF5.jl and `h5dump`, not h5py. Open LOW review findings (R6, R7, R10, R11, R12) are listed in the integration log; the command line entry for the new kinds was not run (the campaigns call `run_campaign` and `run_pilot`).

## Rounds (D94)

Spec: `experiments-spec.md`, "Amendment D94" (E32 to E43, R22 to R27, clarifications of stage 05). With the random source the panel file gains per orbit `round`, the screening point id and `in_D`, per round cell the pool size, size, drawn count and status, and a pool sidecar `<out>.pools.toml` (class bits only, no screening value); a panel input without the D94 keys is byte-identical to before (I28). The plan reads the round log (`--round-log`) and releases each_round arms up to R_plan and closure arms only after `close`; before closure only root attributes and ref arm datasets are read (E43). Limits: the count in the last log entry is not checked here (V1, for the test-side Analysis amendment); the real test panel (D94s) waits for study amendment P69; one machine only.
