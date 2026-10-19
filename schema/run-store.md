# Run Store Schema (HDF5), version 2

Layer 3, canonical. Shared contract between Julia (`Experiments`) and Python (`python/nbody/`). Any change bumps `schema_version`. Version 2 is the refinement by the study-runner amendment D53 (`_config/decision-log.md`); version 1 files are not produced.

## Files

- One HDF5 file per chunk: `<data root>/[study]/chunk-<NNNNN>.h5` (zero-padded chunk ordinal of the plan).
- Chunks are independent and resumable. `complete = true` means the file was written to the end, whatever its `status`. A file is written under a temporary name, closed, read back, then renamed; no file with a final name is incomplete.

## Root attributes (provenance)

| Attribute | Type | Example |
|-----------|------|---------|
| `schema_version` | int | `2` |
| `study` | string | `p1-le-benchmark` |
| `chunk` | int | `17` |
| `chunk_key` | string | canonical text of run type and cell |
| `chunk_hash` | string | SHA-256 of key, sorted item hashes and the values that decide the status |
| `plan_sha256` | string | hash of the run plan |
| `panel_sha256` | string | plain SHA-256 of the panel file |
| `git_sha` | string | 40-hex commit of a clean tree |
| `code_version` | string | identifier of library source plus dependency lock |
| `prereg_sha256` | string | hash of the frozen preregistration |
| `config_sha256` | string | hash of the canonical text of the run block |
| `julia_version` / `python_version` | string | writer runtime |
| `packages` | string | manifest hash or `pip freeze` hash |
| `cpu`, `threads` | string, int | host info |
| `backend_request` | string | `auto`, `gpu` or `cpu` as requested |
| `backend` | string | `gpu` or `cpu` (resolved, never `auto`) |
| `device` | string | GPU name or `cpu` |
| `backend_fallback_reason` | string | empty, or why `auto`/`gpu` ran on CPU |
| `precision` | string | `Float64`, `Double64`, `BigFloat(256)` |
| `method` | string | method id of the run type |
| `horizon`, `tau`, `vectors` | int | horizon in steps, renormalisation interval, deviation vectors k |
| `coordinate_mode` | string | `lifted` or `wrapped` |
| `direction` | string | initial deviation direction, closed-form text |
| `perturbation` | string | L2 only: perturbation size |
| `grid_sha256` | string | hash of the recording grid |
| `status` | string | `ok` or `failed` |
| `attempt` | int | 1 or 2 (rerun on CPU) |
| `setup_check` | string | `ok`, `failed`, `no reference` or `not applicable` (extended precision against Float64, R5 item 1d) |
| `first_attempt_reasons`, `failure_reasons` | string list | empty when none |
| `pilot` | bool | `false` for study files |
| `arms` | string list | arms that contain the chunk |
| `renorm_interval_tau`, `deviation_vectors_k`, `seed`, `chunk_size`, `gpu_chunks`, `steps_per_launch`, `launches_per_chunk` | int | tangent-core provenance fields |
| `orthonormalizer`, `implementation` | string | Gram-Schmidt provenance (P3) |
| `created_utc` | string | ISO 8601, e.g. `2026-10-12T14:03:00Z` |
| `complete` | bool | `true` |

## Groups

```
/orbits/<orbit_id>/            one group per initial condition
    attrs: system, params, ic, class_label, seed,
           item_hash, design_class, flags, steps_recorded
    t                          int64, recorded step counts, length n
    state                      n x 2, binary64; the recorded orbit point (lifted, or wrapped in wrapped mode)
    winding                    n x 2, int64; wrapped mode only
    lyap/<method>              n x k finite-time exponents per method, binary64
    indicators/<name>          n samples per indicator
    diagnostics/<name>         n running maxima, binary64
        attrs: threshold, first_exceedance_step
    cost/<method>              steps, vectors, wall seconds, thread seconds, resident bytes
/grid/                         grid scans only (not used by the study runner)
    attrs: system, params, axes
    <name>                     final values, shape (nx, ny)
```

`flags` lists any of `diagnostics exceeded`, `drift`, `excluded`. `design_class` is the a-priori class (the former `class_label`, which stays for N-body files). A failed chunk file holds the root attributes and the orbit groups that exist.

N-body files use `system = "nbody:<name>"` and the same groups.

## Changes from v1

- `schema_version` is 2.
- New root attributes: `chunk_key`, `chunk_hash`, `plan_sha256`, `panel_sha256`, `code_version`, `backend_request`, `method`, `horizon`, `tau`, `vectors`, `coordinate_mode`, `direction`, `perturbation`, `grid_sha256`, `status`, `attempt`, `setup_check`, `first_attempt_reasons`, `failure_reasons`, `pilot`, `arms`, the tangent-core provenance fields and `orthonormalizer`, `implementation`.
- `config_sha256` is the hash of the canonical text of the run block.
- `complete` means written to the end (atomic write), whatever `status` is; v1 meant only "chunk finished".
- New orbit attributes: `item_hash`, `design_class`, `flags`, `steps_recorded`.
- New orbit datasets: `state`, `winding`; `diagnostics/<name>` holds running maxima with attributes `threshold` and `first_exceedance_step`; `cost/<method>` gains thread seconds and resident bytes.
- `/grid/` is unchanged and unused by the study runner.
