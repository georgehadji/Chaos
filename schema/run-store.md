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
| `dirty` | bool | the working tree of the writing session was not clean (an override of the run, recorded in the manifest) |
| `session` | int | 1-based index of the writing session in the manifest's session list |
| `control_stop` | bool | `true` only for a file written when a control exceeded its bound and the campaign stopped; such a file is not a result, a resume refuses it unless the retry is asked for |
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
    lyap/<method>              n x k finite-time exponents per method, binary64; absent for the items of an indicator spec
    indicators/<id>/           datasets of an indicator spec (Amendment D86, below)
    diagnostics/<name>         n running maxima, binary64
        attrs: threshold, first_exceedance_step
    cost/<method>              steps, vectors, wall seconds, thread seconds, resident bytes
/provenance/<member>/         indicator specs only: one group per chain member (Amendment D86, below)
/pooled/<id>/                  indicator specs only: stretching counts pooled over the orbits of the chunk (Amendment D86, below)
/grid/                         grid scans only (not used by the study runner)
    attrs: system, params, axes
    <name>                     final values, shape (nx, ny)
```

`flags` lists any of `diagnostics exceeded`, `drift`, `excluded`, `breakdown_tau`. `design_class` is the a-priori class (the former `class_label`, which stays for N-body files). A failed chunk file holds the root attributes and the orbit groups that exist.

Axis order. The tables above write the Julia shape `(n, k)`. HDF5.jl stores a Julia array column-major with the dimensions reversed, so the file extents are `(k, n)` (`h5ls` and h5py show `(k, n)` for `lyap/<method>` and `(2, n)` for `state`); a reader in another language indexes the recorded step first in Julia and last in the file. Julia readers through `read_chunk` see the shapes of the tables.

Failed chunk files. Rows after the failure point are blanked, not removed: the recorded step count `t` is 0 and the values are `NaN` (state, exponents, diagnostics) from the first row that was not computed to the end of the array, so every dataset keeps its full length. `status = "failed"` and `failure_reasons` say why; an orbit that failed twice has the flag `excluded`. A file written by a control-bound stop (`control_stop = true`) has no `excluded` flag.

N-body files use `system = "nbody:<name>"` and the same groups.

## Changes from v1

- `schema_version` is 2.
- New root attributes: `chunk_key`, `chunk_hash`, `plan_sha256`, `panel_sha256`, `code_version`, `backend_request`, `method`, `horizon`, `tau`, `vectors`, `coordinate_mode`, `direction`, `perturbation`, `grid_sha256`, `status`, `attempt`, `setup_check`, `first_attempt_reasons`, `failure_reasons`, `pilot`, `arms`, the tangent-core provenance fields and `orthonormalizer`, `implementation`.
- `config_sha256` is the hash of the canonical text of the run block.
- `complete` means written to the end (atomic write), whatever `status` is; v1 meant only "chunk finished".
- New orbit attributes: `item_hash`, `design_class`, `flags`, `steps_recorded`.
- New orbit datasets: `state`, `winding`; `diagnostics/<name>` holds running maxima with attributes `threshold` and `first_exceedance_step`; `cost/<method>` gains thread seconds and resident bytes.
- `/grid/` is unchanged and unused by the study runner.

## Review fixes before any study data (version stays 2)

Added by the D53 review fixes, while no study file exists: root attributes `dirty`, `session`, `control_stop`; the statements on axis order and on failed files above. The write now also compares the root attributes on the read-back and flushes the temporary file to the disk before the rename.

- Added by Amendment D56 (version stays 2): orbit flag `breakdown_tau` with orbit attributes `breakdown_step` (int) and `breakdown_kind` (string) for an orthonormalisation breakdown in an arm marked `stop_exempt`; orbit attributes `control_scope` (`"bound"` or `"info"`) and `control_deviation` (float) for the controls of an `h1_check` arm when the run uses `control_scope = "exact_start"`.

## Added by Amendment D86 (version stays 2)

Recording of the six indicator kinds by the study runner (decisions D86, D87, D88; spec `experiments-spec.md`, R15). Version decision (P34): the store stays at version 2. The change is additive, as D56 was: a file of a catalog method is written and read as before, and a resume must equal an uninterrupted run (I14). Two rules change for readers: `lyap/` is optional (an orbit of an indicator spec has none) and the root attribute `kind` marks the files of the new kinds (absent in the files of the methods). `<id>` is the spec id; in a file of a spec the root `method` is the spec id, `direction` is absent when the item has none and `perturbation` stays absent. Shapes are Julia shapes (see Axis order).

Root attributes of a chunk file of a spec (taken from the first item of the chunk):

| Attribute | Type | Content |
|-----------|------|---------|
| `kind` | string | `tangent`, `pair`, `stretching`, `frequency`, `phase_advance` or `birkhoff` |
| `spec`, `spec_sha256` | string | canonical text of the spec and its hash (the same `spec_sha256` as in the GPU record, `_config/architecture.md`, "Compute") |
| `chain`, `members` | string list | the chain record (the form of each member) and the member names in order; a window kind holds its own name |
| `observable_names` | string list | `birkhoff` only: one name per observable |
| `windows`, `run_length` | int list, int | window kinds only |
| `frame_source`, `direction_rule`, `separation` | string, string, float | `pair` only: the frame source, the direction rule, the nominal separation |

Group `provenance/<member>`: one attribute per published key of that member's reducer (kept apart so that no key is overwritten); a value HDF5 does not write natively (Double64, BigFloat, a tuple) is stored as its exact text. Window kinds hold one group named by the kind.

Datasets under `orbits/<n>/indicators/<id>/`:

| Dataset | Type, shape | Attributes | Kinds |
|---------|-------------|------------|-------|
| `series` | binary64, n x c | `columns` | `tangent`, `pair`, `stretching` (at the recording grid) |
| `final` | binary64, rows x columns | `columns`, `rows` (window or window-and-observable labels) | window kinds |
| `flags` | int64, rows | `bits` (the bit names, in order; the name at 1-based position j has the value 2^(j-1), so the first name is 1, the second 2, and a word is the sum of the values of its set bits; stored words fit this, for example a `mirror1` word 128 and a `degenerate` + `refined1` + `refined2` word 1538 = 2 + 512 + 1024; fixed by P34 of Amendment D89, no change of the format) | window kinds; the tangent family and stretching have no flag word |
| `scalars` | binary64 | `columns`; for `stretching` also `range` (two binary64, lo then hi: the range actually used, derived or given) | `stretching`, `phase_advance`, `birkhoff` |
| `counts` | uint32, (bins + 2) x 2 | none | `stretching`, only for the orbits of the spec's `counts` subset |

`t`, `state`, `winding`, `diagnostics/<name>` and `cost/<id>` are as for any method; `cost/<id>` of `frequency` has a sixth entry (transform seconds). A `pair` orbit has the attributes `Dx`, `Dx_r` (the realized separation), `mode`, `xi0` (two binary64); a failed pair stores `Dx_r = NaN`. Failed chunks follow the rule above (blank rows from the failure on).

Group `pooled/<id>` (chunk level, `stretching` with a `counts` subset): attribute `orbits` (the number of orbits pooled) always; when all pooled orbits share one cell range (bitwise equal `range`), also `counts` (int64, (bins + 2) x 2: the exact integer sum of the orbits' counts at the last grid step) and the attribute `range` (two binary64). When they do not share it, the group holds `orbits` only (D88). The write compares the pooled and provenance groups on the read-back like the orbit data.

## Added by Amendment D94 (version stays 2)

No change of the file layout. In a plan with rounds (`experiments-spec.md`, "Amendment D94", E41) the cell named in `chunk_key` is `k<i>:<class>:r<round>` (K = 0 cells alike; controls stay `controls`); a plan without rounds keeps the names `k<i>:<class>`. No orbit attribute `round` is added (P70): the round of an orbit lives in the panel file only. Chunk ordinals of a plan with rounds follow E41 (each_round arms, then closure arms; round ascending).
