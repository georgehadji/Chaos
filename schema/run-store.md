# Run Store Schema (HDF5), version 1

Layer 3, canonical. Shared contract between Julia (`Experiments`) and Python (`python/nbody/`). The first `Experiments` module spec may refine it; any change bumps `schema_version`.

## Files

- One HDF5 file per chunk: `<data root>/[study]/chunk-<NNNNN>.h5` (zero-padded index).
- Chunks are independent and resumable; a chunk is complete when its root attribute `complete = true`.

## Root attributes (provenance)

| Attribute | Type | Example |
|-----------|------|---------|
| `schema_version` | int | `1` |
| `study` | string | `p1-le-benchmark` |
| `chunk` | int | `17` |
| `git_sha` | string | 40-hex commit of a clean tree |
| `prereg_sha256` | string | hash of the frozen preregistration |
| `config_sha256` | string | hash of the run config |
| `julia_version` / `python_version` | string | writer runtime |
| `packages` | string | manifest hash or `pip freeze` hash |
| `cpu`, `threads` | string, int | host info |
| `precision` | string | `Float64`, `Double64`, `BigFloat(256)` |
| `created_utc` | string | ISO 8601, e.g. `2026-10-12T14:03:00Z` |
| `complete` | bool | `true` |

## Groups

```
/orbits/<orbit_id>/            one group per initial condition
    attrs: system, params, ic, class_label, seed
    t                          sample times or step counts, length n
    lyap/<method>              n x k finite-time exponents per method
    indicators/<name>          n samples per indicator
    diagnostics/<name>         n samples (pairing residual, orthogonality loss, energy error)
    cost/<method>              steps, vectors, wall seconds
/grid/                         grid scans only
    attrs: system, params, axes
    <name>                     final values, shape (nx, ny)
```

N-body files use `system = "nbody:<name>"` and the same groups.
