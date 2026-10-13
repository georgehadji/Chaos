# Environment

Layer 3, canonical. Filled once by `setup`.

## Compute

- Target: {{COMPUTE_TARGET}}
- Max cores: {{MAX_CORES}}

## Storage

- Data root: {{DATA_ROOT}} (run stores and literature full texts; never committed)

## High precision

- Backend for extended precision: Double64 (DoubleFloats.jl, decision D20; effective eps measured 2^-104.4 to 2^-104.9, precision stage 04). BigFloat is always available for oracle runs.
