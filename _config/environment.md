# Environment

Layer 3, canonical. Filled once by `setup`.

## Compute

- Target: {{COMPUTE_TARGET}}
- Max cores: {{MAX_CORES}}

## Local GPU

- Device of every [local] GPU criterion (decision D46): NVIDIA GeForce RTX 5060 Laptop GPU, compute capability 12.0, 8 GB, driver 573.24; CUDA.jl 6.4.2 in the verification environment (`verify/`); Windows TdrDelay 8 s.
- Host: Intel Core Ultra 9 275HX, 24 logical cores, Windows 11.
- Throughput record: `data/gpu-benchmark.toml` (per machine, git-ignored), written by `verify/gpu.jl`; Float64 GPU 2.47 ns per orbit-step, one wave 9.14 ns (2026-10-17).
- Before D46 the local GPU was a GeForce GTX 1050 (3 GB, compute capability 6.1); its runs stay in the stage 04 reports.

## Storage

- Data root: {{DATA_ROOT}} (run stores and literature full texts; never committed)

## High precision

- Backend for extended precision: Double64 (DoubleFloats.jl, decision D20; effective eps measured 2^-104.4 to 2^-104.9, precision stage 04). BigFloat is always available for oracle runs.
