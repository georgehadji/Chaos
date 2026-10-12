# Risks

Layer 3, canonical. Reviewed at every milestone exit (see `build-order.md`). Likelihood and impact: L (low), M (medium), H (high). Close a risk by marking it `closed` with the date; never delete it.

| ID | Risk | Likelihood | Impact | Mitigation | Watch at |
|----|------|------------|--------|------------|----------|
| K1 | Julia not installed on the main machine | H | H (blocks all code) | first step of M0: install via juliaup (`setup-guides/julia.md`) | closed 2026-10-12: Julia 1.13.1 installed via juliaup, channel `release` |
| K2 | celmech does not build on Windows or with REBOUND 5 | M | M (P6 models) | owner requires native Windows: if the build fails, fix it (compiler flags, patch, or local wheel); WSL2 only as last resort | M8 |
| K3 | ASSIST pins `rebound<5`, conflicting with venv A | H | L | separate venv B or WSL2 | M8 |
| K4 | WHFast512 sources conflict on native Windows AVX512 support | M | L | verify on this CPU; else run under WSL2 or Linux | M8 |
| K5 | Geist1990 CPU-time ratios unverified (full text inaccessible) | H | L | reproduction target only after reading the paper; never a test oracle | M1 |
| K6 | Open citations: Lu1997 unmatched, ReinTamayo2018 assumed, CalvoSanzSerna1993 a or b, LGALI source | M | L | resolve at each module's 01-spec; marked `?` until then | M2, M5 |
| K7 | P5 fails its falsification criteria | M | M (P6 dropped) | P6 is conditional (D11); a failed P5 is still reported (value type NEGATIVE) | M7 |
| K8 | Compute budget: 10^6 to 10^7 orbits over long runs exceeds one machine | M | M | preregistered budget checked against `environment.md`; chunked resumable runs; cluster job arrays | M3 |
| K9 | Prior art undercuts a novelty claim (e.g. existing indicator benchmarks) | M | H | literature-review prior-art check before any NOVEL value type | every study/01 |
| K10 | Scope creep: catalogs list far more methods than the papers need | H | M | build only what `build-order.md` names; other rows stay `catalog only` | every milestone |
| K11 | AI-written code with subtle numerical errors | M | H | tests first (RED), oracles, independent reviewer, human math check at 01-spec | every module |
| K12 | Dependency drift changes numbers between runs | M | M | committed `Manifest.toml`, versions in every run file, CI on Windows and Ubuntu | M0 onward |
| K13 | Run data lost (`data/` is not in git) | L | M | manifests tracked; every result regenerable from manifest and commit; back up the data root | M3 onward |
| K14 | Single researcher: decisions and context live in one head | M | M | decision log, ICM stage outputs, per-run notes | every milestone |
| K15 | The main machine's GPU (GeForce GTX 1050, 3 GB, compute capability 6.1, from `nvidia-smi`) may be unsupported by current CUDA.jl or slower than the CPU in Float64 | M | M (GPU-first gives no speedup here) | `tangent-core` 01-spec verifies CUDA.jl support for this device and its Julia 1.10 compat; benchmark CPU vs GPU and record it; `auto` falls back to CPU | M2 |
| K16 | GPU and CPU runs of a chaotic orbit diverge (FMA contraction, reduction order) | H | M (results differ by backend) | compare estimators at the preregistered tolerance, never trajectories; backend recorded per run; CPU is the reference | M2 onward |
