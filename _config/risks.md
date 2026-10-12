# Risks

Layer 3, canonical. Reviewed at every milestone exit (see `build-order.md`). Likelihood and impact: L (low), M (medium), H (high). Close a risk by marking it `closed` with the date; never delete it.

| ID | Risk | Likelihood | Impact | Mitigation | Watch at |
|----|------|------------|--------|------------|----------|
| K1 | Julia not installed on the main machine | H | H (blocks all code) | first step of M0: install via juliaup (`setup-guides/julia.md`) | closed 2026-10-12: Julia 1.13.1 installed via juliaup, channel `release` |
| K2 | celmech does not build on Windows or with REBOUND 5 | M | M (P6 models) | owner requires native Windows: if the build fails, fix it (compiler flags, patch, or local wheel); WSL2 only as last resort | closed 2026-10-12: both build and pass tests natively with `python/nbody/windows-build/` (`setup-guides/rebound.md`, section "Native Windows build"); rebuild REBOUNDx after each REBOUND upgrade |
| K3 | ASSIST pins `rebound<5`, conflicting with venv A | H | L | separate venv B or WSL2 | M8 |
| K4 | WHFast512 sources conflict on native Windows AVX512 support | M | L | verify on this CPU; else run under WSL2 or Linux | M8 |
| K5 | Geist1990 CPU-time ratios unverified (full text inaccessible) | H | L | reproduction target only after reading the paper; never a test oracle. M1: full text still not obtained (not open access) | M3 |
| K6 | Open citations: Lu1997 unmatched, CalvoSanzSerna1993 a or b, LGALI source | M | L | resolve at each module's 01-spec; marked `?` until then. M1 (2026-10-12): ReinTamayo2018 confirmed by Crossref and arXiv (doi:10.1093/mnras/stx2479), N7 resolved with the owner's approval; both Calvo and Sanz-Serna 1993 papers found, no abstract says which holds the coefficients (S9, integrator 01-spec); no source defines LGALI, three candidates listed in `chaos-indicators-candidates.md` (I4, indicator 01-spec); Lu1997 not searched, still unmatched | M2, M5 |
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
| K17 | M1 searches stopped at the extra-query cap before saturation in several sub-questions; ADS was not used (no token) | M | M (a relevant paper may be missing) | each module's 01-spec reruns a targeted search for its own rows; set `ADS_API_TOKEN` to add ADS | every module 01-spec |
| K18 | Published standard-map lambda(K) values are sparse: mostly fits without error bars, single orbits, unstated precision (see `benchmarks.md`) | H | M (weak external regression targets for P1) | P1 uses its own high-precision reference runs (BigFloat, TaylorIntegration) as ground truth; literature values serve as sanity checks with stated tolerance only | M2, M3 |
