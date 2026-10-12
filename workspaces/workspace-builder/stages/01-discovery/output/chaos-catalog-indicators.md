# Catalog: Chaos Indicators

> Promoted to `/_config/catalogs/indicators.md` on 2026-10-12; that copy is canonical. This file is a Stage 01 record.

Layer 3 draft. Indicators classify finite-time dynamics; they are not LE estimators (claim C3). Keys with `?` are pending verification. Decision rules are typical behaviour from the literature; the thresholds actually used are pre-registered per study, never taken from this table.

## A. Variational, deviation-vector based

| ID | Indicator | Vectors | Regular orbit | Chaotic orbit | Main issue | Source | Phase |
|----|-----------|---------|---------------|---------------|-----------|--------|-------|
| I1 | Finite-time mLCE lambda_1(n) | 1 | -> 0 like ln n / n | -> lambda_1 > 0 | slow; threshold choice | Benettin1980a, Skokos2010 | M2 |
| I2 | SALI (min of norm of sum / difference of two unit deviation vectors) | 2 | fluctuates around positive value | -> 0 like exp(-(lambda_1 - lambda_2) n) | in 2D maps: exp(-2 lambda_1 n), cross-check vs I1 | Skokos2001, Skokos2004 | M4 |
| I3 | GALI_k (volume of k unit deviation vectors, product of singular values) | k = 2..2N | constant or power-law decay (torus dimension) | exponential decay with sum of exponent gaps | k > 2 needs dim >= 4 | Skokos2007, Skokos2008 | M4 (k=2), M5 (k>2) |
| I4 | LGALI / GALI for local dynamics | k | see I3 | see I3 | | Skokos2008? (LGALI source not confirmed, resolve at 01-spec) | M5 |
| I5 | FLI (log of largest deviation norm up to time t) | 1 or full basis | linear growth in ln t | exponential growth | initial vector choice, saturation | Froeschle1997 | M4 |
| I6 | OFLI (orthogonal component to the flow) | 1 | distinguishes periodic vs quasi-periodic | as FLI | flows only (tangent to flow removed) | Fouchard2002 | M4 |
| I7 | OFLI2 (second-order variational) | 1 + 2nd order | sharper periodic detection | | needs 2nd variational equations | Barrio2005 | M5 (catalog only) |
| I8 | MEGNO Y(t) and mean <Y>(t) | 1 | <Y> -> 2 (quasi-periodic), 0 (stable periodic) | <Y> ~ lambda_1 t / 2 | transients; analytic link to FLI | CincottaSimo2000, Cincotta2003, Mestre2011 | M4 |
| I9 | RLI (difference of finite-time LEs of two very close orbits) | 1 per orbit | ~ 0 | grows | initial separation dependence | Sandor2004 | M4 |
| I10 | Stretching numbers, spectra, spectral distance | 1 | narrow invariant spectrum | broad spectrum | window parameters | Voglis1994, Contopoulos2002 | M4 |
| I11 | APLE / other variants | | | | | catalog only | none |

## B. Non-variational comparators (different diagnostic constructs)

| ID | Indicator | Measures | Main issue | Source | Phase |
|----|-----------|----------|-----------|--------|-------|
| I12 | Frequency map analysis (NAFF): fundamental frequency vector nu on windows [0,T] and [T,2T] of z_j = q_j + i p_j; index D_T = log10 ||nu_[0,T] - nu_[T,2T]|| | resonance web, frequency diffusion, transport-relevant weak chaos | window length, frequency resolution, sticky orbits look quasiperiodic for long times | Laskar1990, LaskarFroeschleCelletti1992, Laskar1993 | M4 |
| I12b | Refined Fourier / NAFF variants (windowed, interpolated, higher-precision frequency estimators) | frequency precision per window | leakage, close resonances | Gomez2010a, Gomez2010b, Bartolini1996 | M4 |
| I12c | Power-spectrum test (PSD of q_i, p_i or q_i + i p_i) | discrete lines (regular) vs broadband (chaotic) | leakage without windowing; qualitative unless scored | Simo2014 | M4 |
| I12d | Time-frequency analysis (short-time Fourier, wavelet ridges, instantaneous frequencies) | resonance capture, intermittency, chaotic transitions | time-frequency resolution trade-off | Chandre2003 | M4 |
| I13 | Weighted Birkhoff averages (convergence rate) | super-convergence on tori vs slow in chaos | weight function, digits threshold | Das2017 | M4 |
| I14 | Reversibility error method (REM) | sensitivity via forward-backward error | precision-dependent by design | Panichi2017 | M4 |
| I15 | 0-1 test for chaos | K statistic -> 0 or 1 | sampling, regular orbits in conservative systems | GottwaldMelbourne2004 | M6 |
| I16 | Recurrence quantification (RQA) | recurrence structure | threshold epsilon | Marwan2007 | M6 (catalog only) |
| I17 | Lagrangian descriptors | phase-space structure | integration horizon | Mancho2013 | catalog only |
| I18 | ML classifiers on short orbits | learned classification | data leakage, generalization | Celletti2022 | catalog only |

Implementation: NAFF, FMA (windows, D_T, grids) and time-frequency analysis are written in Julia (FFTW.jl) because the core is Julia. Oracles via HDF5 on the Python side: `nafflib` (xsuite / CERN, 2.1.1, Apache-2.0), `PyNAFF` (1.2.0) and REBOUND's `frequency_analysis` (MFT); celmech also ships an FMFT routine.

Fourier rules: never a raw FFT alone; use windowing (e.g. Hanning) and NAFF-type refinement; D_T is a finite-time diagnostic, cross-checked against a tangent-space indicator (FLI, MEGNO or lambda_1(n)) on the same orbit. For the standard map, FMA is applied to the iterates directly (LaskarFroeschleCelletti1992 applies it to the standard map).

## C. What every indicator run records

Full time series (not final value); detection time to stable pre-registered label; cost (steps, vectors, wall time); sensitivity to initial deviation vector (random orientation ensemble); precision.

## D. Literature recommendations to test, not assume

Maffione2011a (arXiv:1108.2196) recommends FLI + RLI for wide phase-space scans and MEGNO + SALI for single-orbit characterization (abstract, tested on 4D maps); companion studies Maffione2011b and Darriba2012 (indicators vs ODE integrators). P2 tests this on the standard-map panel.
