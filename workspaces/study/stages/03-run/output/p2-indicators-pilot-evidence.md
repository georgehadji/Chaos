# p2-indicators: pilot and K19 evidence (stage 03, before any panel or grid data)

Date 2026-10-29. Not study data (P1 Amendment A4 item 2, applied by "Amendment rule" of `p2-indicators-preregistration.md`). Preregistration hash in force `33445fd7f445c18a9c0d816ead1f248b03fb73266c8a6eb344a4c71ecaa58266` (after Amendment A1), code at commit be6a425. Pilot: `bin/chaoslab.jl pilot` at 23 threads, 4 non-panel starts per K index 1 to 6 (24 starts, ids 3188 to 3211), horizon 1e5, store outside `data/`. Files: `p2-indicators-pilot-report.json`, `p2-indicators-k19-records.toml` (this folder). Nothing was written under `data/`; `data/gpu-benchmark.toml` has the same SHA-256 before and after (ec33c0b7...368a).

## 1. Cost and budget (decision input, no amendment needed)

| Item | Pilot | Preregistration |
|------|-------|-----------------|
| Projected thread time, all 17 arms | 24830 s = 6.9 thread-hours | 54,753 s = 15.2 thread-hours before safety factors, 30.8 with them |
| Projected wall time | 1618 s = 0.45 h at parallel efficiency 0.667 | 1.3 to 2.7 h, wall cap 12 h |
| Plan | 26964 items, 695 chunks, 9.04e10 steps | same |

Measured thread-ns per orbit-step (Float64) against the spec lines: tg 300 (line 1000), pr 274 (13000), st 105 (1100), pa1 300 (1000), bk1 504 (1000), fq1 12609 (13400), L1 168 and L5 182 (207.6). Double64: L1 1520 (6634.6), tg 3150 (6000), fq1 17394 (28100), pa1 2184 (4200), bk1 3304 (5000). Every measured cost is under its line, so no arm is near the first-chunk cap (2 x line). The main drivers are the L1/L5 arms (lyap_g and lyap_p: about 19.1 of the 24.8 thousand thread-seconds) and the E1 transform (12.6 microseconds per orbit-step; in the pilot chunk the transform took 25.4 s against 0.27 s for the orbit pass).

Flag shares in the pilot (24 starts, Float64, W up to 51200): fq1 below_floor 0.42, mirror 0.42, near_zero 0.13; bk1 below_floor 0.42, reldig_undefined 0.38; pa1 below_floor 0.42. These are the regular starts at the floor (K = 0 and K = 0.5: 8 of 24 starts, plus regular Kc starts).

## 2. K19: per-spec device records (23 threads, 1024 orbits per spec, best of 2)

| Spec | Kind | Device ns per orbit-step | CPU (23 threads, wall) ns per orbit-step | Device / CPU | gpu_faster |
|------|------|--------------------------|-------------------------------------------|--------------|------------|
| tg | tangent | 43.9 | 8.80 | 4.99 | false |
| pr | pair | 55.3 | 9.15 | 6.05 | false |
| pa1 | phase_advance | 72.2 | 9.48 | 7.62 | false |
| bk1 | birkhoff | 185.1 | 21.19 | 8.74 | false |
| st, fq1 | stretching, frequency | CPU only (reducer is CPU-only) | 4.77, 326.3 | | |

On this device and thread count no device-eligible spec is faster on the GPU, so `auto` resolves every spec to the CPU and the budget of the preregistration is the CPU budget. The current `data/gpu-benchmark.toml` has no entry for these spec ids either (a missing entry means false), so the backend table is the same with or without installing this file; the file is therefore kept here as the K19 record and not copied into `data/`. Takes vary (the D86 verification saw pair and phase advance faster on the device at other horizons); the record of the real run is whatever its manifest holds.

## 3. Evidence on tables R1 and R2

Pilot orbits by the growth ratios and lambda_1: regular-looking (r < 3.1 at both budgets): 10 orbits (K index 1 and 2, and two of the Kc starts); weakly chaotic (lambda_1 4.6e-3 to 5.0e-3 at 1e5): 2 orbits (3199 at Kc, 3201 at K = pi/2); strongly chaotic (lambda_1 >= 0.13): 12 orbits.

Per orbit values at n = 1e5 (E1, E2, WBA at W = 51200). Columns: lambda_1, growth ratios r (LE, SALI, FLI, MEGNO), SALI, RLI, D_E1, D_E2, dig_min, then the R1 verdicts (C chaotic, r regular) for LE SALI FLI MEGNO RLI E1 E2 WBA.

```
=== budget n = 100000  (E1/E2/WBA use the window W = 51200) ===
k   orbit         lam1    rLE  rSALI   rFLI   rMEG       SALI        RLI     DE1     DE2 digmin | verdicts LE SALI FLI MEG RLI E1 E2 WBA
1   3188      1.12e-04   1.26   1.26   1.26   1.00   2.00e-10   0.00e+00  -14.15  -14.15  14.15 | r r r r r r r r
1   3189      1.12e-04   1.26   1.26   1.26   1.00   2.00e-10   0.00e+00  -14.15  -14.15  14.15 | r r r r r r r r
1   3190      1.12e-04   1.26   1.26   1.26   1.00   2.00e-10   0.00e+00  -14.15  -14.15  14.15 | r r r r r r r r
1   3191      1.12e-04   1.26   1.26   1.26   1.00   2.00e-10   0.00e+00  -14.15  -14.15  14.15 | r r r r r r r r
2   3192      7.05e-05   1.50   1.52   1.48   1.00   1.18e-06   4.78e-15  -14.15  -14.15  14.15 | r r r r r r r r
2   3193      6.75e-05   2.56   2.82   1.44   1.14   2.99e-06   1.56e-13  -14.15  -14.15  14.15 | r r r r r r r r
2   3194      1.13e-04   1.23   1.24   1.25   1.00   2.64e-10   1.20e-13  -14.15  -14.15  13.76 | r r r r r r r r
2   3195      8.90e-05   1.32   1.29   1.28   1.00   5.91e-09   7.71e-14  -14.15  -14.15  13.17 | r r r r r r r r
3   3196      1.34e-01  13.36  13.36  13.36  13.02   0.00e+00   1.55e-02   -1.08   -1.74   1.03 | C C C C C C C C
3   3197      9.82e-05   1.30   1.31   1.29   1.00   4.36e-09   6.07e-14  -14.15  -14.15  13.11 | r r r r r r r r
3   3198      1.10e-04   1.25   1.37   1.25   1.00   8.36e-08   1.23e-13  -14.15  -14.15  13.32 | r r r r r r r r
3   3199      4.96e-03   6.41   6.56   6.35   6.21   0.00e+00   2.50e-03   -9.95   -9.86   3.60 | C C C C C r r C
4   3200      2.65e-01   9.05   9.05   9.04   9.26   0.00e+00   6.52e-02   -0.48   -0.40   0.62 | C C C C C C C C
4   3201      4.59e-03   8.35   8.34   8.29   8.79   0.00e+00   7.71e-04   -8.02   -7.75   2.42 | C C C C C r r C
4   3202      3.04e-01   9.76   9.76   9.76  10.07   0.00e+00   6.24e-03   -3.79   -0.71   1.74 | C C C C C C C C
4   3203      3.05e-01  11.35  11.35  11.34  12.64   0.00e+00   1.75e-02   -3.95   -1.45   1.79 | C C C C C C C C
5   3204      1.17e+00  10.09  10.09  10.09  10.09   0.00e+00   4.54e-03   -0.66   -1.45   1.81 | C C C C C C C C
5   3205      1.17e+00   9.88   9.88   9.88   9.82   0.00e+00   7.44e-03   -1.52   -1.05   1.88 | C C C C C C C C
5   3206      1.18e+00   9.87   9.87   9.87   9.89   0.00e+00   1.97e-02   -0.42   -0.48   1.77 | C C C C C C C C
5   3207      1.17e+00   9.97   9.97   9.97   9.94   0.00e+00   5.95e-03   -0.67   -0.31   2.04 | C C C C C C C C
6   3208      1.63e+00   9.99   9.99   9.99   9.93   0.00e+00   7.12e-03   -1.93   -0.98   2.11 | C C C C C C C C
6   3209      1.62e+00  10.00  10.01  10.00   9.95   0.00e+00   1.27e-02   -0.61   -0.49   2.21 | C C C C C C C C
6   3210      1.61e+00  10.04  10.04  10.04  10.05   0.00e+00   4.31e-03   -0.68   -0.46   2.06 | C C C C C C C C
6   3211      1.61e+00  10.00  10.00  10.00   9.97   0.00e+00   6.81e-03   -0.31   -0.63   2.11 | C C C C C C C C
```

The same table at n = 12800 is in the pilot store; summary: the same verdict pattern, regular r_SALI up to 3.10 (orbit 3193), weak chaos r_FLI as low as 5.37 (orbit 3199).

Window dependence of the frequency and Birkhoff scores (D of E1 and E2, dig_min of WBA); R regular-looking, weak and C as above:

```
D_E1 by window (fq1), then D_E2 (pa1), then WBA dig_min:
orbit    cls      100    200    400    800   1600   3200   6400  12800  25600  51200
3188     R K0   -14.1  -14.1  -14.1  -14.1  -14.1  -14.1  -14.1  -14.1  -14.1  -14.1   E1
                -14.1  -14.1  -14.1  -14.1  -14.1  -14.1  -14.1  -14.1  -14.1  -14.1   E2
                  8.9   11.7   14.1   14.1   14.1   14.1   14.1   14.1   14.1   14.1   WBA dig_min
3192     R       -6.3   -8.2   -9.0  -11.3  -11.4  -12.6  -13.9  -14.1  -14.1  -14.1   E1
                 -5.7   -8.4  -11.9  -14.1  -14.1  -14.1  -14.1  -14.1  -14.1  -14.1   E2
                  4.9    7.2   10.6   14.1   14.1   14.1   14.1   14.1   14.1   14.1   WBA dig_min
3194     R       -8.3   -9.7  -10.2  -11.0  -12.2  -13.6  -14.1  -14.1  -14.1  -14.1   E1
                 -7.8   -9.1  -11.6  -14.1  -14.1  -14.1  -14.1  -14.1  -14.1  -14.1   E2
                  6.7    8.3   11.0   13.8   14.1   14.1   14.1   14.1   13.8   13.8   WBA dig_min
3195     R       -5.8   -6.5   -7.3   -7.7   -8.8  -10.1  -11.7  -12.7  -14.0  -14.1   E1
                 -5.6   -5.9   -7.2   -7.5   -8.7  -10.1  -11.5  -14.1  -14.1  -14.1   E2
                  4.2    4.2    5.4    5.7    6.9    8.3    9.7   13.6   13.9   13.2   WBA dig_min
3197     R       -5.1   -5.9   -7.6  -10.9   -9.9  -12.3  -12.4  -13.9  -14.1  -14.1   E1
                 -4.6   -6.2   -7.0  -10.0  -11.8  -14.1  -14.1  -14.1  -14.1  -14.1   E2
                  2.9    4.5    5.4    9.2   10.2   12.7   14.1   13.4   14.0   13.1   WBA dig_min
3198     R       -5.2   -6.7   -6.7   -7.9   -8.4   -9.5  -11.4  -12.5  -13.3  -14.1   E1
                 -5.1   -6.5   -7.2  -11.2  -10.6  -14.1  -14.1  -14.1  -14.1  -14.1   E2
                  4.3    5.5    5.5    6.6    6.9    8.0    9.9   11.1   13.2   13.3   WBA dig_min
3199     weak    -4.5   -5.3   -6.5   -7.0   -7.7   -7.5   -8.7   -8.9   -9.0   -9.9   E1
                 -4.1   -4.9   -5.9   -6.6   -7.1   -7.3   -8.3   -8.8   -8.9   -9.9   E2
                  2.9    3.4    3.2    2.5    4.2    3.0    3.1    3.3    3.0    3.6   WBA dig_min
3201     weak    -5.8   -5.2   -5.7   -7.8   -7.2   -9.3   -6.8   -7.0   -7.5   -8.0   E1
                 -5.2   -4.9   -5.6   -7.3   -6.9   -8.1   -6.8   -7.0   -6.8   -7.8   E2
                  2.2    2.4    2.8    4.2    4.5    4.1    2.8    3.6    3.1    2.4   WBA dig_min
3196     C       -2.3   -2.0   -2.3   -2.2   -1.7   -1.1   -0.8   -0.8   -0.3   -1.1   E1
                 -2.2   -2.0   -2.2   -2.1   -1.8   -0.7   -0.7   -1.0   -0.7   -1.7   E2
                  1.5    1.2    1.5    1.3    1.1    0.2    0.2    0.5    0.2    1.0   WBA dig_min
3200     C       -2.7   -0.7   -1.0   -0.4   -0.7   -0.8   -0.8   -0.8   -0.8   -0.5   E1
                 -1.6   -0.4   -0.7   -1.6   -0.3   -1.4   -0.4   -0.4   -0.4   -0.4   E2
                  1.0    0.1    0.4    1.1    0.8    1.0    0.7    0.7    0.7    0.6   WBA dig_min
3204     C       -0.6   -0.8   -0.6   -0.3   -2.1   -1.3   -0.6   -1.7   -1.7   -0.7   E1
                 -0.8   -0.8   -0.3   -0.4   -0.6   -0.7   -1.0   -1.1   -1.2   -1.4   E2
                  0.5    0.9    0.8    0.9    1.1    1.9    1.6    1.5    2.3    1.8   WBA dig_min
3208     C       -0.5   -0.9   -0.5   -0.6   -1.6   -0.8   -0.3   -0.6   -0.6   -1.9   E1
                 -0.8   -1.4   -0.5   -1.0   -0.9   -0.5   -0.4   -0.3   -0.4   -1.0   E2
                  1.5    1.1    1.1    1.3    1.5    1.4    2.3    2.2    2.3    2.1   WBA dig_min
```

Findings:
- F-R1-a. Ratio rule (r >= 5 for LE, SALI, FLI, MEGNO): at n = 1e5 and 12800 the pilot has regular r <= 3.10 and chaotic r >= 5.37 (weak chaos) or about 10 (strong). 5 is inside the gap at both budgets; the margin on the chaotic side is 0.37 at n = 12800. No change proposed.
- F-R1-b. LE rule: the regular K = 0 starts have lambda_1 = 1.12e-4 at 1e5, above the floor 1e-4, so the second condition (r >= 5) alone prevents a false "chaotic" (as in P1). No change.
- F-R1-c. RLI >= 1e-7: regular 0 to 1.2e-13, chaotic 7.7e-4 to 6.5e-2. Inside the gap, no change.
- F-R1-d. WBA, dig_min < 6: at W = 51200 regular 13.1 to 14.1, weak chaos 2.4 to 3.6, strong 0.6 to 2.1. Inside the gap at large W. At small windows the regular starts also fall below 6 (W = 100: 2.9 to 8.9; W = 800: 5.7 to 14.1), so WBA calls regular orbits chaotic at the small budgets. This is the finite-time effect that H5 measures, not a number error.
- F-R1-e. E1 and E2, D >= -4: the two weakly chaotic starts sit at D = -4.1 to -10 at every window and never reach -4, so under the rule E1 and E2 call them regular at every budget. Regular starts drop from about -5 to -8 at W = 100 to the floor -14.1 at large W; no fixed D separates weak chaos from regular at all windows (-11 separates them at W >= 3200 but marks regular starts chaotic below W = 800). Strong chaos is above -4, but two K = pi/2 starts (3202 and 3203, lambda_1 0.30) have D_E1 = -3.79 and -3.95 at W = 51200, margins of 0.21 and 0.05 decades; the other strong-chaos starts are at -0.3 to -1.9. The pilot therefore shows a real indicator disagreement (E1 and E2 against LE, SALI, FLI, MEGNO, RLI, WBA on weak chaos) and no better single threshold. Changing the number to make them agree would tune the rule to the outcome. No number change proposed.
- F-R1-f. The "where the number comes from" cell of E1 says -4 is placed "5 decades above the mean-removal bias (about -8.5)". The pilot shows the regular starts at the floor (-14.1) with `mean_removal = true`: the bias is the same in both windows and cancels in D. The statement is unsupported; the number stays.
- F-R2-a. R2 SALI(n) <= 1e-8: regular SALI is 2.0e-10 (K = 0) and up to 5.9e-9 at n = 1e5, below 1e-8, so R2-SALI calls regular orbits chaotic from n of about 1e4 on, as the text of R2 warns (power-law decay). R2 FLI and LE: n lambda_1 is about 11 for regular and about 500 for chaotic at 1e5 against 3 ln n = 34.5. No change (R2 is a sensitivity).

## 4. Proposed amendments (none applied; the owner decides)

1. None to the numbers of table R1 or R2 (recommended: the evidence is inside the gaps, or is a result the study exists to measure).
2. Optional text-only correction (Amendment A2): the E1 and E2 rows of table R1 name the mean-removal bias as the reason for -4 (F-R1-f). The corrected reason would read: -4 lies about 10 decades above the regular floor and below every strong-chaos start of the pilot (the nearest at -3.95 and -3.79, margins of 0.05 and 0.21 decades); the number is unchanged.
3. Budget lines: no change (every measured cost is under its line; the projection is 6.9 thread-hours against 30.8 with safety).
4. K19: no file in `data/` changes; the record above is kept in this folder.

Limits: 24 starts (2 weakly chaotic, 12 strongly chaotic, 10 regular-looking), horizon 1e5 instead of 102400, no reference labels, so "weak chaos" is by lambda_1 and the growth ratios, not by the P1 label rules. The pilot cannot validate the sticky definition T_s or the grid.
