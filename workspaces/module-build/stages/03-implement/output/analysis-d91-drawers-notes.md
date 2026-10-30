# Analysis Amendment D91: drawers F1 to F6 of study p2-indicators (stage 03 notes)

Written by the main session from the report of the drawers subagent (the subagent could not write the file). Spec: `analysis-spec.md`, Amendment D91. Module side: `analysis-d91-impl-notes.md`.

## Files

- `ext/p2_figures.jl` (243 lines: helpers, F1 to F3) and `ext/p2_figures_b.jl` (174 lines: F4 to F6, driver `draw_p2_figures`), both included from `ext/ChaoslabCairoMakieExt.jl`.
- `ext/ChaoslabCairoMakieExt.jl`: two `include` lines, one dispatch line at the top of `draw_figures` (study p2-indicators goes to `draw_p2_figures`), one docstring sentence. The P1 drawers, `HALVES`, the constants and `savefig` are untouched (code review of the diff).
- No function over 34 lines, no file over 400 lines, no `@fastmath`, no em dash; the extension imports only Chaoslab and CairoMakie (I29). No change of `Project.toml`.

## Measured (real tables from `analyze_p2`, store read only)

One page each; smallest text 16.00 pt on the page; scaled to 504 pt:

| Fig | Page (pt) | Bytes | Smallest text scaled |
|---|---|---|---|
| F1 | 996 x 625 | 15,695 | 8.10 |
| F2 | 1000 x 900 | 22,382 | 8.06 |
| F3 | 1000 x 800 | 23,807 | 8.06 |
| F4 | 1000 x 830 | 44,116 | 8.06 |
| F5 | 994 x 1022 | 162,835 | 8.11 |
| F6 | 1000 x 560 | 42,109 | 8.06 |

The six drawers take 28 s on the first call (compile) and 1.1 s warm. The real command `bin/chaoslab.jl analyze` printed "figures written", exit 0, six PDFs, and a tables file equal to the accepted one except the two header lines. F5 read back from the PDF at 8 px per pt: all 41,472 cell centres equal the block state (tolerance 2 of 255), exactly 3 colours, i along x, j upward. Cairo writes each cell as a vector rectangle, so a count of the 144 maps must count paths, not images. PDF colours read back one level off the hex values, so the verify script needs a tolerance of at least one level. Text content was checked from the PDFs for all six (names, sets, budgets, K values, counts, the E42 text twice in F4).

## Design

- Every figure is drawn with one unit equal to one point on a page of at most 1000 pt, all text 16 units. At 1000 pt the smallest text scaled to 504 pt is 8.06, so any smaller font fails D91f (margin 0.06 pt).
- Three states: C orange #E69F00, R dark blue #1B4F9C, U light grey #E2E2E2. Checked by grey-scale and deuteranopia simulation on F3 and F5. Eight indicator colours plus eight marker shapes, each at least 25 CIELAB from each class colour; in grey-scale the indicator colours merge and the marker carries the identity. STR is black dashed.
- F1: transposed 8 x 8 heat maps, pair drawn once above the diagonal; A_XY viridis on [0,1], kappa_XY vik on [-1,1], one colour bar per row; titles in the top row only; the RLI note gives N per set.
- F2: log-x step lines of the running share of t_det with a tick at every budget, bottom row three-grey stacked outcome bars.
- F3: 2 x 2 panels, b_P dashed and N_c dotted, named in the legend with their values; RLI is the only legend entry with its own decided count.
- F4: pooled rows, 2 x 5 panels, "nothing" is a gap, an all-"nothing" panel shows only the E42 text.
- F5: one page-sized canvas axis with 144 images of one square cell each (3.84 pt), no interpolation.
- F6: log-log; p5-p95 and p25-p75 bands, median line with markers; T_s dark grey dashed; the legend flags a count that changes with the step and notes that T_s is a subset of C.
- A value that cannot be drawn raises an `ErrorException` naming the figure and row or cause (for example "F6 row 12: step value 0 cannot be drawn on a logarithmic axis"); saving goes through a temp file removed on error, so no half-written file stays.

## Fixture edge cases drawn (hand-edited copies of the real tables)

No RLI in S2; a K with every AUC "nothing" (mix of `nothing` and the string); one K index in F5; F6 count varying with step; one budget in F1; no `res.params`: all drawn without error, all six files. A zero p50 or step in F6, a NaN kappa in F1, a bad F5 state and a missing F3 table raise with the figure and row named and leave no file for that figure.

## Deviations and ambiguities

1. The F3 fourth set is named "regular" in the block; the drawer keeps the block name for ordering and writes "regular set" in the title.
2. F5: the lattice extent text is a line under each group title, not on the first map (the text does not fit a 92 pt map). Row titles "K index k" are drawn in both groups.
3. F1 note: N given per set (E39 left "N = ..." open).
4. F6: T_s is dark grey with its own marker, to keep one colour one meaning (I28).
5. Y limits of -0.03 to 1.03 in F2, F3, F4 so that markers at 0 and 1 are not clipped.
6. F3 legend annotates only RLI; LE, MEGNO, E1 and WBA also have fewer decided orbits in the regular set (1145 to 1154 of 1159) and are not annotated.
7. K text uses 12 significant digits (the tables-file format).
8. D91g "exactly once where a title is concerned": F1 titles repeat the set per budget column by design (W2 appears in two titles plus the note).

## Not checked

Viewers other than PyMuPDF (possible hairline seams between F5 rectangles); D91k determinism; the O12 plotted-number re-read; memory of the figure step; the P1 pictures (code review only); colour-vision deficiency beyond the simulation on F3 and F5. The owner visual check decides the look.
