# Alt text for the figures

Alt text is plain description of what each figure draws. It does not interpret the results.

## Fig. 1 (p1-le-benchmark-fig-F1.pdf)

A grid of small line plots, four panels per row, one panel for each populated cell of the panel (design class and K index, for example "class C, K index 2"). In each panel the horizontal axis is the number of map steps n on a logarithmic scale from 1 to about ten million, and the vertical axis is the finite-time exponent lambda_1(n). Nine colored lines show the nine estimator methods (L1, L3, L4, L4b, L5, L5b, L6, L9 and L2) for the first drawn orbit of that cell. A gray shaded band marks the range between the smallest and the largest value of the eight methods other than L2. A legend at the right lists the methods and the band.

## Fig. 2 (p1-le-benchmark-fig-F2.pdf)

A grid of small plots, four panels per row, one panel for each populated cell and a last panel for the pooled weak-chaos set W. In each panel the horizontal axis is n on a logarithmic scale and the vertical axis is the relative spread R(n). A blue line shows the median of R(n) over the orbits of the cell, and a light blue band shows the 95 percent bootstrap interval. The vertical axis ranges differ from panel to panel. A legend at the right names the line and the band.

## Fig. 3 (p1-le-benchmark-fig-F3.pdf)

One plot of empirical distribution functions. The horizontal axis is the convergence time n_conv on a logarithmic scale, and the vertical axis is the fraction of orbits, from 0 to 1. Each curve is a step function for one cell. Color gives the design class (C, R, S, T, X) and line style gives the K index (1 to 6). A legend at the right lists the classes and the line styles.

## Fig. 4 (p1-le-benchmark-fig-F4.pdf)

A bar chart titled "label stability". The vertical axis is the share of orbits, from 0 to about 0.4. There are seven bars. Five blue bars are labeled ladder, direction, estimator, diagnostics and excluded, and show the share of orbits forced to unresolved by each test among the orbits the test applies to. Two gray bars are labeled "unresolved before the tests" and "unresolved after the tests", and show the share of all orbits that are unresolved. The ladder bar is the tallest, the direction bar is shorter, the estimator, diagnostics and excluded bars have no visible height, and the gray "after" bar is taller than the gray "before" bar.

## Fig. 5 (p1-le-benchmark-fig-F5.pdf)

Five line plots in two rows, titled det, symp, orth, pair and sum. In each plot the horizontal axis is n on a logarithmic scale and the vertical axis is the largest running maximum over the orbits of that diagnostic, on a logarithmic scale. Colored lines show the methods (L1, L3, L4, L4b, L5, L5b, L6, L9, L2). A dashed black horizontal line in each plot marks the median threshold over the orbits. The lines are flat or rise in small steps. A legend at the right lists the methods and the threshold line.

## Fig. 6 (p1-le-benchmark-fig-F6.pdf)

A grouped bar chart titled "cost". The horizontal axis lists the methods L1, L2, L3, L4, L4b, L5, L5b, L6 and L9. The vertical axis is thread nanoseconds per orbit-step on a logarithmic scale from about 100 to 100000. Within each method, bars are colored by precision: green for Float64, orange for Double64 and blue for BigFloat(256). Not every method has a bar for every precision. A legend at the right names the three precisions.

## Fig. 7 (p1-le-benchmark-fig-F7.pdf)

Four dot plots in a two-by-two grid, titled direction, eps, tau and wrapped. In each plot every blue dot is one item, labeled along the horizontal axis (for example "0.3 L1", "1.0e-12 L2", "8 L5", "wrapped L9"). The vertical axis is the median paired difference of lambda at n_H2 from the primary item. The direction plot has twelve dots (four direction offsets times L1, L5 and L9), the eps plot has three dots (all for L2), the tau plot has eight dots (tau of 2, 4, 8 and 16 for L4 and L5) and the wrapped plot has eight dots (one per method), all at zero. Vertical axis scales differ between plots.
