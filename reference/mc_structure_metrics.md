# Recovery of an undirected support

Compares an estimated set of edges with a known one over the
\\R(R-1)/2\\ pairs of responses and returns the usual counts and rates.
With true positives \\TP\\, false positives \\FP\\, false negatives
\\FN\\ and true negatives \\TN\\,
\$\$\mathrm{TPR}=\frac{TP}{TP+FN},\quad
\mathrm{FPR}=\frac{FP}{FP+TN},\quad
\mathrm{precision}=\frac{TP}{TP+FP},\quad F_1=\frac{2TP}{2TP+FP+FN},\$\$
and the structural Hamming distance is \\\mathrm{SHD}=FP+FN\\. The
comparison ignores direction: an adjacency is reduced to its skeleton
first (a pair is an edge when either entry is nonzero).

## Usage

``` r
mc_structure_metrics(estimated, truth, R = NULL)
```

## Arguments

- estimated, truth:

  an adjacency matrix, a logical (or 0/1) indicator over the pairs in
  upper-triangular order, or a two-column matrix of edges.

- R:

  number of responses, needed only when an argument is an edge list;
  when omitted it is taken from `truth` if that is a square matrix.

## Value

A named numeric vector with `tpr`, `fpr`, `precision`, `f1`, the
structural Hamming distance `shd`, and the counts `tp`, `fp`, `fn` and
`tn`. A rate whose denominator is empty is `NA`: precision is undefined
when nothing was selected, the true positive rate when there was nothing
to find, and the false positive rate when the true graph is complete.
\\F_1\\ is `NA` only when both edge sets are empty; it is zero when
nothing was selected and the truth has edges.

## See also

[`mc_cpdag_metrics`](https://evandeilton.github.io/mcglmgraph/reference/mc_cpdag_metrics.md)
for partially directed structures.

Other structure recovery:
[`mc_cpdag_metrics()`](https://evandeilton.github.io/mcglmgraph/reference/mc_cpdag_metrics.md),
[`mc_true_cpdag()`](https://evandeilton.github.io/mcglmgraph/reference/mc_true_cpdag.md)

## Examples

``` r
truth <- rbind(c(1, 2), c(2, 3))
est <- rbind(c(1, 2), c(1, 3))
mc_structure_metrics(est, truth, R = 3)
#>       tpr       fpr precision        f1       shd        tp        fp        fn 
#>       0.5       1.0       0.5       0.5       2.0       1.0       1.0       1.0 
#>        tn 
#>       0.0 
```
