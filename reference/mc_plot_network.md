# Plot the estimated dependence network

Draws the response variables as nodes on a circular layout. GGM edges
are weighted by partial correlations; DAG arrows are weighted by the
structural coefficients \\\gamma\\ of the ordered modified Cholesky
factorization, and each arrow stops at the border of the node it points
to. Edge width encodes magnitude, saturating at an absolute weight of
one, and color encodes sign (blue positive, red negative). The graphical
parameters are restored on exit.

## Usage

``` r
mc_plot_network(fit, threshold = 0.05, edge_labels = FALSE)
```

## Arguments

- fit:

  an object of class `"mcglm_graph"`.

- threshold:

  an edge is drawn when the absolute value of its weight exceeds
  `threshold`: the partial correlation for a GGM fit, the coefficient
  \\\gamma\\ for a DAG fit. The coefficient is on the scale of the
  Pearson residuals and is not bounded by one, so the same value of
  `threshold` does not select comparably in the two layers.

- edge_labels:

  logical; if `TRUE`, annotate edges with the rounded weight.

## Value

Invisibly, the \\R \times R\\ matrix of weights (partial correlations
for the GGM layer; for the DAG layer, \\\gamma\\ in the position \\(i,
j)\\ of each arrow from \\i\\ to \\j\\, zero elsewhere), with the
response names as dimnames. The plot is drawn on the current graphics
device.

## See also

[`mc_get_graph`](https://evandeilton.github.io/mcglmgraph/reference/mc_get_graph.md),
[`mc_partial_cor`](https://evandeilton.github.io/mcglmgraph/reference/mc_partial_cor.md)

Other graph interpretation:
[`mc_get_graph()`](https://evandeilton.github.io/mcglmgraph/reference/mc_get_graph.md),
[`mc_partial_cor()`](https://evandeilton.github.io/mcglmgraph/reference/mc_partial_cor.md),
[`mc_partial_cor_delta()`](https://evandeilton.github.io/mcglmgraph/reference/mc_partial_cor_delta.md)

## Examples

``` r
set.seed(1); n <- 200; x <- rnorm(n)
d <- data.frame(y1 = rnorm(n), y2 = rnorm(n), y3 = rnorm(n), x = x)
f <- mcglm_graph(list(y1~x, y2~x, y3~x), data = d)
#> Automatic initial values selected.
mc_plot_network(f)
```
