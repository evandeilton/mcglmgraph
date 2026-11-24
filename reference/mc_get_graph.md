# Extract the estimated dependence graph

Thresholds the fitted graph layer into an adjacency matrix and an edge
list.

## Usage

``` r
mc_get_graph(fit, threshold = 1e-04)
```

## Arguments

- fit:

  an object of class `"mcglm_graph"`.

- threshold:

  an edge is kept when the absolute value of its statistic exceeds
  `threshold`. The statistic depends on the layer: for a GGM fit it is
  the partial correlation of
  [`mc_partial_cor`](https://evandeilton.github.io/mcglmgraph/reference/mc_partial_cor.md);
  for a DAG fit it is the raw structural coefficient \\\gamma\\, which
  is on the scale of the Pearson residuals and is not a correlation.
  Only edges of the fitted support can be kept.

## Value

A list with the `adjacency` matrix (\\R \times R\\, 0/1 integer, with
the response names as dimnames), the `edges` matrix and `directed`. For
GGM fits the adjacency is symmetric and `edges` has one row `(i, j)`,
\\i \< j\\, per retained edge. For DAG fits the adjacency has a one in
position \\(i, j)\\ for an arrow from \\i\\ to \\j\\, directed from each
earlier response to its later response in the specified ordering, and
`edges` has the columns `from` and `to`. An empty selection gives an
edge matrix with no rows.

## See also

[`mc_partial_cor`](https://evandeilton.github.io/mcglmgraph/reference/mc_partial_cor.md),
[`mc_plot_network`](https://evandeilton.github.io/mcglmgraph/reference/mc_plot_network.md)

Other graph interpretation:
[`mc_partial_cor()`](https://evandeilton.github.io/mcglmgraph/reference/mc_partial_cor.md),
[`mc_partial_cor_delta()`](https://evandeilton.github.io/mcglmgraph/reference/mc_partial_cor_delta.md),
[`mc_plot_network()`](https://evandeilton.github.io/mcglmgraph/reference/mc_plot_network.md)

## Examples

``` r
set.seed(1); n <- 200; x <- rnorm(n)
d <- data.frame(y1 = rnorm(n), y2 = rnorm(n), y3 = rnorm(n), x = x)
f <- mcglm_graph(list(y1~x, y2~x, y3~x), data = d)
#> Automatic initial values selected.
mc_get_graph(f)$adjacency
#>    y1 y2 y3
#> y1  0  1  1
#> y2  1  0  1
#> y3  1  1  0
```
