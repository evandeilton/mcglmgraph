# Build the between-response covariance and its derivatives under a graph

Maps a vector of graph parameters to the between-response covariance
\\\Sigma_b\\ and to the list of its derivatives with respect to each
parameter, in the list form consumed by the estimating-function engine
of
[`mcglm_graph`](https://evandeilton.github.io/mcglmgraph/reference/mcglm_graph.md).
It takes the place that the correlation-based between-response builder
has in a multivariate covariance generalized linear model (MCGLM). For
the Gaussian graphical model (GGM) layer, \\\Sigma_b = (\sum_d \eta_d
Z_d)^{-1}\\ and \\\partial\Sigma_b/\partial\eta_d = -\Sigma_b Z_d
\Sigma_b\\; for the directed acyclic graph (DAG) layer, \\\Sigma_b =
T^{-1} D T^{-\top}\\, with \\T = I - \Gamma\\ and \\D =
\mathrm{diag}(\delta)\\, and the derivatives are taken with respect to
\\\gamma\\ and \\\delta\\, the innovation variances themselves and not
their logarithms. Matrices are returned in the order of the responses,
whatever the DAG ordering.

## Usage

``` r
mc_build_sigma_between_graph(theta_G, n_resp, spec, inverse = FALSE)
```

## Arguments

- theta_G:

  numeric vector of graph parameters. For `"ggm"`: the \\\eta\\ vector
  (diagonal parameters, then edges in canonical order, length \\R +
  \|E\|\\). For `"dag"`: `c(gamma, delta)`, with the \\\|E\|\\
  coefficients \\\gamma\\, one per edge, sorted lexicographically by the
  positions of its endpoints in the ordering (the canonical edge order
  under the natural ordering), followed by the \\R\\ innovation
  variances \\\delta\\, indexed by position in the ordering. The GGM
  vector must give a positive-definite precision, otherwise an error is
  raised; the DAG vector is expected to have \\\delta \> 0\\, which is
  not checked.

- n_resp:

  integer, number of responses \\R\\.

- spec:

  an object of class `"mcg_graph_spec"` from
  [`mc_between_graph`](https://evandeilton.github.io/mcglmgraph/reference/mc_between_graph.md).

- inverse:

  logical; if `TRUE`, also return the precision matrix.

## Value

A list with `Sigmab` (the \\R \times R\\ covariance), `D_Sigmab` (list
of derivatives, one per element of `theta_G`), `edge_pos` (positions of
the edge parameters in `theta_G`) and, when `inverse = TRUE`, `Kb` (the
precision matrix).

## See also

[`mc_between_graph`](https://evandeilton.github.io/mcglmgraph/reference/mc_between_graph.md),
[`mc_precision_bases`](https://evandeilton.github.io/mcglmgraph/reference/mc_precision_bases.md)

## Examples

``` r
sp <- mc_between_graph("ggm")
out <- mc_build_sigma_between_graph(c(2, 2, 2, -0.4, -0.3, 0), 3, sp)
out$Sigmab
#>           [,1]      [,2]  [,3]
#> [1,] 0.5333333 0.1066667 0.080
#> [2,] 0.1066667 0.5213333 0.016
#> [3,] 0.0800000 0.0160000 0.512
```
