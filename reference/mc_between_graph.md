# Specify the between-response graph layer

Builds the specification object that controls how the between-response
covariance \\\Sigma_b\\ of the Pearson residuals is parametrized in
[`mcglm_graph`](https://evandeilton.github.io/mcglmgraph/reference/mcglm_graph.md).
Two layers are available.

- `"ggm"`, the Gaussian graphical model (GGM): the precision matrix is
  linear in the parameters, \\K_b = \Sigma_b^{-1} = \sum_d \eta_d Z_d\\,
  with one diagonal parameter per response and one parameter per edge
  (see
  [`mc_precision_bases`](https://evandeilton.github.io/mcglmgraph/reference/mc_precision_bases.md)).
  A missing edge sets the corresponding precision entry, and hence the
  partial correlation, to zero.

- `"dag"`, the directed acyclic graph (DAG) layer, through the modified
  Cholesky (LDL) factorization \\\Sigma_b = T^{-1} D T^{-\top}\\, with
  \\T = I - \Gamma\\ unit lower triangular in the given ordering,
  \\\Gamma\\ holding the structural coefficients \\\gamma\_{ji}\\ on the
  edges, and \\D = \mathrm{diag}(\delta_1, \dots, \delta_R)\\ the
  innovation variances. Each response is regressed on its parents.

## Usage

``` r
mc_between_graph(type = c("ggm", "dag"), edges = NULL, ordering = NULL)
```

## Arguments

- type:

  character, the layer: `"ggm"` for the Gaussian graphical model
  (precision-based) or `"dag"` for the directed acyclic graph (modified
  Cholesky).

- edges:

  optional integer matrix with columns \\(i, j)\\ giving the edges as
  response positions; `NULL` means the complete graph. A matrix with no
  rows is the empty graph. The edges are stored in canonical form, each
  pair with \\i \< j\\ and the rows in lexicographic order, which is the
  order of the edge parameters, so the order in which they are supplied
  does not matter. In the DAG layer an edge is an adjacency and its
  direction is set by `ordering`.

- ordering:

  optional integer permutation of `1:R` defining the ordering of the DAG
  layer: an edge points from the response that comes first in `ordering`
  to the one that comes later. The default is the natural order. It has
  no meaning for `type = "ggm"`, whose precision support is unordered,
  and supplying it there is an error rather than a silently ignored
  argument.

## Value

An object of class `"mcg_graph_spec"`: a list with components `type`,
`edges` (canonical edge matrix or `NULL`) and `ordering` (integer vector
or `NULL`).

## References

Dempster, A. P. (1972). Covariance selection. *Biometrics*, 28(1),
157-175. [doi:10.2307/2528966](https://doi.org/10.2307/2528966)

Lauritzen, S. L. (1996). *Graphical Models*. Oxford University Press.

Pourahmadi, M. (1999). Joint mean-covariance models with applications to
longitudinal data: unconstrained parameterisation. *Biometrika*, 86(3),
677-690.
[doi:10.1093/biomet/86.3.677](https://doi.org/10.1093/biomet/86.3.677)

## See also

[`mcglm_graph`](https://evandeilton.github.io/mcglmgraph/reference/mcglm_graph.md),
[`mc_build_sigma_between_graph`](https://evandeilton.github.io/mcglmgraph/reference/mc_build_sigma_between_graph.md)

## Examples

``` r
mc_between_graph("ggm")
#> $type
#> [1] "ggm"
#> 
#> $edges
#> NULL
#> 
#> $ordering
#> NULL
#> 
#> attr(,"class")
#> [1] "mcg_graph_spec"
mc_between_graph("dag", edges = rbind(c(1, 2), c(2, 3)))
#> $type
#> [1] "dag"
#> 
#> $edges
#>      [,1] [,2]
#> [1,]    1    2
#> [2,]    2    3
#> 
#> $ordering
#> NULL
#> 
#> attr(,"class")
#> [1] "mcg_graph_spec"
```
