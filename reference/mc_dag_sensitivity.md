# Block-diagonal sensitivity of the directed layer

Assembles the model-based Pearson sensitivity \\\mathbf J_n\\ of the
directed layer from its closed form. The empirical sandwich and
two-stage covariance functions use the observed Hessian of the working
criterion instead, since a sparse directed acyclic graph (DAG) can omit
second-moment structure. On the accumulated scale, with \\P_j\\ the
parents of response \\j\\ and \\\boldsymbol\gamma_j\\ their
coefficients, \$\$\mathbf J\_{n,\boldsymbol\gamma_j\boldsymbol\gamma_j}=
\frac{n}{\delta_j}\boldsymbol\Sigma\_{b,P_jP_j},\qquad
J\_{n,\delta_j\delta_j}=\frac{n}{2\delta_j^2},\$\$ and every other block
vanishes: between the coefficients of distinct responses, between
coefficients and innovation variances, and between innovation variances
of distinct responses. The generic assembly used by the fitting engine
forms \\P^2\\ traces of \\R \times R\\ products; this form costs one
block per node.

## Usage

``` r
mc_dag_sensitivity(Sigma_b, delta, n, R, edges = NULL)
```

## Arguments

- Sigma_b:

  the fitted \\R \times R\\ covariance of the layer.

- delta:

  the \\R\\ innovation variances.

- n:

  number of experimental units.

- R:

  a single positive integer, the number of responses.

- edges:

  optional integer matrix with columns \\(i, j)\\, \\i \< j\\, the
  directed edges \\i \to j\\. If `NULL`, the complete DAG on the natural
  ordering \\1 \< 2 \< \dots \< R\\ is used. Each edge must occur
  exactly once; duplicated rows, including reversed duplicates, are
  rejected because they do not define distinct structural coordinates. A
  row given as \\(j, i)\\ with \\j \> i\\ is read as the same edge \\i
  \to j\\: the direction always follows the natural order of the
  responses. The rows are sorted lexicographically before use, and
  `gamma` is read in that sorted order, not in the order the rows were
  supplied. Sort the matrix yourself, or read the order back from the
  `edges` component of the result, whenever the correspondence matters.

## Value

The \\(\|E\| + R) \times (\|E\| + R)\\ symmetric, block-diagonal
sensitivity in the parameter order `c(gamma, delta)`, where \\\|E\|\\ is
the number of edges and the `gamma` entries follow the sorted order of
the edges (see `edges` in
[`mc_build_T_matrix`](https://evandeilton.github.io/mcglmgraph/reference/mc_build_T_matrix.md)).

## See also

[`mc_dag_closed_root`](https://evandeilton.github.io/mcglmgraph/reference/mc_dag_closed_root.md)

Other directed layer:
[`mc_build_D_matrix()`](https://evandeilton.github.io/mcglmgraph/reference/mc_build_D_matrix.md),
[`mc_build_T_matrix()`](https://evandeilton.github.io/mcglmgraph/reference/mc_build_T_matrix.md),
[`mc_cholesky_derivatives()`](https://evandeilton.github.io/mcglmgraph/reference/mc_cholesky_derivatives.md),
[`mc_cholesky_to_covariance()`](https://evandeilton.github.io/mcglmgraph/reference/mc_cholesky_to_covariance.md),
[`mc_cholesky_to_precision()`](https://evandeilton.github.io/mcglmgraph/reference/mc_cholesky_to_precision.md),
[`mc_dag_closed_root()`](https://evandeilton.github.io/mcglmgraph/reference/mc_dag_closed_root.md),
[`mc_dag_coefficients()`](https://evandeilton.github.io/mcglmgraph/reference/mc_dag_coefficients.md),
[`mc_dag_unit_diagonal()`](https://evandeilton.github.io/mcglmgraph/reference/mc_dag_unit_diagonal.md)

## Examples

``` r
cc <- mc_cholesky_to_covariance(c(0.7, 0.3, 0.5), c(1, 0.64, 0.49), R = 3)
J <- mc_dag_sensitivity(cc$Sigma_b, c(1, 0.64, 0.49), n = 100, R = 3)
round(J, 6)
#>        [,1]     [,2]     [,3] [,4]     [,5]     [,6]
#> [1,] 156.25   0.0000   0.0000    0   0.0000   0.0000
#> [2,]   0.00 204.0816 142.8571    0   0.0000   0.0000
#> [3,]   0.00 142.8571 230.6122    0   0.0000   0.0000
#> [4,]   0.00   0.0000   0.0000   50   0.0000   0.0000
#> [5,]   0.00   0.0000   0.0000    0 122.0703   0.0000
#> [6,]   0.00   0.0000   0.0000    0   0.0000 208.2466
```
