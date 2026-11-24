# Precision from the modified Cholesky (LDL) parametrization

Computes the between-response precision matrix \\K_b = \Sigma_b^{-1} =
T^\top D^{-1} T\\ of the LDL parametrization directly from `gamma` and
`delta`, without inverting \\\Sigma_b\\. A zero entry of \\K_b\\ is a
missing edge of the concentration graph. For \\i \< k\\,
\$\$\[K_b\]\_{ik}=-\frac{\gamma\_{ki}}{\delta_k}+
\sum\_{j\>k}\frac{\gamma\_{ji}\gamma\_{jk}}{\delta_j},\$\$ so the
concentration graph is contained in the moral graph of the directed
acyclic graph (DAG), in which the parents of a common child are joined,
and coincides with it except where the two terms cancel, a set of
coefficients of Lebesgue measure zero.

## Usage

``` r
mc_cholesky_to_precision(gamma, delta, R, edges = NULL)
```

## Arguments

- gamma:

  finite numeric vector of structural coefficients, exactly one per
  edge.

- delta:

  numeric vector of exactly \\R\\ finite, strictly positive innovation
  variances.

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

The \\R \times R\\ precision matrix \\K_b = T^\top D^{-1} T\\
(symmetrized).

## See also

Other directed layer:
[`mc_build_D_matrix()`](https://evandeilton.github.io/mcglmgraph/reference/mc_build_D_matrix.md),
[`mc_build_T_matrix()`](https://evandeilton.github.io/mcglmgraph/reference/mc_build_T_matrix.md),
[`mc_cholesky_derivatives()`](https://evandeilton.github.io/mcglmgraph/reference/mc_cholesky_derivatives.md),
[`mc_cholesky_to_covariance()`](https://evandeilton.github.io/mcglmgraph/reference/mc_cholesky_to_covariance.md),
[`mc_dag_closed_root()`](https://evandeilton.github.io/mcglmgraph/reference/mc_dag_closed_root.md),
[`mc_dag_coefficients()`](https://evandeilton.github.io/mcglmgraph/reference/mc_dag_coefficients.md),
[`mc_dag_sensitivity()`](https://evandeilton.github.io/mcglmgraph/reference/mc_dag_sensitivity.md),
[`mc_dag_unit_diagonal()`](https://evandeilton.github.io/mcglmgraph/reference/mc_dag_unit_diagonal.md)

## Examples

``` r
mc_cholesky_to_precision(c(0.7, 0.3, 0.5), c(1, 0.64, 0.49), R = 3)
#>            [,1]       [,2]       [,3]
#> [1,]  1.9492985 -0.7876276 -0.6122449
#> [2,] -0.7876276  2.0727041 -1.0204082
#> [3,] -0.6122449 -1.0204082  2.0408163
```
