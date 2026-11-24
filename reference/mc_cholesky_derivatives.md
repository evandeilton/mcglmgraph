# Analytical LDL derivatives of the covariance

Returns the derivatives of \\\Sigma_b = T^{-1} D (T^{-1})^\top\\ with
respect to the innovation variances \\\delta_k\\ and the structural
coefficients \\\gamma\_{ji}\\:
\$\$\frac{\partial\Sigma_b}{\partial\delta_k}=
T^{-1}E\_{kk}(T^{-1})^\top,\qquad
\frac{\partial\Sigma_b}{\partial\gamma\_{ji}}=
T^{-1}E\_{ji}\Sigma_b+\Sigma_b E\_{ij}(T^{-1})^\top,\$\$ where
\\E\_{ab}\\ is the matrix with a one in position \\(a, b)\\ and zeros
elsewhere. The second form is symmetric, its second term being the
transpose of the first. These are the matrices that enter the estimating
functions of the directed layer.

## Usage

``` r
mc_cholesky_derivatives(gamma, delta, R, edges = NULL)
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

A list with `D_delta`, a list of \\R\\ symmetric matrices (one per
innovation variance, in the order of `delta`), and `D_gamma`, a list
with one symmetric matrix per edge, in the sorted order of the edges
(see `edges` in
[`mc_build_T_matrix`](https://evandeilton.github.io/mcglmgraph/reference/mc_build_T_matrix.md)).

## See also

Other directed layer:
[`mc_build_D_matrix()`](https://evandeilton.github.io/mcglmgraph/reference/mc_build_D_matrix.md),
[`mc_build_T_matrix()`](https://evandeilton.github.io/mcglmgraph/reference/mc_build_T_matrix.md),
[`mc_cholesky_to_covariance()`](https://evandeilton.github.io/mcglmgraph/reference/mc_cholesky_to_covariance.md),
[`mc_cholesky_to_precision()`](https://evandeilton.github.io/mcglmgraph/reference/mc_cholesky_to_precision.md),
[`mc_dag_closed_root()`](https://evandeilton.github.io/mcglmgraph/reference/mc_dag_closed_root.md),
[`mc_dag_coefficients()`](https://evandeilton.github.io/mcglmgraph/reference/mc_dag_coefficients.md),
[`mc_dag_sensitivity()`](https://evandeilton.github.io/mcglmgraph/reference/mc_dag_sensitivity.md),
[`mc_dag_unit_diagonal()`](https://evandeilton.github.io/mcglmgraph/reference/mc_dag_unit_diagonal.md)

## Examples

``` r
d <- mc_cholesky_derivatives(c(0.7, 0.3, 0.5), c(1, 0.64, 0.49), R = 3)
length(d$D_gamma)
#> [1] 3
```
