# Covariance from the modified Cholesky (LDL) parametrization

Computes the between-response covariance \$\$\Sigma_b = T^{-1} D
(T^{-1})^\top,\$\$ where \\T\\ is the unit lower-triangular matrix of
structural coefficients, \\T\_{ji} = -\gamma\_{ji}\\, and \\D =
\mathrm{diag}(\delta)\\ (Pourahmadi, 1999). Because \\\det T = 1\\,
\\\Sigma_b\\ is positive definite for any `gamma` and any `delta > 0`,
so the parametrization is unconstrained apart from the positivity of
`delta`.

## Usage

``` r
mc_cholesky_to_covariance(gamma, delta, R, edges = NULL)
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

A list with `Sigma_b`, the \\R \times R\\ covariance (symmetrized),
`T_inv`, the inverse of \\T\\, and `T`.

## References

Pourahmadi, M. (1999). Joint mean-covariance models with applications to
longitudinal data: unconstrained parameterisation. *Biometrika* 86(3),
677-690.
[doi:10.1093/biomet/86.3.677](https://doi.org/10.1093/biomet/86.3.677)

## See also

Other directed layer:
[`mc_build_D_matrix()`](https://evandeilton.github.io/mcglmgraph/reference/mc_build_D_matrix.md),
[`mc_build_T_matrix()`](https://evandeilton.github.io/mcglmgraph/reference/mc_build_T_matrix.md),
[`mc_cholesky_derivatives()`](https://evandeilton.github.io/mcglmgraph/reference/mc_cholesky_derivatives.md),
[`mc_cholesky_to_precision()`](https://evandeilton.github.io/mcglmgraph/reference/mc_cholesky_to_precision.md),
[`mc_dag_closed_root()`](https://evandeilton.github.io/mcglmgraph/reference/mc_dag_closed_root.md),
[`mc_dag_coefficients()`](https://evandeilton.github.io/mcglmgraph/reference/mc_dag_coefficients.md),
[`mc_dag_sensitivity()`](https://evandeilton.github.io/mcglmgraph/reference/mc_dag_sensitivity.md),
[`mc_dag_unit_diagonal()`](https://evandeilton.github.io/mcglmgraph/reference/mc_dag_unit_diagonal.md)

## Examples

``` r
mc_cholesky_to_covariance(c(0.7, 0.3, 0.5), c(1, 0.64, 0.49), R = 3)$Sigma_b
#>      [,1]  [,2]   [,3]
#> [1,] 1.00 0.700 0.6500
#> [2,] 0.70 1.130 0.7750
#> [3,] 0.65 0.775 1.0725
```
