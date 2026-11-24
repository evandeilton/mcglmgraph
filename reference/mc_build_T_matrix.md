# Build the structural matrix T of the LDL parametrization

Assembles the unit lower-triangular matrix \\T\\ of the modified
Cholesky (LDL) parametrization of a directed acyclic graph (DAG) on
\\R\\ responses, in which the between-response covariance is \\\Sigma_b
= T^{-1} D (T^{-1})^\top\\ with \\D\\ diagonal and positive (Pourahmadi,
1999). Each edge \\i \to j\\ contributes one structural coefficient
\\\gamma\_{ji}\\, the regression coefficient of response \\j\\ on its
parent \\i\\, and the responses are listed in a topological order of the
graph (parents before children).

## Usage

``` r
mc_build_T_matrix(gamma, R, edges = NULL)
```

## Arguments

- gamma:

  finite numeric vector of structural coefficients, exactly one per
  edge.

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

The \\R \times R\\ unit lower-triangular matrix \\T\\ with \\T\_{ji} =
-\gamma\_{ji}\\ for each edge \\i \to j\\ and zero elsewhere off the
diagonal. Its determinant is one.

## References

Pourahmadi, M. (1999). Joint mean-covariance models with applications to
longitudinal data: unconstrained parameterisation. *Biometrika* 86(3),
677-690.
[doi:10.1093/biomet/86.3.677](https://doi.org/10.1093/biomet/86.3.677)

## See also

[`mc_build_D_matrix`](https://evandeilton.github.io/mcglmgraph/reference/mc_build_D_matrix.md),
[`mc_cholesky_to_covariance`](https://evandeilton.github.io/mcglmgraph/reference/mc_cholesky_to_covariance.md)

Other directed layer:
[`mc_build_D_matrix()`](https://evandeilton.github.io/mcglmgraph/reference/mc_build_D_matrix.md),
[`mc_cholesky_derivatives()`](https://evandeilton.github.io/mcglmgraph/reference/mc_cholesky_derivatives.md),
[`mc_cholesky_to_covariance()`](https://evandeilton.github.io/mcglmgraph/reference/mc_cholesky_to_covariance.md),
[`mc_cholesky_to_precision()`](https://evandeilton.github.io/mcglmgraph/reference/mc_cholesky_to_precision.md),
[`mc_dag_closed_root()`](https://evandeilton.github.io/mcglmgraph/reference/mc_dag_closed_root.md),
[`mc_dag_coefficients()`](https://evandeilton.github.io/mcglmgraph/reference/mc_dag_coefficients.md),
[`mc_dag_sensitivity()`](https://evandeilton.github.io/mcglmgraph/reference/mc_dag_sensitivity.md),
[`mc_dag_unit_diagonal()`](https://evandeilton.github.io/mcglmgraph/reference/mc_dag_unit_diagonal.md)

## Examples

``` r
mc_build_T_matrix(c(0.4, 0.2, 0.3), R = 3)
#>      [,1] [,2] [,3]
#> [1,]  1.0  0.0    0
#> [2,] -0.4  1.0    0
#> [3,] -0.2 -0.3    1
```
