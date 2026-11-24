# Build the diagonal innovation-variance matrix D

Returns the diagonal matrix \\D = \mathrm{diag}(\delta_1, \dots,
\delta_R)\\ of innovation variances of the LDL parametrization,
\\\Sigma_b = T^{-1} D (T^{-1})^\top\\. The innovation variance
\\\delta_j\\ is the variance of response \\j\\ that remains after
regression on its parents.

## Usage

``` r
mc_build_D_matrix(delta, R)
```

## Arguments

- delta:

  numeric vector of exactly \\R\\ finite, strictly positive innovation
  variances.

- R:

  a single positive integer, the number of responses.

## Value

The \\R \times R\\ diagonal matrix \\D\\.

## See also

[`mc_build_T_matrix`](https://evandeilton.github.io/mcglmgraph/reference/mc_build_T_matrix.md)

Other directed layer:
[`mc_build_T_matrix()`](https://evandeilton.github.io/mcglmgraph/reference/mc_build_T_matrix.md),
[`mc_cholesky_derivatives()`](https://evandeilton.github.io/mcglmgraph/reference/mc_cholesky_derivatives.md),
[`mc_cholesky_to_covariance()`](https://evandeilton.github.io/mcglmgraph/reference/mc_cholesky_to_covariance.md),
[`mc_cholesky_to_precision()`](https://evandeilton.github.io/mcglmgraph/reference/mc_cholesky_to_precision.md),
[`mc_dag_closed_root()`](https://evandeilton.github.io/mcglmgraph/reference/mc_dag_closed_root.md),
[`mc_dag_coefficients()`](https://evandeilton.github.io/mcglmgraph/reference/mc_dag_coefficients.md),
[`mc_dag_sensitivity()`](https://evandeilton.github.io/mcglmgraph/reference/mc_dag_sensitivity.md),
[`mc_dag_unit_diagonal()`](https://evandeilton.github.io/mcglmgraph/reference/mc_dag_unit_diagonal.md)

## Examples

``` r
mc_build_D_matrix(c(1, 0.8, 0.6), 3)
#>      [,1] [,2] [,3]
#> [1,]    1  0.0  0.0
#> [2,]    0  0.8  0.0
#> [3,]    0  0.0  0.6
```
