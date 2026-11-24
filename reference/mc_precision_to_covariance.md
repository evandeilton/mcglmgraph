# Convert a precision matrix into the between-response covariance

Inverts a symmetric positive-definite precision matrix, \\\Sigma_b =
K_b^{-1}\\, and reports its log-determinant and condition number. An
error is raised when `K` is not positive definite in the sense of
[`mc_check_spd`](https://evandeilton.github.io/mcglmgraph/reference/mc_check_spd.md),
that is, when its smallest eigenvalue does not exceed `1e-10`.

## Usage

``` r
mc_precision_to_covariance(K)
```

## Arguments

- K:

  symmetric positive-definite precision matrix.

## Value

A list with `Sigma_b` (the covariance \\K^{-1}\\), `logdet`
(log-determinant of \\K\\) and `cond` (condition number of \\K\\, the
ratio of its largest to its smallest eigenvalue).

## See also

Other GGM layer:
[`mc_build_precision()`](https://evandeilton.github.io/mcglmgraph/reference/mc_build_precision.md),
[`mc_check_spd()`](https://evandeilton.github.io/mcglmgraph/reference/mc_check_spd.md),
[`mc_precision_bases()`](https://evandeilton.github.io/mcglmgraph/reference/mc_precision_bases.md),
[`mc_precision_derivatives()`](https://evandeilton.github.io/mcglmgraph/reference/mc_precision_derivatives.md),
[`mc_spd_line_search()`](https://evandeilton.github.io/mcglmgraph/reference/mc_spd_line_search.md)

## Examples

``` r
K <- mc_build_precision(c(2, 2, 2, -0.4, -0.3, 0), mc_precision_bases(3))
mc_precision_to_covariance(K)$Sigma_b
#>           [,1]      [,2]  [,3]
#> [1,] 0.5333333 0.1066667 0.080
#> [2,] 0.1066667 0.5213333 0.016
#> [3,] 0.0800000 0.0160000 0.512
```
