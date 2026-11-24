# Check positive-definiteness of a matrix

A symmetric matrix is declared positive definite when the smallest
eigenvalue of its symmetric part \\(K + K^\top)/2\\ exceeds `tol`.

## Usage

``` r
mc_check_spd(K, tol = 1e-10)
```

## Arguments

- K:

  square numeric matrix; only its symmetric part is used.

- tol:

  tolerance on the smallest eigenvalue (default `1e-10`).

## Value

A list with `is_spd` (logical, `TRUE` when the smallest eigenvalue
exceeds `tol`) and `min_eig` (the smallest eigenvalue).

## See also

Other GGM layer:
[`mc_build_precision()`](https://evandeilton.github.io/mcglmgraph/reference/mc_build_precision.md),
[`mc_precision_bases()`](https://evandeilton.github.io/mcglmgraph/reference/mc_precision_bases.md),
[`mc_precision_derivatives()`](https://evandeilton.github.io/mcglmgraph/reference/mc_precision_derivatives.md),
[`mc_precision_to_covariance()`](https://evandeilton.github.io/mcglmgraph/reference/mc_precision_to_covariance.md),
[`mc_spd_line_search()`](https://evandeilton.github.io/mcglmgraph/reference/mc_spd_line_search.md)

## Examples

``` r
mc_check_spd(diag(3))
#> $is_spd
#> [1] TRUE
#> 
#> $min_eig
#> [1] 1
#> 
```
