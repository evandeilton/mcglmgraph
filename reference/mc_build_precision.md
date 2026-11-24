# Assemble the between-response precision matrix

Computes \\K_b = \sum_d \eta_d Z_d\\ for a parameter vector \\\eta\\ and
the bases \\Z_d\\ of
[`mc_precision_bases`](https://evandeilton.github.io/mcglmgraph/reference/mc_precision_bases.md).
The result is symmetrized; positive definiteness is not checked here
(see
[`mc_check_spd`](https://evandeilton.github.io/mcglmgraph/reference/mc_check_spd.md)).

## Usage

``` r
mc_build_precision(eta, bases)
```

## Arguments

- eta:

  numeric vector of length \\q\_{\mathcal G} = R + \|E\|\\, the number
  of graph parameters: the diagonal parameters first, then the edge
  parameters, in the order of `bases`.

- bases:

  an object of class `"mcg_bases"` from
  [`mc_precision_bases`](https://evandeilton.github.io/mcglmgraph/reference/mc_precision_bases.md).

## Value

The \\R \times R\\ symmetric precision matrix \\K_b\\.

## See also

Other GGM layer:
[`mc_check_spd()`](https://evandeilton.github.io/mcglmgraph/reference/mc_check_spd.md),
[`mc_precision_bases()`](https://evandeilton.github.io/mcglmgraph/reference/mc_precision_bases.md),
[`mc_precision_derivatives()`](https://evandeilton.github.io/mcglmgraph/reference/mc_precision_derivatives.md),
[`mc_precision_to_covariance()`](https://evandeilton.github.io/mcglmgraph/reference/mc_precision_to_covariance.md),
[`mc_spd_line_search()`](https://evandeilton.github.io/mcglmgraph/reference/mc_spd_line_search.md)

## Examples

``` r
b <- mc_precision_bases(3)
mc_build_precision(c(2, 2, 2, -0.4, -0.3, 0), b)
#>      [,1] [,2] [,3]
#> [1,]  2.0 -0.4 -0.3
#> [2,] -0.4  2.0  0.0
#> [3,] -0.3  0.0  2.0
```
