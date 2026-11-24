# Analytical derivatives of the covariance with respect to eta

Returns the derivatives of \\\Sigma_b = K_b^{-1}\\ with respect to the
precision parameters, \$\$\frac{\partial \Sigma_b}{\partial \eta_d} =
-\Sigma_b Z_d \Sigma_b, \qquad d = 1, \dots, D,\$\$ which follow from
\\\partial K_b^{-1} = -K_b^{-1}(\partial K_b)K_b^{-1}\\ and \\\partial
K_b/\partial\eta_d = Z_d\\. Each matrix is symmetrized.

## Usage

``` r
mc_precision_derivatives(Sigma_b, bases)
```

## Arguments

- Sigma_b:

  the between-response covariance matrix \\\Sigma_b\\.

- bases:

  an object of class `"mcg_bases"`.

## Value

A list of \\D\\ symmetric \\R \times R\\ matrices in the order of
`bases`.

## See also

Other GGM layer:
[`mc_build_precision()`](https://evandeilton.github.io/mcglmgraph/reference/mc_build_precision.md),
[`mc_check_spd()`](https://evandeilton.github.io/mcglmgraph/reference/mc_check_spd.md),
[`mc_precision_bases()`](https://evandeilton.github.io/mcglmgraph/reference/mc_precision_bases.md),
[`mc_precision_to_covariance()`](https://evandeilton.github.io/mcglmgraph/reference/mc_precision_to_covariance.md),
[`mc_spd_line_search()`](https://evandeilton.github.io/mcglmgraph/reference/mc_spd_line_search.md)

## Examples

``` r
b <- mc_precision_bases(3)
S <- mc_precision_to_covariance(mc_build_precision(c(2,2,2,-.4,-.3,0), b))$Sigma_b
length(mc_precision_derivatives(S, b))
#> [1] 6
```
