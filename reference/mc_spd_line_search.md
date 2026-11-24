# Backtracking line search preserving positive-definiteness

Given a Newton direction \\d\\, finds the largest step \$\$\alpha =
\max\\2^{-k}: k = 0, 1, \dots, k\_{\max},\\ \lambda\_{\min}(K_b(\eta +
2^{-k} d)) \> 10^{-10}\\,\$\$ with \\\lambda\_{\min}\\ the smallest
eigenvalue and \\k\_{\max}\\ given by `max_halvings`. The positivity
test uses the numerical threshold \\10^{-10}\\ in place of zero. Only
positive definiteness is checked, not a decrease of a merit function. By
Weyl's inequality, \\\lambda\_{\min}(K_b + \alpha\\\Delta K) \ge
\lambda\_{\min}(K_b) - \alpha\\\\\Delta K\\\_2\\, with \\\Delta K =
\sum_d d_d Z_d\\, so an admissible step exists for small enough
\\\alpha\\ whenever \\K_b(\eta)\\ passes the test; it may still lie
below \\2^{-k\_{\max}}\\, in which case the search fails. The fitting
engine of
[`mcglm_graph`](https://evandeilton.github.io/mcglmgraph/reference/mcglm_graph.md)
calls the same search with 30 halvings.

## Usage

``` r
mc_spd_line_search(eta_old, delta_eta, bases, max_halvings = 30L)
```

## Arguments

- eta_old:

  current parameter vector \\\eta\\, with \\K_b(\eta)\\ positive
  definite.

- delta_eta:

  the Newton direction \\d\\, of the same length.

- bases:

  an object of class `"mcg_bases"` from
  [`mc_precision_bases`](https://evandeilton.github.io/mcglmgraph/reference/mc_precision_bases.md).

- max_halvings:

  maximum number \\k\_{\max}\\ of step halvings (default 30), so that at
  most `max_halvings + 1` steps are tried.

## Value

A list with components

- `eta_new`:

  a one-column matrix with the last point tried, \\\eta + \alpha d\\:
  the accepted point when the search succeeds.

- `step_size`:

  the accepted step \\\alpha\\ when the search succeeds.

- `halvings`:

  the number of halvings performed.

- `is_spd`:

  logical, `TRUE` when an admissible step was found.

When no admissible step is found, `is_spd` is `FALSE`, `eta_new` is the
point tried with the smallest step, \\2^{-k\_{\max}}\\, which fails the
test and must not be used, `halvings` is `max_halvings + 1` and
`step_size` is \\2^{-(k\_{\max}+1)}\\, half the step that produced
`eta_new`.

## References

Horn, R. A. and Johnson, C. R. (2013). *Matrix Analysis*, 2nd ed.
Cambridge University Press.
[doi:10.1017/CBO9781139020411](https://doi.org/10.1017/CBO9781139020411)

## See also

[`mc_check_spd`](https://evandeilton.github.io/mcglmgraph/reference/mc_check_spd.md)

Other GGM layer:
[`mc_build_precision()`](https://evandeilton.github.io/mcglmgraph/reference/mc_build_precision.md),
[`mc_check_spd()`](https://evandeilton.github.io/mcglmgraph/reference/mc_check_spd.md),
[`mc_precision_bases()`](https://evandeilton.github.io/mcglmgraph/reference/mc_precision_bases.md),
[`mc_precision_derivatives()`](https://evandeilton.github.io/mcglmgraph/reference/mc_precision_derivatives.md),
[`mc_precision_to_covariance()`](https://evandeilton.github.io/mcglmgraph/reference/mc_precision_to_covariance.md)

## Examples

``` r
b <- mc_precision_bases(3)
mc_spd_line_search(c(2,2,2,0,0,0), c(-.5,-.5,-.5,-.4,-.4,-.4), b)$halvings
#> [1] 0
```
