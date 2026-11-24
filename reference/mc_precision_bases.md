# Build the base matrices of the GGM precision predictor

Constructs the base matrices \\Z_d\\ of the linear-in-precision
predictor of the Gaussian graphical model (GGM) layer, \$\$K_b =
\sum\_{d=1}^{D} \eta_d Z_d,\qquad \Sigma_b = K_b^{-1},\$\$ where \\K_b\\
is the between-response precision matrix and \\D = R + \|E\|\\. There
are \\R\\ diagonal bases \\Z\_{kk} = e_k e_k^\top\\ and one symmetric
base \\Z\_{ij} = e_i e_j^\top + e_j e_i^\top\\ per edge \\(i, j)\\, so
that \\\eta\_{kk} = \[K_b\]\_{kk}\\ and \\\eta\_{ij} = \[K_b\]\_{ij}\\.
The pair \\(i, j)\\ is conditionally independent given the other
responses, in the Gaussian sense, exactly when it is outside the edge
set (Dempster, 1972; Lauritzen, 1996).

## Usage

``` r
mc_precision_bases(R, edges = NULL)
```

## Arguments

- R:

  integer, number of responses (\\R \ge 2\\).

- edges:

  optional integer matrix with two columns giving the edges \\(i, j)\\,
  \\i \< j\\, as response positions in `1:R`. If `NULL`, the complete
  graph is used; a matrix with no rows is the empty graph. The edges are
  stored sorted in lexicographic order, which is the order of the edge
  parameters.

## Value

An object of class `"mcg_bases"`: a list with components `bases` (list
of \\D\\ symmetric \\R \times R\\ matrices), `cube` (the same bases as
an \\R \times R \times D\\ array, the form consumed by the C++ kernels),
`names` (parameter names `eta_kk` and `eta_ij`), `type` (`"diag"` or
`"edge"` for each base), `edges` (the canonical edge matrix) and `R`.

## References

Dempster, A. P. (1972). Covariance selection. *Biometrics*, 28(1),
157-175. [doi:10.2307/2528966](https://doi.org/10.2307/2528966)

Lauritzen, S. L. (1996). *Graphical Models*. Oxford University Press.

## See also

[`mc_build_precision`](https://evandeilton.github.io/mcglmgraph/reference/mc_build_precision.md),
[`mc_between_graph`](https://evandeilton.github.io/mcglmgraph/reference/mc_between_graph.md)

Other GGM layer:
[`mc_build_precision()`](https://evandeilton.github.io/mcglmgraph/reference/mc_build_precision.md),
[`mc_check_spd()`](https://evandeilton.github.io/mcglmgraph/reference/mc_check_spd.md),
[`mc_precision_derivatives()`](https://evandeilton.github.io/mcglmgraph/reference/mc_precision_derivatives.md),
[`mc_precision_to_covariance()`](https://evandeilton.github.io/mcglmgraph/reference/mc_precision_to_covariance.md),
[`mc_spd_line_search()`](https://evandeilton.github.io/mcglmgraph/reference/mc_spd_line_search.md)

## Examples

``` r
b <- mc_precision_bases(3)
length(b$bases)            # 3 diagonal + 3 edges = 6
#> [1] 6
mc_precision_bases(3, edges = rbind(c(1, 2), c(2, 3)))$names
#> [1] "eta_11" "eta_22" "eta_33" "eta_12" "eta_23"
```
