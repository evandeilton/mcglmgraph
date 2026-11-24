# L1-penalized estimation of the precision graph

Fits the GGM layer with an L1 penalty on the edge parameters,
\$\$\hat{\boldsymbol\eta}^{(\lambda)} = \arg\min\_{\boldsymbol\eta}
\Bigl\\-\ell^p_e(\boldsymbol\eta) +
\lambda\sum\_{e}\|\eta_e\|\Bigr\\,\$\$ where \\\ell^p_e\\ is the
Gaussian working log-criterion of the residual layer (see
[`mc_graph_ploglik`](https://evandeilton.github.io/mcglmgraph/reference/mc_graph_ploglik.md)),
\\K_b=\sum_d\eta_d Z_d\\ is the linearly parametrized precision matrix,
and the sum runs over the edge parameters only: the diagonal parameters,
which set the scale of each response, are not penalized. The penalty is
on the accumulated scale of the criterion, so the smallest value that
zeroes every edge is \\\lambda\_{\max}=n\max\_{i\<j}\|\hat
S\_{e,ij}\|\\, with \\\hat S_e\\ the sample covariance of the Pearson
residuals: for \\\lambda\ge\lambda\_{\max}\\ the solution is the empty
graph. The problem is solved by the proximal gradient method (iterative
shrinkage-thresholding): a gradient step
\\\boldsymbol\eta+s\boldsymbol\psi\_{\boldsymbol\eta}\\, with
\\\boldsymbol\psi\_{\boldsymbol\eta}=\partial\ell^p_e/\partial\boldsymbol\eta\\
the Pearson estimating function, followed by soft-thresholding of the
edge coordinates at \\s\lambda\\. Backtracking halves the step until the
precision is positive definite and the sufficient-decrease condition of
Beck and Teboulle (2009) holds, under which the penalized objective does
not increase, so an admissible matrix alone is not sufficient for
accepting an iteration.

## Usage

``` r
mc_graph_penalized(
  fit,
  lambda,
  max_iter = 300L,
  tol = 1e-06,
  precondition = FALSE,
  kkt_tol = 1e-04,
  start = NULL
)
```

## Arguments

- fit:

  an object of class `"mcglm_graph"` (GGM layer).

- lambda:

  one finite non-negative penalty parameter, on the accumulated scale of
  the criterion (see the objective above); `lambda = 0` gives the
  unpenalized fit of the complete graph.

- max_iter:

  maximum number of proximal iterations (default 300).

- tol:

  tolerance on the relative change of the parameter vector,
  \\\max_d\|\eta_d^{+}-\eta_d\|/(1+\|\eta_d\|)\\ (default \\10^{-6}\\).

- precondition:

  logical; must be `FALSE`. `TRUE` is refused with an error, since the
  preconditioned step is not a valid proximal step. See the section
  "Preconditioned step".

- kkt_tol:

  tolerance for the per-unit KKT residual required to certify
  convergence, in addition to `tol` on the relative parameter change
  (default \\10^{-4}\\).

- start:

  optional starting value of the parameter vector, diagonal parameters
  first and then one per edge of the complete graph, as returned in
  `eta`; it must give a positive definite precision. The path of
  [`mc_select_lambda`](https://evandeilton.github.io/mcglmgraph/reference/mc_select_lambda.md)
  uses it for warm starts.

## Value

A list with `eta` (parameter vector, diagonal parameters first and then
one per edge of the complete graph), `edges` (two-column matrix of the
pairs whose coefficient is not exactly zero), `n_selected` (their
number), `Sigma_b` and `K_b` (fitted covariance and precision),
`iterations`, `converged` (logical), `status` (`"converged"`,
`"kkt_not_certified"`, `"maximum_iterations"` or
`"objective_line_search_failed"`), `objective` (final value of the
penalized objective), `objective_path` (its value along the iterations),
`kkt` (KKT residual, accumulated scale), `kkt_unit` (per unit, that is,
divided by \\n\\) and `kkt_tol` (the tolerance used).

## Convergence certificate

A small relative change of the parameter vector does not by itself
certify a stationary point of a non-smooth objective. Convergence is
therefore declared only when the relative change is below `tol` *and*
the per-unit KKT residual is below `kkt_tol`: when the first holds and
the second does not, the iteration goes on. The optimality conditions
are \\\hat\Sigma\_{b,kk}=\hat S\_{e,kk}\\ (diagonal), \\n(\hat
S\_{e,ij}-\hat\Sigma\_{b,ij})+\lambda\\\mathrm{sgn}(\hat K\_{b,ij})=0\\
for a nonzero edge and \\\|n(\hat
S\_{e,ij}-\hat\Sigma\_{b,ij})\|\le\lambda\\ for a zero edge. The KKT
residual \\r\_{\mathrm{KKT}}\\ is their largest violation over the
coordinates: \\\tfrac{n}{2}\|\hat S\_{e,kk}-\hat\Sigma\_{b,kk}\|\\ on
the diagonal, \\\|n(\hat
S\_{e,ij}-\hat\Sigma\_{b,ij})+\lambda\\\mathrm{sgn}(\hat K\_{b,ij})\|\\
on a nonzero edge and \\\max\\\|n(\hat
S\_{e,ij}-\hat\Sigma\_{b,ij})\|-\lambda,\\0\\\\ on a zero edge. The
certificate is \\r\_{\mathrm{KKT}}/n\le\\ `kkt_tol`. If `max_iter` is
reached in that state the status is `"kkt_not_certified"`, otherwise
`"maximum_iterations"`, and in both cases `converged` is `FALSE`.

## Step size

With \\f=-\ell^p_e\\ the smooth part of the objective, whose gradient is
\\-\boldsymbol\psi\_{\boldsymbol\eta}\\, every iteration starts from the
step \\s=\\\max_d J\_{n,dd}+10^{-8}\\^{-1}\\ and halves it until the
proximal point \\\boldsymbol\eta^{+}\\ gives a positive definite
precision and \$\$f(\boldsymbol\eta^{+})\le f(\boldsymbol\eta)-
\boldsymbol\psi\_{\boldsymbol\eta}^\top(\boldsymbol\eta^{+}-\boldsymbol\eta)
+\frac{1}{2s}\\\boldsymbol\eta^{+}-\boldsymbol\eta\\^2 .\$\$ The
starting step is a scale and not the inverse of a Lipschitz constant;
the condition is what guarantees the decrease.

## Starting value

The iteration starts from the diagonal fit that is compatible with the
empty graph, \\\eta\_{kk}=1/\widehat S\_{e,kk}\\ and every edge at zero,
which is also the fit from which the largest penalty of the path is
computed. Both need \\\widehat S\_{e,kk}\>0\\ for every response; a
response whose Pearson residuals have no variance is refused with an
error. The problem has a unique solution, so the start changes the
number of iterations and not the target; certified solutions from
different starts agree up to the tolerance.

## Candidate set

The path is always built over the complete graph on the \\R\\ responses,
whatever support the supplied fit carries: the penalty, not the input
support, decides which edges survive. Passing a fit that was restricted
to a sparse support therefore does not restrict the candidates, and the
returned `edges` may include pairs absent from that fit. To keep a
support fixed, set it through the `edges` argument of
[`mc_between_graph`](https://evandeilton.github.io/mcglmgraph/reference/mc_between_graph.md)
and refit, rather than penalizing.

## Preconditioned step

A preconditioned step, which multiplies the estimating function by
\\J^{-1}\\ and then applies coordinatewise soft-thresholding, is not the
proximal operator of the L1 penalty in the non-diagonal metric induced
by \\J\\: it can stop at a point that satisfies the parameter-change
test while grossly violating the KKT conditions. `precondition = TRUE`
therefore raises an error, and the fit always uses iterative
shrinkage-thresholding, whose certificate is the KKT residual above.

## References

Beck, A. and Teboulle, M. (2009). A fast iterative
shrinkage-thresholding algorithm for linear inverse problems. *SIAM
Journal on Imaging Sciences* 2(1), 183–202.
[doi:10.1137/080716542](https://doi.org/10.1137/080716542)

Friedman, J., Hastie, T. and Tibshirani, R. (2008). Sparse inverse
covariance estimation with the graphical lasso. *Biostatistics* 9(3),
432–441.
[doi:10.1093/biostatistics/kxm045](https://doi.org/10.1093/biostatistics/kxm045)

Yuan, M. and Lin, Y. (2007). Model selection and estimation in the
Gaussian graphical model. *Biometrika* 94(1), 19–35.
[doi:10.1093/biomet/asm018](https://doi.org/10.1093/biomet/asm018)

## See also

[`mc_select_lambda`](https://evandeilton.github.io/mcglmgraph/reference/mc_select_lambda.md)
to choose `lambda` by an information criterion.

Other structure selection:
[`mc_edge_scores()`](https://evandeilton.github.io/mcglmgraph/reference/mc_edge_scores.md),
[`mc_graph_stability()`](https://evandeilton.github.io/mcglmgraph/reference/mc_graph_stability.md),
[`mc_graph_stepwise()`](https://evandeilton.github.io/mcglmgraph/reference/mc_graph_stepwise.md),
[`mc_select_lambda()`](https://evandeilton.github.io/mcglmgraph/reference/mc_select_lambda.md)

## Examples

``` r
# \donttest{
set.seed(1); n <- 300; x <- rnorm(n)
y1 <- rnorm(n); y2 <- 0.6 * y1 + rnorm(n); y3 <- rnorm(n)
d <- data.frame(y1 = y1, y2 = y2, y3 = y3, x = x)
f <- mcglm_graph(list(y1 ~ x, y2 ~ x, y3 ~ x), data = d)
#> Automatic initial values selected.
pen <- mc_graph_penalized(f, lambda = 20)
pen$edges
#>      [,1] [,2]
#> [1,]    1    2
pen$status
#> [1] "converged"
# }
```
