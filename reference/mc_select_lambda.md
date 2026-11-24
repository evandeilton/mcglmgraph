# Select the penalty parameter along an L1 path

Builds a geometric grid of penalty values down from the smallest value
that zeroes every edge, \\\lambda\_{\max}=n\max\_{i\<j}\|\hat
S\_{e,ij}\|\\, to `lambda_min_ratio` times that value. At each grid
point the penalized model of
[`mc_graph_penalized`](https://evandeilton.github.io/mcglmgraph/reference/mc_graph_penalized.md)
is fitted, its support is refitted without penalization, and the
information criterion is evaluated on that refit. The selected penalty
is \$\$\hat\lambda=\arg\min\_{\lambda}\\
\mathrm{crit}\bigl(\hat{\mathcal G}\_\lambda\bigr),\$\$ with
\\\hat{\mathcal G}\_\lambda\\ the support of the penalized solution at
\\\lambda\\; ties go to the larger penalty. The grid starts at
\\\max(\lambda\_{\max},10^{-3})\\: the lower limit keeps the geometric
grid defined when the residual covariances are numerically zero, and it
does not change the first solution, which is the empty graph at any
penalty of at least \\\lambda\_{\max}\\.

## Usage

``` r
mc_select_lambda(
  fit,
  n_lambda = 30L,
  lambda_min_ratio = 0.001,
  criterion = c("bic", "aic", "qic", "gbic", "clbic", "ebic"),
  refit = TRUE,
  workers = 1L,
  warm_start = FALSE,
  include_zero = FALSE,
  ebic_gamma = 0.5,
  cluster = NULL,
  small_sample = FALSE,
  max_iter = 300L,
  tol = 1e-06,
  kkt_tol = 1e-04
)
```

## Arguments

- fit:

  an object of class `"mcglm_graph"` (GGM layer).

- n_lambda:

  number of grid points (default 30).

- lambda_min_ratio:

  ratio of the smallest to the largest penalty (default \\10^{-3}\\).

- criterion:

  `"bic"` (default), `"aic"`, `"qic"`, `"gbic"` or its alias `"clbic"`,
  or `"ebic"`.

- refit:

  logical; if `TRUE` (default) the returned fit is the unpenalized refit
  of the selected support. If `FALSE`, the penalized solution is
  attached to the fit instead; its covariance is kept for diagnostics
  but it is flagged as not certified for edge-wise inference.

- workers:

  number of parallel workers for the penalty grid. The grid points are
  independent fits, so the selected penalty does not depend on this
  value. Sequential by default; see
  [`mc_graph_stability`](https://evandeilton.github.io/mcglmgraph/reference/mc_graph_stability.md)
  for why it is opt-in. Not used with `warm_start = TRUE`.

- warm_start:

  logical; run the grid in sequence from the largest penalty, each fit
  starting at the previous solution.

- include_zero:

  logical; add the unpenalized complete graph as the candidate at
  \\\lambda=0\\.

- ebic_gamma:

  the parameter of the extended BIC, used only with
  `criterion = "ebic"`.

- cluster:

  optional independent-cluster identifier, one per residual row, passed
  to
  [`mc_compute_criterion`](https://evandeilton.github.io/mcglmgraph/reference/mc_compute_criterion.md):
  Pearson contributions are summed within clusters in QIC and GBIC, and
  the number of clusters replaces the number of rows in the logarithmic
  penalty of BIC, GBIC and EBIC.

- small_sample:

  logical; apply the cluster factor \\G/(G-1)\\, with \\G\\ the number
  of clusters, to the empirical variability of the clustered criteria
  (default `FALSE`); ignored without `cluster`.

- max_iter, tol, kkt_tol:

  passed to
  [`mc_graph_penalized`](https://evandeilton.github.io/mcglmgraph/reference/mc_graph_penalized.md)
  at every grid point.

## Value

An object of class `"mcg_selection"` (see
[`mc_graph_stepwise`](https://evandeilton.github.io/mcglmgraph/reference/mc_graph_stepwise.md)
for the common components) with `method = "penalized"`, the selected
`lambda`, a logical `warm_start`, and the table `lambda_path`, a data
frame with one row per grid point (plus one for \\\lambda=0\\ when
`include_zero = TRUE`) and columns `lambda`, `n_edges`, `criterion`,
`penalized_converged`, `penalized_status`, `refit_converged` and
`refit_status`.

## Independent fits, warm starts and the unpenalized end

By default every grid point is an independent fit from the same starting
value, the diagonal fit \\\eta\_{kk}=1/\widehat S\_{e,kk}\\, which also
gives the largest penalty of the path; both need \\\widehat
S\_{e,kk}\>0\\ for every response, and a response whose Pearson
residuals have no variance is refused with an error. The points can be
evaluated in any order and with any number of workers, and the selected
penalty does not depend on either. With `warm_start = TRUE` the grid is
run from the largest penalty down, each fit starting at the solution of
the one before, which is the usual way of computing a path; the fits are
then sequential and `workers` is not used. The problem has a unique
solution, so the two routes share the target and differ in the number of
iterations. Certified solutions satisfy the optimality conditions only
up to the tolerance, so their supports can differ in pairs whose
zero-edge optimality condition, the scaled gap between the residual
covariance and the fitted covariance, is close to \\\lambda\\, within a
margin proportional to the tolerance and to the conditioning of the
curvature. With `include_zero = TRUE` the unpenalized fit of the
complete graph is added to the candidates as the end of the path,
\\\lambda=0\\; the geometric grid never reaches it. By default the grid
is independent and does not include that end.

## Which grid points compete

The criterion is evaluated on the unpenalized refit of the support of
every grid point whose penalized fit was certified, by the parameter
change and by the KKT residual, and whose refit converged. A point
without the certificate stays in `lambda_path` with its status and does
not enter the comparison. An error is raised when no grid point
qualifies.

## See also

[`mc_graph_penalized`](https://evandeilton.github.io/mcglmgraph/reference/mc_graph_penalized.md),
[`mc_compute_criterion`](https://evandeilton.github.io/mcglmgraph/reference/mc_compute_criterion.md).

Other structure selection:
[`mc_edge_scores()`](https://evandeilton.github.io/mcglmgraph/reference/mc_edge_scores.md),
[`mc_graph_penalized()`](https://evandeilton.github.io/mcglmgraph/reference/mc_graph_penalized.md),
[`mc_graph_stability()`](https://evandeilton.github.io/mcglmgraph/reference/mc_graph_stability.md),
[`mc_graph_stepwise()`](https://evandeilton.github.io/mcglmgraph/reference/mc_graph_stepwise.md)

## Examples

``` r
# \donttest{
set.seed(1); n <- 300; x <- rnorm(n)
y1 <- rnorm(n); y2 <- 0.6 * y1 + rnorm(n); y3 <- rnorm(n)
d <- data.frame(y1 = y1, y2 = y2, y3 = y3, x = x)
f <- mcglm_graph(list(y1 ~ x, y2 ~ x, y3 ~ x), data = d)
#> Automatic initial values selected.
sel <- mc_select_lambda(f, n_lambda = 10)
sel$edges
#>      [,1] [,2]
#> [1,]    1    2
sel$lambda_path[, c("lambda", "n_edges", "criterion", "penalized_status")]
#>         lambda n_edges criterion penalized_status
#> 1  156.2001035       0  2565.181        converged
#> 2   72.5016656       1  2474.515        converged
#> 3   33.6522922       1  2474.515        converged
#> 4   15.6200104       1  2474.515        converged
#> 5    7.2501666       3  2484.557        converged
#> 6    3.3652292       3  2484.557        converged
#> 7    1.5620010       3  2484.557        converged
#> 8    0.7250167       3  2484.557        converged
#> 9    0.3365229       3  2484.557        converged
#> 10   0.1562001       3  2484.557        converged
# }
```
