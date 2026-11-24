# Fit a graph-structured multivariate covariance generalized linear model

Fits a multivariate covariance generalized linear model (MCGLM) and then
estimates a graph-structured covariance for its Pearson residuals, in a
two-stage residual analysis. The between-response covariance is
parametrized either by a Gaussian graphical model (GGM, linear in the
precision matrix) or by a directed acyclic graph (DAG, modified Cholesky
factorization).

## Usage

``` r
mcglm_graph(
  linear_pred,
  between_graph = mc_between_graph("ggm"),
  data,
  matrix_pred = NULL,
  control = list(),
  ...
)

mcglmgraph(
  linear_pred,
  between_graph = mc_between_graph("ggm"),
  data,
  matrix_pred = NULL,
  control = list(),
  ...
)

# S3 method for class 'mcglm_graph'
print(x, ...)

# S3 method for class 'mcglm_graph'
summary(
  object,
  variability = c("model", "empirical", "two_stage"),
  cluster = NULL,
  small_sample = FALSE,
  ...
)

# S3 method for class 'summary.mcglm_graph'
print(x, ...)

# S3 method for class 'mcglm_graph'
coef(object, type = NULL, ...)

# S3 method for class 'mcglm_graph'
plot(x, type = c("network", "residuals"), ...)
```

## Arguments

- linear_pred:

  a list of model formulas, one per response; at least two responses are
  required.

- between_graph:

  a specification from
  [`mc_between_graph`](https://evandeilton.github.io/mcglmgraph/reference/mc_between_graph.md)
  (default: the GGM layer on the complete graph).

- data:

  a data frame containing the variables in `linear_pred`. It must hold
  complete cases: a missing or non-finite response value, a missing
  covariate, or a non-finite `offset`, `Ntrial` or `weights` entry
  passed through `...`, is refused with an error that names the response
  and the number of rows affected, because the graph layer needs one
  complete vector of Pearson residuals per unit. Removing or imputing
  incomplete rows is left to the user.

- matrix_pred:

  optional list of within-response covariance bases, as in
  [`mcglm`](https://rdrr.io/pkg/mcglm/man/mcglm.html), one element per
  response. When `NULL`, an identity structure
  ([`mc_id`](https://rdrr.io/pkg/mcglm/man/mc_id.html)) is used for
  every response.

- control:

  a named list of control parameters for the between-response engine;
  unknown names are an error. The components are:

  `max_iter`

  :   maximum number of Newton iterations (default 100).

  `tol`

  :   tolerance on the relative parameter change \\\max_d
      \|\theta_d^{(t+1)} - \theta_d^{(t)}\| / (1 + \|\theta_d^{(t)}\|)\\
      between consecutive iterations, computed on \\\eta\\ or on
      \\(\gamma, \delta)\\ (default `1e-8`).

  `score_tol`

  :   bound on the maximum absolute value of the estimating function
      divided by \\n\\ at the solution (default `sqrt(tol)`).

  `rcond_tol`

  :   smallest admissible reciprocal condition number of the sensitivity
      matrix; below it a ridge is used in linear solves and the fit is
      not certified (default `.Machine$double.eps^(2/3)`).

  `verbose`

  :   logical; report the iterations (default `FALSE`).

  `dag_solver`

  :   how the DAG layer is solved: `"closed"` (default) takes the
      node-by-node root of
      [`mc_dag_closed_root`](https://evandeilton.github.io/mcglmgraph/reference/mc_dag_closed_root.md),
      and `"newton"` runs the scoring iteration shared with the GGM
      layer. The two solve the same equations, the first exactly and
      without iterating, and agree up to the convergence tolerance of
      the iteration; the closed route is much faster when the number of
      responses is large. It has no effect on the GGM layer, which has
      no closed root.

- ...:

  for `mcglm_graph()` and `mcglmgraph()`: further arguments passed to
  [`mcglm`](https://rdrr.io/pkg/mcglm/man/mcglm.html) (such as `link`,
  `variance`, `power_fixed`, `offset`, `Ntrial`, `weights`). For the
  methods: ignored by the `print` methods and by `plot` with
  `type = "residuals"`; passed on to the inherited mcglm method by
  `summary` and `coef`; passed to
  [`mc_plot_network`](https://evandeilton.github.io/mcglmgraph/reference/mc_plot_network.md)
  by `plot` with `type = "network"`.

- x:

  an object of class `"mcglm_graph"`.

- object:

  an object of class `"mcglm_graph"`.

- variability:

  for `summary`: the covariance of the graph parameters behind the
  standard errors, `"model"` (default), `"empirical"` for the sandwich
  or `"two_stage"`, as in
  [`mc_vcov_graph`](https://evandeilton.github.io/mcglmgraph/reference/mc_vcov_graph.md).

- cluster, small_sample:

  for `summary`: passed to
  [`mc_vcov_graph`](https://evandeilton.github.io/mcglmgraph/reference/mc_vcov_graph.md);
  they apply to `variability = "empirical"` and `"two_stage"` only and
  are refused with `"model"`.

- type:

  for `coef`: optional, one of `"regression"`, `"graph"` or
  `"dispersion"`, or a `type` accepted by the mcglm method. For `plot`:
  `"network"` for the dependence network or `"residuals"` for residual
  histograms.

## Value

An object of class `c("mcglm_graph", "mcglm")`. In addition to all the
standard `mcglm` components, it carries a list `$graph` with the
components

- `type`, `spec`: the layer (`"ggm"` or `"dag"`) and the canonical
  [`mc_between_graph`](https://evandeilton.github.io/mcglmgraph/reference/mc_between_graph.md)
  specification, with the number of responses recorded;

- `par`, `vcov`, `J`, `npar`, `edge_pos`: the named numeric vector of
  graph parameters (`eta_kk` then `eta_ij` for the GGM layer;
  `gamma_ji`, the coefficient of response \\i\\ in the equation of
  response \\j\\, then `delta_j` for the DAG layer), their model-based
  covariance \\J_n^{-1}\\, the accumulated sensitivity matrix \\J_n\\,
  the number of parameters, and the integer positions of the edge
  parameters in `par`. The matrix `vcov` carries the attributes
  `"inversion_status"` and `"inversion_ridge"`;

- `Sigma_b`, `K_b`: the fitted \\R \times R\\ between-response
  covariance and precision, in the order of the responses;

- `Dlist`: the list of `npar` matrices \\\dot\Sigma\_{b,d}\\, the
  derivatives of `Sigma_b` with respect to each element of `par`;

- `Shat`, `residuals`, `n`, `R`, `response_names`: the empirical
  covariance of the Pearson residuals \\\hat S_e\\, the \\n \times R\\
  matrix of residuals (not centered), the number of units and responses,
  and the response labels;

- `iter`, `converged`, `status`, `solver`: the number of iterations
  (zero for the closed-form DAG solution), the logical convergence
  indicator, a status string (`"converged"` or `"closed_form"` for a
  certified fit, otherwise the reason why it is not, such as
  `"max_iter"`, `"spd_line_search_failed"`, `"linear_solve_failed"`,
  `"score_not_converged"` or `"sensitivity_ill_conditioned"`), and the
  solver used (`"newton"` or `"closed"`);

- `score_norm`, `score_norm_unit`, `rcond_J`, `score_tol`, `rcond_tol`,
  `convergence_gates`, `stopping_rule`: the largest absolute coordinate
  of the estimating function at the solution (accumulated and divided by
  \\n\\), the reciprocal condition number of \\J_n\\, the two
  thresholds, the named logical vector of the three conditions
  (`parameter_change`, `score`, `sensitivity`) and a character
  description of the rule, so that a fit stopped by the parameter-change
  rule can be checked against the size of the estimating function it
  leaves;

- `inversion_status`, `inversion_ridge`, `linear_solve_status`,
  `linear_solve_ridge`, `inference_ok`: whether the inversion of \\J_n\\
  and the linear solves of the iterations were exact (`"exact"`, or
  `"exact_blockwise"` for the closed-form DAG solution) or
  `"ridge_regularized"`, the ridge added in each (zero when exact), and
  a logical that is `TRUE` when the inversion was exact and the three
  convergence conditions hold, in which case graph-layer inference is
  certified;

- `meta`, `init_info`, `control`: a list of builder metadata (layer,
  numbers of responses, edges and parameters, and the bases of the GGM
  layer or the ordering and edge matrices of the DAG layer), the
  diagnostics of the starting value (a list for the DAG layer solved by
  Newton iterations, `NULL` otherwise) and the control list used.

A graph layer that did not converge is returned with a warning.

## Details

**First stage.** The mean and within-response dispersion parameters are
estimated by [`mcglm`](https://rdrr.io/pkg/mcglm/man/mcglm.html) with
the quasi-score and Pearson estimating functions of Bonat and Jorgensen
(2016) [doi:10.1111/rssc.12145](https://doi.org/10.1111/rssc.12145) ,
under the usual correlation-based between-response structure of that
package. The \\n \times R\\ matrix \\E\\ of Pearson residuals, with rows
\\r_i^\top\\, is then formed, and \$\$\hat S_e =
\frac{1}{n}\sum\_{i=1}^{n}(r_i - \bar r)(r_i - \bar r)^\top\$\$ is its
empirical covariance, centered at the column means \\\bar r\\ and with
divisor \\n\\.

**Second stage.** The between-response covariance \\\Sigma_b\\ with
precision \\K_b = \Sigma_b^{-1}\\ is estimated from \\\hat S_e\\ by the
second-order Pearson estimating function \$\$\psi_d(\theta) =
-\frac{n}{2}\left\\ \mathrm{tr}(\Sigma_b^{-1} \dot\Sigma\_{b,d}) -
\mathrm{tr}(\Sigma_b^{-1} \dot\Sigma\_{b,d} \Sigma_b^{-1} \hat S_e)
\right\\,\qquad \dot\Sigma\_{b,d} =
\partial\Sigma_b/\partial\theta_d,\$\$ whose sensitivity matrix, in the
accumulated scale, has entries \$\$\[J_n\]\_{da} =
\frac{n}{2}\\\mathrm{tr}(\Sigma_b^{-1} \dot\Sigma\_{b,d} \Sigma_b^{-1}
\dot\Sigma\_{b,a}).\$\$ The vector of graph parameters is \\\theta =
\eta\\ for the GGM layer, with \\K_b = \sum_d \eta_d Z_d\\ (see
[`mc_precision_bases`](https://evandeilton.github.io/mcglmgraph/reference/mc_precision_bases.md)),
and \\\theta = (\gamma, \delta)\\ for the DAG layer, with \\\Sigma_b =
T^{-1} D T^{-\top}\\ (see
[`mc_between_graph`](https://evandeilton.github.io/mcglmgraph/reference/mc_between_graph.md)).
The order of the parameters within \\\theta\\ is given in
[`mc_build_sigma_between_graph`](https://evandeilton.github.io/mcglmgraph/reference/mc_build_sigma_between_graph.md).
Both layers keep \\R\\ free scale parameters, the diagonal of \\K_b\\ or
the innovation variances \\\delta\\: \\\Sigma_b\\ is the covariance of
the Pearson residuals and is not normalized to a correlation matrix.
Partial correlations and the support of the precision are unaffected by
that scale.

In the GGM layer, \\\dot\Sigma\_{b,d} = -\Sigma_b Z_d \Sigma_b\\ reduces
the two expressions to \$\$\psi\_{\eta_d} =
\frac{n}{2}\\\mathrm{tr}\\Z_d(\Sigma_b - \hat S_e)\\, \qquad
\[J_n\]\_{da} = \frac{n}{2}\\\mathrm{tr}(Z_d \Sigma_b Z_a \Sigma_b),\$\$
so that the root matches \\\Sigma_b\\ to \\\hat S_e\\ on the diagonal
and on the edges of the graph, which are the covariance selection
equations of Dempster (1972). The equations are solved by Newton
iterations \\\eta \leftarrow \eta + \alpha J_n^{-1}\psi\\, started at
the corresponding entries of \\\hat S_e^{-1}\\, with the step \\\alpha\\
halved from 1, at most 30 times, until the smallest eigenvalue of
\\K_b\\ exceeds \\10^{-10}\\ (see
[`mc_spd_line_search`](https://evandeilton.github.io/mcglmgraph/reference/mc_spd_line_search.md));
when no such step is found the iteration stops at the current iterate
with status `"spd_line_search_failed"`. The DAG layer has a closed-form
root, node by node (see
[`mc_dag_closed_root`](https://evandeilton.github.io/mcglmgraph/reference/mc_dag_closed_root.md)),
and can also be solved by Newton iterations on \\(\gamma, \log\delta)\\.

A fit is declared converged when three conditions hold: the relative
change of the parameters in the last iteration falls below `tol`, the
largest absolute coordinate of the estimating function, divided by
\\n\\, is at most `score_tol`, and the reciprocal condition number of
the sensitivity matrix exceeds `rcond_tol`. The closed-form DAG solution
takes no step, so there the first condition is replaced by the
requirements that no innovation variance is floored and no parent block
is ill-conditioned.

The graph layer is stored in `$graph` and does not replace the
covariance components of the first-stage `"mcglm"` object. Consequently,
inherited mcglm methods such as `vcov`, `gof`, `pAIC`, `pBIC` and
`plogLik` refer to the first-stage fit; graph-specific results are
exposed through the `mc_*` helpers. Because the graph layer is estimated
from fitted residuals, its default (model-based) covariance does not
account for the estimation of the first stage; see
[`mc_vcov_graph`](https://evandeilton.github.io/mcglmgraph/reference/mc_vcov_graph.md)
for the corrected forms.

## Methods (by generic)

- `print(mcglm_graph)`: Print the fit, highlighting the graph layer: the
  layer, the number of responses, units and edges, the convergence of
  the between-response engine and the between-response covariance
  \\\boldsymbol\Sigma_b\\. Returns `x` invisibly.

- `summary(mcglm_graph)`: Summary combining the mcglm summary with the
  graph layer (graph parameters, partial correlations and information
  criteria). Returns an object of class `"summary.mcglm_graph"`, a list
  with the first-stage summary in `first_stage` and its printed lines in
  `first_stage_text`; `layer`, `R`, `n` and `n_edges`; `converged`,
  `status` and `inference_ok`; the `variability` used; the data frame
  `coefficients` of the graph parameters, with columns `Estimate`,
  `Std.Error`, `z` (their ratio) and `p.value` (two-sided, standard
  normal reference, for the null value zero of each parameter taken
  alone), the last three `NA` when the fit is not certified for
  inference; `partial_cor`, the matrix of partial correlations; and
  `criteria`, the AIC, BIC and QIC of the graph layer. Its print method
  shows the first stage and then the graph layer, and returns its
  argument invisibly.

- `coef(mcglm_graph)`: Extract coefficients. With no `type`, or with a
  `type` of
  [`coef.mcglm`](https://rdrr.io/pkg/mcglm/man/coef.mcglm.html)
  (`"beta"`, `"tau"`, `"power"`, `"correlation"`), the inherited mcglm
  method answers, with its arguments `std.error` and `response` passed
  through `...`, and returns its data frame; otherwise `type`, which may
  be abbreviated, selects the regression, graph or dispersion
  coefficients, returned as plain numeric vectors.

- `plot(mcglm_graph)`: Diagnostic and network plots. With
  `type = "network"` it draws the dependence network with
  [`mc_plot_network`](https://evandeilton.github.io/mcglmgraph/reference/mc_plot_network.md)
  and returns invisibly the matrix of edge weights that function
  returns; with `type = "residuals"` it draws one histogram of the
  Pearson residuals per response, restores the graphical parameters on
  exit and returns `x` invisibly.

## Functions

- `print(summary.mcglm_graph)`: Print the summary: the first stage as
  mcglm prints it, then the graph layer.

## Two spellings of the same function

`mcglmgraph()` and `mcglm_graph()` are the same object: the second name
is bound to the first, not wrapped around it, so there is no forwarding
layer and the two cannot drift apart. `mcglmgraph()` matches the package
name and is the preferred spelling for new code; `mcglm_graph()` is kept
for existing scripts. Neither is deprecated.

The class of the returned object is `"mcglm_graph"` under both
spellings, and it is the class on which the S3 methods are registered.

## References

Bonat, W. H. and Jorgensen, B. (2016). Multivariate covariance
generalized linear models. *Journal of the Royal Statistical Society:
Series C (Applied Statistics)*, 65(5), 649-675.
[doi:10.1111/rssc.12145](https://doi.org/10.1111/rssc.12145)

Dempster, A. P. (1972). Covariance selection. *Biometrics*, 28(1),
157-175. [doi:10.2307/2528966](https://doi.org/10.2307/2528966)

Lauritzen, S. L. (1996). *Graphical Models*. Oxford University Press.

## See also

[`mc_between_graph`](https://evandeilton.github.io/mcglmgraph/reference/mc_between_graph.md)
to specify the layer,
[`mc_partial_cor`](https://evandeilton.github.io/mcglmgraph/reference/mc_partial_cor.md),
[`mc_partial_cor_delta`](https://evandeilton.github.io/mcglmgraph/reference/mc_partial_cor_delta.md)
and
[`mc_get_graph`](https://evandeilton.github.io/mcglmgraph/reference/mc_get_graph.md)
to interpret it,
[`mc_vcov_graph`](https://evandeilton.github.io/mcglmgraph/reference/mc_vcov_graph.md)
for the covariance of the graph parameters, and
[`mc_build_sigma_between_graph`](https://evandeilton.github.io/mcglmgraph/reference/mc_build_sigma_between_graph.md)
for the map from the parameters to \\\Sigma_b\\.

## Examples

``` r
set.seed(1); n <- 200; x <- rnorm(n)
y1 <- 1 + 0.5 * x + rnorm(n)
y2 <- 0.4 * y1 + 0.3 * x + rnorm(n)
y3 <- 0.2 * y2 + rnorm(n)
dat <- data.frame(y1, y2, y3, x)
fit <- mcglm_graph(list(y1 ~ x, y2 ~ x, y3 ~ x),
                   between_graph = mc_between_graph("ggm"), data = dat)
#> Automatic initial values selected.
fit
#> Multivariate Covariance GLM with graph-structured between-response covariance
#> Layer: GGM | responses: 3 | n: 200 | edges: 3
#> Between-response engine converged: TRUE (1 iterations)
#> 
#> Between-response covariance (Sigma_b):
#>        y1     y2     y3
#> y1 0.9900 0.3303 0.0632
#> y2 0.3303 0.9900 0.1647
#> y3 0.0632 0.1647 0.9900
#> 
#> This is also a valid 'mcglm' object: use coef(), vcov(), gof(), ...
coef(fit)            # inherited mcglm method
#>      Estimates Parameters        Type Response
#> 1   1.04145787     beta10        beta        1
#> 2   0.47692216     beta11        beta        1
#> 3   0.37196615     beta20        beta        2
#> 4   0.56927839     beta21        beta        2
#> 5  -0.02812858     beta30        beta        3
#> 6   0.19431290     beta31        beta        3
#> 7   0.33358711      rho12 correlation       NA
#> 8   0.06386570      rho13 correlation       NA
#> 9   0.16631396      rho23 correlation       NA
#> 10  1.02622172      tau11         tau        1
#> 11  1.28803308      tau21         tau        2
#> 12  1.22256601      tau31         tau        3
mc_partial_cor(fit)  # graph layer
#>             y1        y2          y3
#> y1 1.000000000 0.3281972 0.009020603
#> y2 0.328197212 1.0000000 0.154134931
#> y3 0.009020603 0.1541349 1.000000000
```
