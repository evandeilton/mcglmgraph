# Changelog

## mcglmgraph 0.8.5

- [`summary()`](https://rdrr.io/r/base/summary.html) of a fit and
  [`mc_dag_coefficients()`](https://evandeilton.github.io/mcglmgraph/reference/mc_dag_coefficients.md)
  take `method`, the route of the two-stage correction, and report the
  route taken in `two_stage_method`.
- [`mc_wald_graph()`](https://evandeilton.github.io/mcglmgraph/reference/mc_wald_graph.md)
  and
  [`mc_test_adequacy()`](https://evandeilton.github.io/mcglmgraph/reference/mc_test_adequacy.md)
  refuse a `method` other than `"auto"` when `type` is not
  `"two_stage"`, as
  [`mc_vcov_graph()`](https://evandeilton.github.io/mcglmgraph/reference/mc_vcov_graph.md)
  does.
- [`mc_vcov_graph()`](https://evandeilton.github.io/mcglmgraph/reference/mc_vcov_graph.md):
  `sensitivity = "expected"` is refused where it has no effect, with
  `type` other than `"two_stage"` or with `method = "shortcut"`, and
  gives a warning when `method = "auto"` takes the shortcut.
- [`mc_pc_skeleton()`](https://evandeilton.github.io/mcglmgraph/reference/mc_pc_skeleton.md)
  skips a conditioning set that gives no valid test, as
  [`mc_separating_sets()`](https://evandeilton.github.io/mcglmgraph/reference/mc_separating_sets.md)
  does, instead of stopping: the set removes no edge, and the skipped
  tests are counted in `diagnostics$n_unusable` and reported in one
  warning of class `"mcg_unusable_tests"`.
- [`mc_graph_pipeline()`](https://evandeilton.github.io/mcglmgraph/reference/mc_graph_pipeline.md)
  has the gate `complete_separation`, which fails when the skeleton
  search skipped a test, so such an output is not labeled a CPDAG.
  [`mc_cpdag_stability()`](https://evandeilton.github.io/mcglmgraph/reference/mc_cpdag_stability.md)
  counts those resamples in `B_incomplete` and reports them in one
  warning.

## mcglmgraph 0.8.4

- [`mc_test_independence()`](https://evandeilton.github.io/mcglmgraph/reference/mc_test_independence.md)
  and
  [`mc_partial_cor_delta()`](https://evandeilton.github.io/mcglmgraph/reference/mc_partial_cor_delta.md)
  take `method`, the route of the two-stage correction, as
  [`mc_vcov_graph()`](https://evandeilton.github.io/mcglmgraph/reference/mc_vcov_graph.md)
  does.
- `method = "auto"` of the two-stage correction takes the full form,
  when it is available, for a fit in which a response has a variance
  function other than `"constant"`. Every function that uses the
  two-stage variance follows it; `method = "shortcut"` gives the
  shortcut for such a fit.
- [`mc_partial_cor_delta()`](https://evandeilton.github.io/mcglmgraph/reference/mc_partial_cor_delta.md):
  the rows of the `"gradient"` attribute for pairs outside the support
  hold 0, the gradient of a constant, instead of `NA`.
- [`mc_wald_graph()`](https://evandeilton.github.io/mcglmgraph/reference/mc_wald_graph.md)
  refuses a contrast whose variance is at most 1e-6 of its model-based
  variance, which covers a parameter the first stage fixes under either
  route of the two-stage variance.

## mcglmgraph 0.8.3

- [`mc_spd_line_search()`](https://evandeilton.github.io/mcglmgraph/reference/mc_spd_line_search.md):
  the default `max_halvings` is 30, the value
  [`mcglm_graph()`](https://evandeilton.github.io/mcglmgraph/reference/mcglm_graph.md)
  uses internally.
- [`mc_partial_cor_delta()`](https://evandeilton.github.io/mcglmgraph/reference/mc_partial_cor_delta.md):
  a pair outside the specified support (an absent GGM edge, or a pair
  outside the moral graph of a DAG) returns `pcor = 0` and `NA` for
  `se`, `statistic`, `p_value`, `lower` and `upper`, and its row of the
  `"gradient"` attribute is `NA`. Rows of present pairs are unchanged.
- [`mc_partial_cor_delta()`](https://evandeilton.github.io/mcglmgraph/reference/mc_partial_cor_delta.md):
  a pair whose partial correlation is undefined (a nonpositive product
  of precision diagonals) returns `NA` in `pcor` and in the inference
  columns instead of 0. A variance that is negative beyond rounding,
  because the covariance of the parameters is not positive semidefinite,
  gives `NA` rather than 0.
- Renamed, with no alias: the attribute `centre` of
  `mc_vcov_graph(type = "two_stage")` is now `center`, and the column
  `prop_centre` of the v-structure tables of
  [`mc_find_vstructures()`](https://evandeilton.github.io/mcglmgraph/reference/mc_find_vstructures.md)
  and
  [`mc_graph_pipeline()`](https://evandeilton.github.io/mcglmgraph/reference/mc_graph_pipeline.md)
  is now `prop_center`.

## mcglmgraph 0.8.2

- [`coef()`](https://rdrr.io/r/stats/coef.html) on an `mcglm_graph`
  object accepts the `type` values of `coef.mcglm()` (`"beta"`, `"tau"`,
  `"power"`, `"correlation"`), so
  [`summary()`](https://rdrr.io/r/base/summary.html) prints the
  first-stage tables before the graph layer and
  [`anova()`](https://rdrr.io/r/stats/anova.html) works. A failure of
  the first-stage summary raises a warning.

- [`mc_graph_bootstrap()`](https://evandeilton.github.io/mcglmgraph/reference/mc_graph_bootstrap.md),
  `mc_graph_stability(refit = "full")` and
  `mc_cpdag_stability(refit = "full")` return, in every resample, the
  first stage refitted on that resample. A `statistic` that reads the
  regression or dispersion parameters gets their bootstrap spread.

- A refit of the graph layer on a resample or on another support
  refreshes every convergence and conditioning diagnostic of the
  `$graph` block.

- [`mc_select_lambda()`](https://evandeilton.github.io/mcglmgraph/reference/mc_select_lambda.md)
  and
  [`mc_graph_penalized()`](https://evandeilton.github.io/mcglmgraph/reference/mc_graph_penalized.md):
  supports from independent and warm-started grids can differ for pairs
  near the penalty threshold, since certified solutions agree only up to
  the tolerance.

- [`mc_graph_penalized()`](https://evandeilton.github.io/mcglmgraph/reference/mc_graph_penalized.md)
  accepts a step when the sufficient-decrease condition of Beck and
  Teboulle (2009) holds, so the penalized objective does not increase
  and the iteration reaches its parameter tolerance along the grid of
  [`mc_select_lambda()`](https://evandeilton.github.io/mcglmgraph/reference/mc_select_lambda.md)
  and the path of `mc_graph_stability(q = )`.

- [`mc_select_lambda()`](https://evandeilton.github.io/mcglmgraph/reference/mc_select_lambda.md)
  takes `max_iter`, `tol` and `kkt_tol` and passes them to the penalized
  fit of every grid point.

- The support of a penalized fit is the set of edge coefficients that
  are not exactly zero.

- [`mcglm_graph()`](https://evandeilton.github.io/mcglmgraph/reference/mcglm_graph.md)
  warns when the graph layer did not converge, and says so when the
  number of units does not exceed the number of responses.

- [`mc_wald_graph()`](https://evandeilton.github.io/mcglmgraph/reference/mc_wald_graph.md),
  and with it
  [`mc_test_independence()`](https://evandeilton.github.io/mcglmgraph/reference/mc_test_independence.md)
  and
  [`mc_test_adequacy()`](https://evandeilton.github.io/mcglmgraph/reference/mc_test_adequacy.md),
  refuses a contrast whose covariance is singular on its range, such as
  a parameter the first stage fixes under `type = "two_stage"`.

- A response transformed on the left side of its formula, as
  `log(y) ~ x`, is accepted and labeled by its expression.

- An identity within-response structure stored as a sparse matrix is
  recognized at any size.

- [`mc_graph_stability()`](https://evandeilton.github.io/mcglmgraph/reference/mc_graph_stability.md)
  checks `criterion` on entry.

- The examples generate their data with base R.

- [`summary()`](https://rdrr.io/r/base/summary.html) returns an object
  of class `"summary.mcglm_graph"` with a print method, takes
  `variability` for the standard errors of the graph parameters and
  shows the status of a graph layer that did not converge.

- [`mc_partial_cor_delta()`](https://evandeilton.github.io/mcglmgraph/reference/mc_partial_cor_delta.md)
  takes `interval = "fisher"` and
  [`mc_dag_coefficients()`](https://evandeilton.github.io/mcglmgraph/reference/mc_dag_coefficients.md)
  takes `r2_interval = "logit"`, for bounds that stay inside the range
  of the quantity.

- [`mc_edge_scores()`](https://evandeilton.github.io/mcglmgraph/reference/mc_edge_scores.md)
  computes only the block asked for in `which`, and the efficient score
  of an absent edge in closed form, which shortens
  [`mc_graph_stepwise()`](https://evandeilton.github.io/mcglmgraph/reference/mc_graph_stepwise.md).

- [`mc_plot_network()`](https://evandeilton.github.io/mcglmgraph/reference/mc_plot_network.md)
  stops each arrow of a DAG at the border of its node, so the direction
  is visible, and bounds the edge width.

- The skeleton and the enumeration of separating sets apply one kernel
  for the conditional-association test, which refuses a block of the
  covariance that is not positive definite.

## mcglmgraph 0.8.1

- `mc_vcov_graph(type = "two_stage", method = "full")` with the default
  `sensitivity = "observed"` differentiates the first-stage equations as
  `mcglm` solved them, with its bias correction when the fit used it.
  Standard errors of parameters the first stage nearly fixes, such as
  the innovation variance of a root node, match their sampling spread.
- The two-stage shortcut divides the propagation factor of each response
  by its average squared Pearson residual, the value the first stage
  fixes. A parameter the first stage fixes exactly gets a zero standard
  error, and where the shortcut is exact it matches the full form with
  the observed sensitivity.
- The two-stage shortcut is offered only when the prior weights are
  constant; otherwise `method = "auto"` takes the full form and
  `method = "shortcut"` is refused.
- The help page of
  [`mc_vcov_graph()`](https://evandeilton.github.io/mcglmgraph/reference/mc_vcov_graph.md)
  states when `sensitivity = "expected"` can be far off and recommends
  the observed sensitivity.

## mcglmgraph 0.8.0

- [`mc_between_graph()`](https://evandeilton.github.io/mcglmgraph/reference/mc_between_graph.md)
  stores the edges in canonical order (`i < j`, rows sorted), the order
  of the edge parameters. Edge scores, Wald statistics and the
  candidates ranked by `mc_graph_stepwise(direction = "forward")` or
  `"both"` do not depend on the order in which the edges are supplied,
  and a pair given as `(j, i)` is accepted.
- The DAG layer fits the empty support, `edges = matrix(0L, 0, 2)`, with
  `delta_j = S_jj` for every response; print, summary, the covariances,
  [`mc_dag_coefficients()`](https://evandeilton.github.io/mcglmgraph/reference/mc_dag_coefficients.md)
  and the criteria handle it.
- A near-singular parent block in the closed-form DAG root (reciprocal
  condition number of the block rescaled to a unit diagonal at or below
  `rcond_tol`) is reported in `status`, which names the affected
  responses, and the fit is not certified for inference. The GES local
  score refuses such a parent set. The test, the fallback ridge and the
  innovation-variance floor are relative to the scale of each response,
  so the root is equivariant to a rescaling of the responses.
- [`mcglm_graph()`](https://evandeilton.github.io/mcglmgraph/reference/mcglm_graph.md)
  refuses missing or non-finite responses, missing covariates and
  non-finite `offset`, `weights` or `Ntrial` values with an error naming
  the response and the number of rows, and refuses a fit whose Pearson
  residuals are not finite; complete cases are required.
- The Wald statistic of an edge whose estimated variance is not positive
  is `NA` instead of `0`, and
  [`mc_graph_stepwise()`](https://evandeilton.github.io/mcglmgraph/reference/mc_graph_stepwise.md)
  does not take that edge as the weakest one; it warns instead.
- `method = "auto"` of the two-stage correction, on a DAG with an
  immorality (where the shortcut is only first-order for the parameters
  of the collider and its descendants), takes the full form when it is
  available and, when it is not, the shortcut with a warning. Other fits
  keep their route.

## mcglmgraph 0.7.0

- [`mc_graph_penalized()`](https://evandeilton.github.io/mcglmgraph/reference/mc_graph_penalized.md),
  [`mc_select_lambda()`](https://evandeilton.github.io/mcglmgraph/reference/mc_select_lambda.md)
  and stability selection refuse, with a clear error, a response whose
  Pearson residuals have no variance: the diagonal start and the largest
  penalty of the path require `S_kk > 0` for every response. Previously
  the start held an infinite parameter and the failure surfaced later as
  a decomposition or convergence error.
- Documentation:
  [`mc_vcov_graph()`](https://evandeilton.github.io/mcglmgraph/reference/mc_vcov_graph.md),
  [`mc_dag_sensitivity()`](https://evandeilton.github.io/mcglmgraph/reference/mc_dag_sensitivity.md)
  and
  [`mc_dag_coefficients()`](https://evandeilton.github.io/mcglmgraph/reference/mc_dag_coefficients.md)
  name the curvature of the empirical forms uniformly as the observed
  Hessian of the working criterion. No computation changed.

## mcglmgraph 0.6.1

- For a sparse DAG, the empirical sandwich, full and shortcut two-stage
  covariances, and Godambe penalties now use the observed negative
  Hessian of the residual working criterion. The model-based Pearson
  Gram matrix remains available in fit$`graph`$J and is still used by
  Fisher scoring and the model-based covariance. With misspecified
  second moments, the observed parent block uses the residual sample
  covariance rather than the covariance implied by the DAG. Point
  estimates and the GGM computations are unchanged.
- QIC and GBIC can aggregate Pearson scores by independent cluster
  through cluster and small_sample in mc_compute_criterion(),
  mc_graph_stepwise() and mc_select_lambda(). BIC, GBIC and EBIC use the
  number of independent clusters in their logarithmic penalty when
  cluster is supplied. Existing calls retain their previous defaults.
- The documentation now states that the QIC optimism derivation assumes
  a genuinely known first stage. Same-sample estimation of marginal
  parameters needs an additional joint expansion. It also states the
  per-unit or per-cluster score-centering condition for empirical
  Godambe variability; aggregate centering alone does not make that
  variability consistent. These statistical conditions cannot be
  certified from a fitted object.
- A regression test compares the sparse-DAG Hessian to numerical
  differentiation under deliberate covariance misspecification and
  checks its use in empirical sandwich, two-stage covariance and
  criteria.

## mcglmgraph 0.6.0

The correction for the first stage reaches every fit, the tests of the
graph layer take any contrast, and the documentation says under which
conditions each result holds. The default routes return what 0.5.1
returned.

### The full two-stage correction

- `mc_vcov_graph(type = "two_stage")` gains
  `method = c("auto", "shortcut", "full")`. The full form rebuilds the
  first stage at its estimate with the routines of `mcglm`, splits its
  quasi-score and Pearson functions over units, and forms the two-stage
  influence function `H^-1 (psi_2 - B A^-1 psi_1)`, with `B` by the
  chain rule on the Pearson residual and no condition on the mean or on
  the residual covariance. It covers what the shortcut refuses: a
  Poisson-Tweedie power that is estimated, more than one covariance
  component, and within-response structures that are not the identity.
  `"auto"` takes the shortcut where it applies, so a number computed
  before is computed again.
- The sensitivity of the first stage is the observed one by default, the
  Jacobian of its estimating functions at the estimate by central
  differences of their totals, rebuilt with `mcglm`;
  `sensitivity = "expected"` takes the model-based form instead. The two
  agree when the first stage is correctly specified, and only the
  observed one keeps the correction right when it is not: on a gamma
  design whose third response misses a covariate, the corrected standard
  errors stay within 5% of the sampling spread with it and drift by up
  to 20% without it. The observed form inherits the conditioning of the
  first stage: on a weakly identified Poisson-Tweedie fit it can inflate
  the correction several times, which is why the attribute
  `sensitivity_rcond` reports it.
- Either way the block of the covariance parameters on the regression
  parameters is kept. The `joint_inv_sensitivity` that `mcglm` reports
  sets it to zero, because the covariance it builds at the end of the
  fit carries no derivative in the means. That block is not zero when
  the variance depends on the mean, and without it the regression part
  of the first stage stops cancelling: on a gamma design with n = 20000
  the diagonal coefficients come out 56% to 99% too wide. With it, the
  full form and the shortcut agree to 1e-4 wherever both apply.
- With a structure that is not the identity the units are the clusters,
  which have to nest it: every matrix of the matrix linear predictor
  must vanish between distinct clusters, and that is checked.
- [`mc_two_stage_available()`](https://evandeilton.github.io/mcglmgraph/reference/mc_two_stage_available.md)
  reports both routes, `shortcut_ok` and `full_ok` with their reasons,
  and the one `"auto"` takes. `ok` and `why` keep their meaning, the
  shortcut’s.
- Checked against the definitions: the unit contributions, the
  sensitivity and `B` by `numDeriv`, the power columns included. Checked
  against the sampling spread: on 400 replicates of a clustered count
  design, with random-effect blocks nested in 40 clusters, the full-form
  standard errors are 0.92 to 0.99 of the Monte Carlo deviation, where
  the naive sandwich is up to 2.2 times off; with an estimated
  Poisson-Tweedie power, 0.90 to 0.99 against 3.7 to 5.4.

### Resampling by clusters

- [`mc_graph_bootstrap()`](https://evandeilton.github.io/mcglmgraph/reference/mc_graph_bootstrap.md),
  [`mc_cpdag_stability()`](https://evandeilton.github.io/mcglmgraph/reference/mc_cpdag_stability.md)
  and `mc_graph_stability(refit = "full")` refit both stages for
  structures that are not the identity when `cluster` is given. Each
  drawn cluster keeps its own block of the structure, and a cluster
  drawn twice gives two blocks that do not touch, where the submatrix of
  the rows drawn would have joined the copies. Without `cluster` such a
  fit is refused, and the message no longer suggests `cluster` where it
  could not help.
- [`mc_graph_stability()`](https://evandeilton.github.io/mcglmgraph/reference/mc_graph_stability.md)
  takes `cluster` and then draws half of the clusters.

### Tests of the graph layer

- [`mc_wald_graph()`](https://evandeilton.github.io/mcglmgraph/reference/mc_wald_graph.md)
  tests any linear hypothesis `L theta = c` on the graph parameters,
  with the three covariances, degrees of freedom equal to the rank of
  `L` and a refusal of a right-hand side that redundant rows contradict.
- [`mc_test_adequacy()`](https://evandeilton.github.io/mcglmgraph/reference/mc_test_adequacy.md)
  refits the saturated graph on the same residuals and tests the omitted
  edges jointly, and marks the result as descriptive when the support
  was selected on the same data.
- [`mc_test_independence()`](https://evandeilton.github.io/mcglmgraph/reference/mc_test_independence.md)
  is
  [`mc_wald_graph()`](https://evandeilton.github.io/mcglmgraph/reference/mc_wald_graph.md)
  with the contrast that selects every edge, gains `type`, and returns
  the same statistic to the last bit.

### Residual degrees of freedom

- [`mc_residual_df()`](https://evandeilton.github.io/mcglmgraph/reference/mc_residual_df.md)
  returns `nu_e = n - h` in the exact case, every response Gaussian with
  identity link, constant variance and a common design of rank `h`, and
  says why otherwise. The separation routines take
  `df_residual = "auto"` where a fit is at hand and refuse a fit for
  which it is not exact. The default remains the approximation by `n`;
  with the intercept alone the two coincide.

### Directed layer

- [`mc_dag_coefficients()`](https://evandeilton.github.io/mcglmgraph/reference/mc_dag_coefficients.md)
  gives intervals for the standardized coefficient and for `R^2` by the
  delta method, with analytical gradients through `Sigma_b`. That is the
  new default, `interval = "delta"`; `"fixed_factor"` keeps the interval
  of 0.5.1, which carried the variances in the factor as if known.

### The L1 path and stability selection

- [`mc_select_lambda()`](https://evandeilton.github.io/mcglmgraph/reference/mc_select_lambda.md)
  takes `warm_start` and `include_zero`, and
  [`mc_graph_penalized()`](https://evandeilton.github.io/mcglmgraph/reference/mc_graph_penalized.md)
  takes `start`.
- [`mc_compute_criterion()`](https://evandeilton.github.io/mcglmgraph/reference/mc_compute_criterion.md)
  and the selectors accept `"ebic"`, the extended BIC of Foygel and
  Drton, with `ebic_gamma`.
- [`mc_graph_stability()`](https://evandeilton.github.io/mcglmgraph/reference/mc_graph_stability.md)
  takes `q`: every subsample selects the first `q` edges along the L1
  path, and the bound of Meinshausen and Buehlmann on the expected
  number of false edges is reported. `sampling = "complementary"` draws
  the pairs of Shah and Samworth.
- The L1 iteration starts from the diagonal fit of the empty graph and
  does not stop at a small step while the KKT residual fails. The
  certified solutions do not change.

### Documentation

- The conditions under which the first stage acts as a diagonal
  rescaling, a correctly specified mean and a residual covariance that
  does not depend on the covariates, are stated where the package relies
  on them, and so is what the checks of a fit cannot verify.
- [`mc_graph_stepwise()`](https://evandeilton.github.io/mcglmgraph/reference/mc_graph_stepwise.md)
  refuses a directed fit instead of converting it.
- The small-sample flag of
  [`mc_vcov_graph()`](https://evandeilton.github.io/mcglmgraph/reference/mc_vcov_graph.md)
  says whether `G/(G-1)` was applied.

## mcglmgraph 0.5.1

A one-line correction with a disproportionate effect, and the test that
pins it.

- `.is_identity_structure()` read the dimension of a within-response
  structure through [`dim()`](https://rdrr.io/r/base/dim.html). `Matrix`
  sits in Imports, so in a session that had not yet attached it the S4
  method was not dispatched and
  [`dim()`](https://rdrr.io/r/base/dim.html) answered `NULL` on the
  `ddiMatrix` the first stage builds, which sent the function down its
  `FALSE` branch before it could reach the `ddiMatrix` case written for
  exactly that object. The dimension and the unit-diagonal flag are now
  read from the slots with [`attr()`](https://rdrr.io/r/base/attr.html),
  which needs no namespace. Nothing else changes, and no dependency is
  added.
- What it fixed for the user. `mc_graph_stability(refit = "full")`
  stopped on the first call of a session, with a message about a matrix
  linear predictor that resampling would not preserve, and succeeded on
  the second call with the same input.
  [`mc_two_stage_available()`](https://evandeilton.github.io/mcglmgraph/reference/mc_two_stage_available.md)
  reported a spurious extra reason on the first call. Both came from the
  same line.
- A regression test runs the check in a fresh process through `callr`,
  since any earlier call in the same session loads the namespace and
  hides the fault.

## mcglmgraph 0.5.0

Three additions. Nothing that existed changes: the default fitting
route, the skeleton search and the standard orientation rule return
exactly what 0.4.0 returned, and the test suite pins that.

### Inference on the graph layer

- [`mc_vcov_graph()`](https://evandeilton.github.io/mcglmgraph/reference/mc_vcov_graph.md)
  exports the covariance of the graph parameters in its two forms, the
  model-based `J_n^-1` and the sandwich `J_n^-1 V_n J_n^-1`, with
  clustering and the small-sample factor. The sandwich was only
  reachable through an internal function, so a reported standard error
  could not be reproduced from the public interface.
  [`mc_test_independence()`](https://evandeilton.github.io/mcglmgraph/reference/mc_test_independence.md)
  now consumes this very matrix, so the accessor and the test cannot
  drift apart.
- The sandwich uses the conditioning tolerance the fit was solved under,
  instead of the default. It could otherwise accept a sensitivity the
  fitting engine had refused, or refuse one it had accepted.

### Order-independent orientation

- [`mc_separating_sets()`](https://evandeilton.github.io/mcglmgraph/reference/mc_separating_sets.md)
  enumerates every set that separates a pair of responses on a given
  skeleton, with the partial correlation and the p-value of each. The
  skeleton search records only the first such set, which is all the
  standard rule needs; two sets can separate the same pair with very
  different strength, and which one is kept is decided by the order of
  enumeration.
- [`mc_graph_pipeline()`](https://evandeilton.github.io/mcglmgraph/reference/mc_graph_pipeline.md)
  and
  [`mc_find_vstructures()`](https://evandeilton.github.io/mcglmgraph/reference/mc_find_vstructures.md)
  take `orient_rule`, with the conservative rule of Ramsey, Zhang and
  Spirtes (2006) and the majority rule of Colombo and Maathuis (2014)
  beside the standard one. Combined with the order-independent skeleton
  the package already used, they make the whole procedure independent of
  the labeling of the responses. Both are native and use the same
  conditional-association test and the same degrees-of-freedom contract
  as the skeleton, so a difference between two runs is attributable to
  the rule and not to a change of test.
- [`mc_orient_meek()`](https://evandeilton.github.io/mcglmgraph/reference/mc_orient_meek.md)
  takes `ambiguous`. A triple that a rule leaves undecided withholds R1
  and R3 where they would settle it by propagation. R2 forbids a cycle
  and involves no unshielded triple, so it is unaffected.
- On a chain whose extremes have marginal correlation small enough for
  the order-zero test to accept, the standard rule declares a collider
  that does not exist, and the two new rules leave the triple
  unoriented.

### The unit-diagonal constraint

- [`mc_dag_unit_diagonal()`](https://evandeilton.github.io/mcglmgraph/reference/mc_dag_unit_diagonal.md)
  returns the innovation variances under which the covariance of the
  directed layer is a correlation matrix, and says whether the
  coefficients admit them. The undirected layer fixes its diagonal at
  one at the cost of nothing; here the constraint confines the
  coefficients to the open set where every innovation variance stays
  positive, defined by as many nonlinear inequalities as there are
  responses. The set contains the origin and is a proper subset whenever
  there is an edge, which is why the fitting routines do not impose it.
  It also carries the rescaling to the one representative of the orbit
  with unit innovations.

### Two gates in C++

- The check for a directed cycle and the Dor and Tarsi check for a
  consistent extension are combinatorics on the adjacency, with no
  floating point in them, and both are evaluated once per resample by
  the structure-frequency routines. They moved to C++, and agree with
  the plain-R versions on 300 random partially directed graphs.
- Two pieces stay in R on purpose. The per-unit Pearson contributions
  were already reassociated into one BLAS call per parameter, and sit on
  the path that decides every selection: moving them would buy little
  and risk changing a criterion in its last bit, which is enough to flip
  an edge at a tie. The conditional-association test and the loop of the
  skeleton search stay for the same reason.

### Greedy equivalence search, native

- [`mc_ges()`](https://evandeilton.github.io/mcglmgraph/reference/mc_ges.md)
  searches over Markov equivalence classes by the two-phase procedure of
  Chickering: a forward phase from the empty class applying the
  insertion that raises the score most, then a backward phase applying
  the best deletion. No external library is involved. It answers the
  same question as
  [`mc_graph_pipeline()`](https://evandeilton.github.io/mcglmgraph/reference/mc_graph_pipeline.md)
  by other means: one trades every edge against a penalty and never
  tests a hypothesis, the other tests separations at a fixed level and
  then orients, and the two need not agree.
- The search is well posed only because the working criterion is score
  equivalent, two members of a class fitting the same covariance with
  the same parameter count. The test suite checks that on enumerated
  classes before checking anything else.
- The operators, the local scores and the extension of a partially
  directed graph to a DAG are in C++. Each operator’s result is mapped
  back to a CPDAG, so the search never leaves the space it is defined
  on.
- On thirty random DAGs with n = 3000, the search recovers the true
  class in 23 cases against 22 for the constraint-based route.

### Structural intervention distance

- [`mc_sid()`](https://evandeilton.github.io/mcglmgraph/reference/mc_sid.md)
  counts the ordered pairs of responses for which the estimated graph
  would send one to the wrong adjustment set. The Hamming distance
  counts the edges two graphs disagree about and says nothing about what
  the disagreement costs; two estimates at the same Hamming distance can
  still give the right adjustment set for every pair or for none. The
  `SID` package is not a dependency: the distance is native, with
  d-separation decided by ancestral moralization in C++.
- A partially directed estimate gives an interval rather than a number,
  over the members of its class, which is what a CPDAG admits. The
  members come from
  [`mc_class_members()`](https://evandeilton.github.io/mcglmgraph/reference/mc_class_members.md).
- Validity is decided by the generalized adjustment criterion and not by
  the back-door one. The back-door criterion is sufficient and not
  necessary: it refuses any set containing a descendant of the
  treatment, while a descendant lying on no causal path to the response
  does no harm. The difference is not cosmetic. Under the back-door
  reading the distance counted pairs the adjustment in fact recovers;
  under the generalized one, on 1992 pairs of random graphs, the count
  agrees exactly, in both directions, with whether the population
  regression returns the total causal effect.

### What a class identifies, and what it does not

- [`mc_class_members()`](https://evandeilton.github.io/mcglmgraph/reference/mc_class_members.md)
  enumerates the DAGs a CPDAG stands for, by orienting one undirected
  edge at a time, propagating with Meek’s rules, and abandoning a branch
  as soon as it creates a cycle or an immorality the class does not
  have. No member is produced twice, since the two branches fix opposite
  directions on the same edge.
- [`mc_class_coefficients()`](https://evandeilton.github.io/mcglmgraph/reference/mc_class_coefficients.md)
  fits the directed layer on the residual covariance under every member
  and reports, per arc, the range the class allows and whether the head
  is free of undirected edges, which is the condition under which a
  single value can be reported at all. In a five-response example with
  two members, the compelled arc takes 0.700 in one member and 0.431 in
  the other, and its sibling 0.500 against 0.780, with the closed form
  `(g - c h) / (1 + c^2)` matching.
- Every member is checked to fit the same covariance, which is the
  property that makes the class, and not the DAG, what second moments
  identify.

### Out-of-sample comparison

- [`mc_graph_cv()`](https://evandeilton.github.io/mcglmgraph/reference/mc_graph_cv.md)
  splits the units into folds, fits the graph layer outside each fold
  and scores it inside, so candidate supports are compared without a
  penalty term and without scoring a support on the sample that chose
  it. What is scored is the Gaussian working criterion per held-out
  unit, on the scale of a log density, so differences are comparable
  across fold sizes.

### Structure metrics inside the package

- [`mc_structure_metrics()`](https://evandeilton.github.io/mcglmgraph/reference/mc_structure_metrics.md)
  compares an estimated graph with a known one: true and false positive
  rates, precision, F1 and structural Hamming distance.
  [`mc_cpdag_metrics()`](https://evandeilton.github.io/mcglmgraph/reference/mc_cpdag_metrics.md)
  does the same for CPDAGs (CPDAG distance, skeleton F1 and SHD,
  v-structure recall and precision), and
  [`mc_true_cpdag()`](https://evandeilton.github.io/mcglmgraph/reference/mc_true_cpdag.md)
  returns the CPDAG of a known precision matrix.
- The measures keep apart what the skeleton recovers and what the
  orientation does, and a rate whose denominator is empty is `NA` rather
  than a number: precision is undefined when nothing was selected,
  recall when there was nothing to find, over-orientation when the truth
  directs every edge.
- The pair states, the v-structures and the confusion counts are in C++.

### Frequencies over resamples

- [`mc_graph_stability()`](https://evandeilton.github.io/mcglmgraph/reference/mc_graph_stability.md)
  gains `refit = "full"`, which draws each half-sample from the rows of
  the data and refits the first stage on it, so the selection
  frequencies carry the variability of both stages. The default stays
  `"residual"`, which resamples the stored residuals with the first
  stage held fixed and is what earlier versions did, so a recorded
  result keeps reproducing.
- [`mc_cpdag_stability()`](https://evandeilton.github.io/mcglmgraph/reference/mc_cpdag_stability.md)
  records, for every pair of responses, how often the pair came out
  absent, undirected, or directed each way. An orientation the data
  barely support shows up as a pair whose two directions split the
  frequency, which a single run cannot reveal. It carries background
  knowledge and the orientation rule through to every resample, and
  reports the share of outputs that passed the validity gates and the
  average number of ambiguous triples.
- Units are drawn with replacement there, keeping the sample size, which
  is the bootstrap convention; the selection frequencies keep
  half-samples, the convention of the stability-selection literature.
  The difference is stated in both help pages rather than left to be
  inferred.

### Resampling both stages

- [`mc_graph_bootstrap()`](https://evandeilton.github.io/mcglmgraph/reference/mc_graph_bootstrap.md)
  resamples units, or whole clusters, and refits the first-stage MCGLM
  and the graph layer in every resample. Nothing in the package
  resampled the first stage before:
  [`mc_graph_stability()`](https://evandeilton.github.io/mcglmgraph/reference/mc_graph_stability.md)
  reuses the stored residuals, which is the comparison and not the
  propagation. The support of the graph is held fixed, so this is an
  inferential check on a given structure and not a post-selection
  estimator.
- It is the route for the fits the analytical shortcut refuses, the
  Poisson-Tweedie ones with an estimated power among them, and it is the
  yardstick the shortcut is checked against. On a Gaussian design the
  bootstrap standard errors track `type = "two_stage"` and depart from
  the naive sandwich exactly on the diagonal coefficients, in the
  direction and by the amount the closed form predicts. With
  `stage = "residual"` it resamples the residuals alone and reproduces
  the naive sandwich, which is what isolates the contribution of the
  first stage.
- The first-stage call is rebuilt from the fitted object, so a fit whose
  within-response structure is not the identity is refused rather than
  resampled in a way that would not preserve it.
- The resamples are drawn in the parent before any dispatch, so a given
  seed yields the same resamples at any worker count.

### A parameter vector that is a vector

- The parameter block of a GGM fit was a one-column matrix, because that
  is what the line search returns, while that of a DAG fit was a plain
  vector. `coef(fit, type = "graph")` inherited the difference. Both are
  vectors now.

### Propagating the first stage

- `mc_vcov_graph(type = "two_stage")` corrects the variance of the graph
  parameters for the estimation of the residuals. Stacking the two sets
  of estimating functions gives a block triangular sensitivity whose
  inversion gives the influence function of the second stage; the cross
  sensitivity factorizes as `B = Q Ubar` with
  `Q[d, r] = -(1/2) [Z_d Sigma_b](r, r)`, and when every response
  carries a single dispersion parameter that multiplies its variance
  function the correction needs second-stage quantities only. Until now
  the package reported the two conditional forms and nothing else.
- [`mc_two_stage_available()`](https://evandeilton.github.io/mcglmgraph/reference/mc_two_stage_available.md)
  says whether a fit qualifies and, when it does not, why. A fit that
  does not qualify raises an error rather than returning the uncorrected
  matrix under the corrected name. Poisson-Tweedie is additive in its
  dispersion parameter and qualifies only with the power fixed at one.
- What the correction reaches is checked rather than asserted. The test
  suite rebuilds `Q` from the canonical directions of the precision,
  verifies that `H^-1 Q[, r]` is the derivative of the parameters along
  the diagonal rescaling of the precision, by finite differences, and
  that the gradients of the partial correlations annihilate it, which is
  why those are unaffected. With Gaussian residuals the ratio of
  corrected to naive standard error of a diagonal coefficient matches
  its closed form `sqrt(2 (K_kk - 1) / K_kk)`, and the correction to an
  absent edge vanishes.
- The per-unit Pearson contributions were split out of the empirical
  variability, so the correction and the sandwich share one piece of
  arithmetic rather than two copies of it.

### Background knowledge on the orientation

- [`mc_orient_meek()`](https://evandeilton.github.io/mcglmgraph/reference/mc_orient_meek.md)
  gains Meek’s fourth rule and the arguments `required` and `forbidden`.
  R1 to R3 are complete when the only oriented edges are those of the
  v-structures, which is why the CPDAG of a class never needed R4; it
  becomes necessary the moment edges are oriented from outside, because
  the graph is then no longer the CPDAG of any class. Adding it changes
  nothing without background knowledge: the CPDAG is a fixed point of
  the four rules, and the test suite checks that on random classes.
- [`mc_graph_pipeline()`](https://evandeilton.github.io/mcglmgraph/reference/mc_graph_pipeline.md)
  takes the same two arguments, by response name or by index, and
  reports what was imposed. With knowledge the output is no longer the
  CPDAG of an equivalence class but the maximally oriented graph, the
  union of the members of the class that satisfy the knowledge, and the
  `label` says so. Knowledge that contradicts a compelled edge, itself,
  or the skeleton is refused rather than silently dropped.
- The rules are checked against their definition and not against a
  second implementation: for random DAGs the test suite enumerates every
  orientation of the skeleton, keeps those Markov equivalent to the
  original that satisfy the knowledge, and requires the union to be
  exactly what the rules return.
- The four rules moved to C++. R4 adds a search over pairs, which would
  have made a sweep cost the fourth power of the number of responses in
  R.

### Closed-form directed layer

- [`mc_dag_closed_root()`](https://evandeilton.github.io/mcglmgraph/reference/mc_dag_closed_root.md)
  solves the Pearson equations of the directed layer node by node. With
  `directed = TRUE` the arcs are read as parent and child, any acyclic
  set is accepted and a compatible ordering is derived, which is what
  fitting an estimated DAG requires.
- [`mc_dag_sensitivity()`](https://evandeilton.github.io/mcglmgraph/reference/mc_dag_sensitivity.md)
  assembles the sensitivity from its block form, one block per node,
  instead of forming `P^2` traces of `R x R` products.
- [`mcglm_graph()`](https://evandeilton.github.io/mcglmgraph/reference/mcglm_graph.md)
  takes that route by default for the directed layer, including a
  blockwise inverse of the sensitivity. It solves the same equations
  exactly: on the cases in the test suite the two routes agree to 1e-15
  in the parameters and in their variance, and the closed route is about
  forty times faster at ten responses or more.
  `control = list(dag_solver = "newton")` restores the scoring
  iteration. A fit taken by the closed route reports
  `status = "closed_form"` instead of `"converged"`, since it took no
  step whose size a convergence rule could read; what certifies it is
  the score gate, which it passes by construction.
- [`mc_dag_coefficients()`](https://evandeilton.github.io/mcglmgraph/reference/mc_dag_coefficients.md)
  reports the scale-free quantities, the standardized coefficient and
  the proportion of variance accounted for by the parents, which are the
  ones invariant to the diagonal rescaling the first stage induces. The
  raw coefficients were all the package exposed.

### Documentation

- The edge matrix is sorted lexicographically before use and `gamma` is
  read in that sorted order, not in the order the rows were supplied.
  That was the behavior and it was not written down.

## mcglmgraph 0.4.0

Performance release. Nothing in the statistical contract changes: every
result below is the same arithmetic, reassociated, and the parallel
routines return exactly what the sequential ones return for a given
seed.

### Faster

- The empirical variability of the estimating functions used to evaluate
  its quadratic form one observation at a time, inside a double loop
  over `n` and the parameter count: one BLAS call per row of data.
  Written as one matrix product per parameter it costs `O(n R^2)` in a
  single call. Measured 31x to 45x on the sizes this package is used at,
  agreeing to about 1e-14. It sits on the critical path of every
  Godambe-adjusted criterion, so the stepwise search under GBIC and the
  L1 path both inherit the gain: stability selection with ten responses,
  `n = 500` and `B = 60` went from 234.5 s to 63.9 s sequentially.

### Parallel, opt-in

- [`mc_graph_stability()`](https://evandeilton.github.io/mcglmgraph/reference/mc_graph_stability.md)
  and
  [`mc_select_lambda()`](https://evandeilton.github.io/mcglmgraph/reference/mc_select_lambda.md)
  accept `workers`. The same case above runs in 10.1 s with sixteen
  workers, twenty-three times the original. Both default to
  `workers = 1`.
- The default is deliberate. A simulation study usually parallelizes at
  the level of the simulation cell, and a second level of workers
  underneath oversubscribes the machine. Inside a PSOCK worker the
  argument is ignored with a warning, so the nesting cannot happen by
  accident.
- Reproducibility is preserved by construction: every subsample is drawn
  in the parent before any dispatch, so a given seed yields the same
  subsamples at any worker count and `workers = 16` returns exactly what
  `workers = 1` returns. Forking is used where available; Windows falls
  back to sequential.

## mcglmgraph 0.3.4

The strict input contract introduced in 0.3.3 covered only the low-level
LDL wrappers; the path that
[`mcglm_graph()`](https://evandeilton.github.io/mcglmgraph/reference/mcglm_graph.md)
actually takes reached the kernels through separate, permissive helpers.
A non-integer edge index was silently truncated there, so a caller could
ask for one edge and be fitted another with no signal. There is now a
single validator.

### Corrections

- `.validate_edges()` is the one edge contract of the package.
  [`mc_precision_bases()`](https://evandeilton.github.io/mcglmgraph/reference/mc_precision_bases.md),
  [`mc_between_graph()`](https://evandeilton.github.io/mcglmgraph/reference/mc_between_graph.md),
  the LDL wrappers and the fitting path all funnel through it, so the
  same malformed input is refused the same way whichever door it came
  through. Non-integer, non-finite, out-of-range, self-looping and
  duplicated edges are errors, not silent coercions.
- The GGM layer no longer deduplicates edges silently: a repeated pair
  does not define a second parameter, so it is rejected rather than
  dropped. This matches the behavior the DAG layer already had.
- [`mc_precision_bases()`](https://evandeilton.github.io/mcglmgraph/reference/mc_precision_bases.md)
  validates `R` as a single positive integer instead of coercing it,
  which turned `2.7` into `2` and a length-two vector into an internal
  condition-length error.
- [`mc_between_graph()`](https://evandeilton.github.io/mcglmgraph/reference/mc_between_graph.md)
  validates at construction time and rejects `ordering` when
  `type = "ggm"`: the precision support is unordered, and the argument
  was previously accepted and ignored without warning.

### Documentation

- [`mc_graph_penalized()`](https://evandeilton.github.io/mcglmgraph/reference/mc_graph_penalized.md)
  states that the regularization path is always built over the complete
  graph, whatever support the supplied fit carries.

## mcglmgraph 0.3.3

This release strengthens low-level validation, graph semantics,
convergence diagnostics and reproducibility contracts.

### Corrections

- DAG partial-correlation support follows concentration/moral support
  instead of equating absence of a directed arc with a structural zero.
- The LDL API rejects duplicated edges, inconsistent dimensions,
  non-finite values and non-positive innovation variances before
  entering C++ kernels.
- CPDAG helpers distinguish internal consistency from certification of
  an essential graph, use an explicit residual-degrees-of-freedom
  contract and expose regularization or insufficient-degrees-of-freedom
  states.
- Optional score-based orientation returns fields that describe the same
  graph, while retaining the PC skeleton separately for inspection.
- Inference and post-selection objects expose, recompute or invalidate
  numerical diagnostics rather than silently carrying stale values.
- The unused preconditioned proximal-Newton C++ kernel was removed; the
  public penalized route remains the Euclidean ISTA implementation
  certified by the parameter-change and KKT gates.
- Gaussian-copula documentation now requires an explicit latent-scale
  step; inversion of observed Pearson covariance is not presented as
  sufficient.

### New features

- [`mcglmgraph()`](https://evandeilton.github.io/mcglmgraph/reference/mcglm_graph.md)
  is now a second name for
  [`mcglm_graph()`](https://evandeilton.github.io/mcglmgraph/reference/mcglm_graph.md),
  spelled the way the package itself is. It is not a wrapper: the new
  name is bound to the existing object, so there is no forwarding layer
  and no way for the two spellings to drift apart. Neither name is
  deprecated, and the documentation is shared. The class of the returned
  object remains `"mcglm_graph"` under both spellings, because that is
  the identifier recorded inside saved fits and the one the S3 methods
  are registered on.

- [`mc_partial_cor_delta()`](https://evandeilton.github.io/mcglmgraph/reference/mc_partial_cor_delta.md)
  adds standard errors, Wald statistics, p-values and confidence
  intervals to the between-response partial correlations, by the delta
  method. The gradient has a closed form and the variance is the
  quadratic form `g' V g`.
  [`mc_partial_cor()`](https://evandeilton.github.io/mcglmgraph/reference/mc_partial_cor.md)
  returns the point estimate only.

  - `variability = "model"` uses `J^-1`; `"empirical"` uses the sandwich
    form `J^-1 (n Vhat) J^-1`.
  - It applies to both layers. The precision derivatives are obtained
    from `Kdot_d = -K Sigmadot_d K`, an identity that in the GGM layer
    returns exactly the canonical base `Z_d` of the statement, and in
    the DAG layer supplies the directions in `(gamma, delta)` without
    changing the formula.
  - Pairs outside the support are structural zeros: the gradient
    vanishes and the standard error is exactly zero, flagged by
    `present = FALSE`.
  - The documentation states the scope: the variance is conditional on
    the estimated residuals and does not propagate first-stage
    uncertainty; applied to a support chosen by search, it remains
    descriptive.

- New RcppArmadillo kernel
  [`cpp_pcor_delta()`](https://evandeilton.github.io/mcglmgraph/reference/cpp_pcor_delta.md).

### Tests

The analytical, numerical and graph-contract suites include adversarial
cases.

- Numerical verification of the analytical derivatives is no longer spot
  checking. `test-analytic-derivatives.R` sweeps **every** coordinate of
  both layers, in three configurations of dimension and support each:
  complete GGM with three responses, chain with four and sparse with
  five; complete DAG with three, chain with four and hub with five.
  Previously a single coordinate per layer was checked, which would let
  through a sign error confined to a diagonal base or to one specific
  edge orientation.
- The estimating function is now checked against the numerical gradient
  of the Gaussian working criterion, in both layers, with `Shat`
  different from `Sigma` so that the gradient is not trivially zero.
- The sensitivity is now checked against minus the numerical Jacobian of
  the estimating function, evaluated at the point where `Shat` coincides
  with the fitted covariance, which is where the identity holds exactly.
- Structural tests of the derivatives: the GGM diagonal derivatives are
  rank one and negative semidefinite; the `delta` derivatives of the DAG
  layer are rank one and positive semidefinite.
- The saturated equivalence between the GGM and DAG layers gained a unit
  test.
- `tests/testthat/setup.R` collects the configurations and the
  numerical-differentiation helpers.

## mcglmgraph 0.3.2

This version fixes public contracts that promised more formal validity
than the procedures delivered.

- [`mc_graph_pipeline()`](https://evandeilton.github.io/mcglmgraph/reference/mc_graph_pipeline.md)
  no longer labels its output as a CPDAG unconditionally. Orientation by
  v-structures now collects the requirements before writing them: when
  two colliders demand opposite directions on the same edge, neither
  orientation is applied and the conflict is recorded. The output is
  submitted to four gates: absence of conflict, acyclicity, preservation
  of the written colliders, and existence of a consistent DAG extension,
  exposed in `$gates`, `$valid_cpdag` and `$label`. If any gate fails,
  the object is called a second-order partially directed structure.
- [`mc_edge_scores()`](https://evandeilton.github.io/mcglmgraph/reference/mc_edge_scores.md)
  now scores absent edges by the efficient score, which projects out the
  nuisance parameters, instead of dividing `psi_e^2` by `J_ee`. Under
  the working identity the variance is the Schur complement
  `J_ee - J_el J_ll^-1 J_le`; the argument `variability = "empirical"`
  uses the fourth-order variability. The previous statistic remains in
  `stat_unadjusted`, without chi-squared calibration.
- `mc_graph_penalized(precondition = TRUE)` was disabled. The
  preconditioned step is not the proximal operator in the metric of `J`
  and could declare convergence while violating the KKT conditions; on a
  chain with R = 6, n = 80 and lambda = 5 it returned an empty graph
  with a KKT residual close to 30.
- Convergence of the L1 route now requires a double tolerance: relative
  change in the parameters and KKT residual per unit. The fields `kkt`,
  `kkt_unit` and `kkt_tol` accompany the fit, and the
  `kkt_not_certified` state distinguishes stopping by stability from
  stopping by stationarity.
- The graph-layer fit records `score_norm`, `score_norm_unit`, `rcond_J`
  and `stopping_rule`, so that a stopping rule based on parameter change
  can be checked against the score norm it leaves behind.
- [`mc_compute_criterion()`](https://evandeilton.github.io/mcglmgraph/reference/mc_compute_criterion.md)
  no longer silently turns GBIC into BIC when the inversion of the
  sensitivity fails. The fallback is graded (exact solve, documented
  spectral pseudo-inverse, and only then parameter counting) and travels
  in the `criterion_status` attribute, alongside `effective_dimension`.
- [`mc_pc_skeleton()`](https://evandeilton.github.io/mcglmgraph/reference/mc_pc_skeleton.md)
  and
  [`mc_graph_pipeline()`](https://evandeilton.github.io/mcglmgraph/reference/mc_graph_pipeline.md)
  accept `df_residual`, the effective sample size used in the degrees of
  freedom of the Fisher approximation.
- A complete worked example of the workflow is a regression test.

## mcglmgraph 0.3.1

- The stepwise search and the penalty selection explicitly discard
  candidates that did not converge or whose criterion is not finite, and
  record the failures.
- Frequency aggregation distinguishes attempted, valid and failed
  resamples, using the valid ones as the primary denominator and
  preserving a conservative frequency over all attempts.
- DAG initialization deterministically regularizes singular sample
  covariances and records the magnitude of the spectral adjustment.
- Fit objects propagate state, iterations and convergence metadata.

## mcglmgraph 0.3.0

- The `gbic` criterion now explicitly names the working BIC adjusted by
  the Godambe dimension. The name `clbic` remains as a legacy synonym.
- Penalized estimation now accepts an iteration only when the precision
  is positive definite and the penalized objective does not increase.
  The return value includes the convergence state and the objective
  trajectory.
- [`mc_graph_stability()`](https://evandeilton.github.io/mcglmgraph/reference/mc_graph_stability.md)
  is documented as a diagnostic aggregation by selection frequency. The
  function and the former components are preserved for compatibility,
  without attributing formal error control to the procedure.
- The documentation distinguishes the residual Gaussian criterion from a
  node-by-node pseudo-likelihood or a composite likelihood.

## mcglmgraph 0.2.0

### Fixes

- The residual layer explicitly requests
  `residuals(..., type = "pearson")`; the previous call received the
  default raw residuals of `mcglm`.
- The `ordering` argument of the DAG layer is now applied by a
  consistent permutation of the responses, the edges and the fitted
  matrices.
- The GGM line search keeps the last admissible iteration when it fails
  to find a positive definite proposal.
- The PC-stable skeleton examines the frozen neighborhoods of both
  endpoints of each edge. On the optional route, a blacklist prevents
  `bnlearn` from adding edges outside the skeleton; since score-based
  search can remove edges, that route is not described as the
  orientation of a fixed skeleton.
- The documentation distinguishes the methods inherited from the
  first-stage `mcglm` fit from the specific results stored in the
  `$graph` component.

This release aligns the estimators with the MCGLM-Graph method.

### Estimation changes

- DAG layer: the Newton step for the innovation variances now operates
  on the unconstrained reparametrization `xi = log(delta)`, with
  `delta = exp(xi)`. Positive definiteness of `Sigma_b` is guaranteed by
  construction for any `(gamma, xi)`, dispensing with the hard
  projection `pmax(delta, eps)` that introduced a non-differentiable
  update at the boundary. The estimating function and the sensitivity
  are mapped to the `xi` scale by the diagonal Jacobian
  `d(delta)/d(xi) = diag(delta)`, and the reported inference remains on
  the interpretable `(gamma, delta)` scale.

### Documentation

- In the user-facing interpretations, “conditional independence” was
  replaced by “second-order conditional association”, which coincides
  with full conditional independence under joint normality. In Gaussian
  copula models the graphical interpretation belongs to the latent scale
  (`mc_partial_cor`, `mc_test_independence`, package overview).
- `mc_test_independence`: title and label of the test updated to
  “conditional association”.
- CPDAG partial-correlation test: recorded the caveat that Fisher’s z
  transformation is a working approximation under non-Gaussian families,
  with the Wald/delta-method route as a robust alternative.
