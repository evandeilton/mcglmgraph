# Wald test of between-response conditional association

Tests the null hypothesis that all edge parameters of the fitted graph
layer are zero, under which the between-response covariance
\\\boldsymbol\Sigma_b\\ of the Pearson residuals is diagonal: the
absence of second-order conditional association among the responses
given the covariates. Under joint normality of the analyzed residuals
this coincides with mutual conditional independence; in Gaussian-copula
models the same reading requires a precision matrix fitted on the latent
Gaussian scale.

## Usage

``` r
mc_test_independence(
  fit,
  robust = FALSE,
  cluster = NULL,
  small_sample = FALSE,
  type = NULL,
  method = c("auto", "shortcut", "full")
)
```

## Arguments

- fit:

  an object of class `"mcglm_graph"` with at least one edge.

- robust:

  logical, default `FALSE`; if `TRUE`, the test uses the robust sandwich
  covariance \\\mathbf H_n^{-1}\mathbf V_n\mathbf H_n^{-1}\\, with
  \\\mathbf H_n\\ the observed curvature of the graph block, built from
  the empirical fourth-moment variability, instead of the model-based
  \\\mathbf J_n^{-1}\\ (which assumes \\\mathbf V=\mathbf J\\). It is
  equivalent to `type = "sandwich"`, which also offers `"two_stage"`.

- cluster:

  optional vector with one identifier per observational row, identifying
  independent clusters. Requires an empirical covariance
  (`robust = TRUE` or an empirical `type`): Pearson contributions are
  then summed within clusters before the empirical variability is
  formed.

- small_sample:

  logical, default `FALSE`; multiply the cluster variability by
  \\G/(G-1)\\, where \\G\\ is the number of clusters. Requires an
  empirical covariance and has no effect without `cluster`.

- type:

  optional covariance, `"model"`, `"sandwich"` or `"two_stage"`, as in
  [`mc_vcov_graph`](https://evandeilton.github.io/mcglmgraph/reference/mc_vcov_graph.md).
  When given it takes precedence over `robust`, which must then not
  contradict it. The default `NULL` uses `"sandwich"` if `robust` is
  `TRUE` and `"model"` otherwise.

- method:

  for `type = "two_stage"`, the route of the correction, `"auto"`
  (default), `"shortcut"` or `"full"`, as in
  [`mc_vcov_graph`](https://evandeilton.github.io/mcglmgraph/reference/mc_vcov_graph.md).
  It has no meaning for the other two types and is refused with them
  when it is not `"auto"`.

## Value

An object of class `"htest"`: the Wald `statistic`, the degrees of
freedom in `parameter` (the number of edges), the `p.value`, and the
`method` and `data.name` labels.

## Details

It is
[`mc_wald_graph`](https://evandeilton.github.io/mcglmgraph/reference/mc_wald_graph.md)
with the contrast \\\mathbf L\\ that selects the \\q\\ edge parameters
\\\boldsymbol\theta_E\\ (the off-diagonal precision parameters of the
GGM layer, the structural coefficients of the DAG layer) and \\\mathbf
c=\mathbf 0\\:
\$\$W=\hat{\boldsymbol\theta}\_E^\top\\\widehat{\mathrm{Var}}(
\hat{\boldsymbol\theta}\_E)^{-1}\\\hat{\boldsymbol\theta}\_E,\$\$ with
the covariance of
[`mc_vcov_graph`](https://evandeilton.github.io/mcglmgraph/reference/mc_vcov_graph.md)
selected by `type`. Under the null, with the number of responses fixed
and \\n\to\infty\\, \\W\\ is asymptotically \\\chi^2_q\\, \\q\\ being
the number of edges; the p-value is the upper tail of that distribution.
The default covariance is the model-based one, for which that reference
also requires the working identity \\\mathbf V=\mathbf J\\ (a correctly
specified second moment and Gaussian fourth moments);
`type = "sandwich"` does not use it, and `type = "two_stage"` also
accounts for the estimation of the first stage, by the route `method`
selects. The edges tested are those of the support of the fit, and the
graph is treated as fixed: for a support chosen by a search over graphs,
the result is descriptive. It is a standalone function so that the
inherited `anova.mcglm` method remains available on the object. The
function refuses fits whose numerical diagnostics mark graph-layer
inference as uncertified.

## References

Godambe, V. P. (1960). An optimum property of regular maximum likelihood
estimation. *Annals of Mathematical Statistics* 31(4), 1208–1211.
[doi:10.1214/aoms/1177705693](https://doi.org/10.1214/aoms/1177705693)

## See also

[`mc_wald_graph`](https://evandeilton.github.io/mcglmgraph/reference/mc_wald_graph.md)
for other linear hypotheses,
[`mc_test_adequacy`](https://evandeilton.github.io/mcglmgraph/reference/mc_test_adequacy.md)
for the adequacy of a selected support, and
[`mc_vcov_graph`](https://evandeilton.github.io/mcglmgraph/reference/mc_vcov_graph.md)
for the covariance forms.

Other inference on the graph layer:
[`mc_test_adequacy()`](https://evandeilton.github.io/mcglmgraph/reference/mc_test_adequacy.md),
[`mc_two_stage_available()`](https://evandeilton.github.io/mcglmgraph/reference/mc_two_stage_available.md),
[`mc_vcov_graph()`](https://evandeilton.github.io/mcglmgraph/reference/mc_vcov_graph.md),
[`mc_wald_graph()`](https://evandeilton.github.io/mcglmgraph/reference/mc_wald_graph.md)

## Examples

``` r
# \donttest{
set.seed(1); n <- 300
K <- diag(3) * 2; K[1,2] <- K[2,1] <- -0.7
Y <- matrix(rnorm(n * 3), n) %*% chol(solve(K))
d <- data.frame(y1 = Y[,1], y2 = Y[,2], y3 = Y[,3], x = rnorm(n))
f <- mcglm_graph(list(y1~x, y2~x, y3~x), data = d)
#> Automatic initial values selected.
mc_test_independence(f)
#> 
#>  Wald test of between-response conditional association
#> 
#> data:  GGM layer with 3 edges
#> X-squared = 33.437, df = 3, p-value = 2.604e-07
#> 
mc_test_independence(f, robust = TRUE)
#> 
#>  Wald test of between-response conditional association
#> 
#> data:  GGM layer with 3 edges
#> X-squared = 28.885, df = 3, p-value = 2.368e-06
#> 
mc_test_independence(f, type = "two_stage")
#> 
#>  Wald test of between-response conditional association
#> 
#> data:  GGM layer with 3 edges
#> X-squared = 28.282, df = 3, p-value = 3.169e-06
#> 
# }
```
