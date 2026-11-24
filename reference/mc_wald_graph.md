# Wald test of a linear hypothesis on the graph layer

Tests \\H_0:\mathbf L\boldsymbol\theta=\mathbf c\\ on the parameters of
the graph layer by the generalized Wald statistic, \$\$W=(\mathbf
L\hat{\boldsymbol\theta}-\mathbf c)^\top \left\[\mathbf
L\\\widehat{\mathrm{Var}}(\hat{\boldsymbol\theta})\\ \mathbf
L^\top\right\]^{-1}(\mathbf L\hat{\boldsymbol\theta}-\mathbf c),\$\$
referred to a chi-squared distribution with as many degrees of freedom
as the rank of \\\mathbf L\\. The covariance
\\\widehat{\mathrm{Var}}(\hat{\boldsymbol\theta})\\ is any of the three
forms of
[`mc_vcov_graph`](https://evandeilton.github.io/mcglmgraph/reference/mc_vcov_graph.md):
model-based, sandwich, or sandwich corrected for the first stage.

## Usage

``` r
mc_wald_graph(
  fit,
  L,
  rhs = 0,
  type = c("sandwich", "model", "two_stage"),
  cluster = NULL,
  small_sample = FALSE,
  method = c("auto", "shortcut", "full")
)
```

## Arguments

- fit:

  an object of class `"mcglm_graph"`.

- L:

  the contrast: a numeric matrix with one column per graph parameter, in
  the order of `coef(fit, type = "graph")` (or with those names as
  column names, in any order); a numeric vector, read as one row; or a
  character vector of parameter names, which selects them.

- rhs:

  the right-hand side \\\mathbf c\\, a single number or one value per
  row of `L`.

- type:

  the covariance: `"sandwich"` (default), `"model"` or `"two_stage"`, as
  in
  [`mc_vcov_graph`](https://evandeilton.github.io/mcglmgraph/reference/mc_vcov_graph.md).

- cluster:

  optional vector with one identifier per observational row. With
  `type = "sandwich"` or `"two_stage"`, the contributions are summed
  within clusters before the variability is formed, which is the form to
  use when the unit of dependence is not the row of the data matrix. It
  must have no missing value and at least two distinct values; the
  default `NULL` takes the rows as the units. It is refused with
  `type = "model"`.

- small_sample:

  logical, default `FALSE`; multiply the clustered variability by
  \\G/(G-1)\\, with \\G\\ the number of clusters. It has no effect
  without `cluster` and is refused with `type = "model"`.

- method:

  for `type = "two_stage"`, the route of the correction: `"shortcut"`,
  `"full"`, or `"auto"` (default). `"auto"` takes the full form when it
  is available and either a response has a variance function other than
  `"constant"` or the fitted DAG has an immorality; in every other case
  it takes the shortcut when it applies and the full form when it does
  not. The full form rebuilds the first stage and is slower. It has no
  meaning for the other two types.

## Value

An object of class `"htest"`: the Wald `statistic` \\W\\ (named
`"X-squared"`); `parameter`, the degrees of freedom `df`, equal to the
rank of `L`; the `p.value`; the `estimate` \\\mathbf
L\hat{\boldsymbol\theta}\\ and the `null.value` \\\mathbf c\\, both
named by the rows of `L`; the `method` and `data.name` labels; the
numeric matrix `contrast` used, with one row per row of `L` and one
column per graph parameter; and the covariance `type`.

## Details

Let \\q=\mathrm{rank}(\mathbf L)\\. Under \\H_0\\, with the number of
responses fixed and the number of units \\n\to\infty\\, the asymptotic
normality of \\\sqrt n(\hat{\boldsymbol\theta}-\boldsymbol\theta_0)\\
with covariance \\\mathbf J^{-1}\mathbf V\mathbf J^{-\top}\\ gives
\\\sqrt n(\mathbf L\hat{\boldsymbol\theta}-\mathbf c)\xrightarrow{d}
\mathcal N(\mathbf 0,\boldsymbol\Xi)\\, with \\\boldsymbol\Xi=\mathbf
L\mathbf J^{-1}\mathbf V\mathbf J^{-\top} \mathbf L^\top\\ positive
definite when \\\mathbf J\\ is nonsingular and \\\mathbf V\\ is positive
definite; \\W\\ then converges to \\\chi^2_q\\ once \\\boldsymbol\Xi\\
is replaced by a consistent estimator. For a single edge,
\\W=\hat\theta_e^2/ \widehat{\mathrm{Var}}(\hat\theta_e)\\ and \\q=1\\.
The sensitivity \\\mathbf J\\ is the expected one in the GGM layer,
where it coincides with the curvature, and the observed curvature of the
working criterion in the DAG layer. With `cluster`, the contributions
are summed within independent clusters before the variability is formed,
and the asymptotic statement is in the number of clusters.

## Scale

The statistic can be written with a per-unit curvature and variability,
\\W=n(\mathbf L\hat{\boldsymbol\theta}-\mathbf c)^\top\[\mathbf L
\mathbf H^{-1}\mathbf V\mathbf H^{-\top}\mathbf L^\top\]^{-1} (\mathbf
L\hat{\boldsymbol\theta}-\mathbf c)\\, with \\\mathbf H\\ the observed
curvature of the graph block, equal to the expected sensitivity
\\\mathbf J\\ in the GGM layer. The covariance returned by
[`mc_vcov_graph`](https://evandeilton.github.io/mcglmgraph/reference/mc_vcov_graph.md)
is already on the scale of \\\hat{\boldsymbol\theta}\\, \\\mathbf
H_n^{-1}(n\mathbf V) \mathbf H_n^{-1}=n^{-1}\mathbf H^{-1}\mathbf
V\mathbf H^{-1}\\, so the two expressions are the same number.

## A contrast that is not of full rank

Redundant rows add no restriction. The statistic is then formed on the
range of \\\mathbf L\\\widehat{\mathrm{Var}}\\\mathbf L^\top\\ and the
degrees of freedom are the rank of \\\mathbf L\\; a right-hand side that
the redundant rows contradict is refused. The statistic is invariant to
any invertible recombination of the rows, applied to \\\mathbf L\\ and
\\\mathbf c\\ alike.

## A contrast without variability

A hypothesis is refused when the covariance of the contrast is singular
on its range, judged against the model-based variance of the same
contrast. That happens with `type = "two_stage"` for a parameter the
first stage fixes exactly, such as the innovation variance of a root
node of a DAG with Gaussian intercept-only margins, whose corrected
variance is zero. A contrast counts as having no variability when its
variance is at most \\10^{-6}\\ of its model-based variance, on the
scale of the model-based standard deviations; the margin covers the
numerical residue that the full form leaves where the exact variance is
zero.

## What the reference distribution requires

The chi-squared reference holds for a hypothesis fixed before looking at
the data, under the regularity of the second-order estimating functions
and, with `type = "model"`, the working identity \\\mathbf V=\mathbf
J\\. Applied to parameters of a support chosen by a search over graphs,
it is descriptive. Hypotheses on the first stage, the regression and
dispersion parameters, are outside this function, which acts on the
graph layer only; the htmcglm package tests those.

## References

Godambe, V. P. (1960). An optimum property of regular maximum likelihood
estimation. *Annals of Mathematical Statistics* 31(4), 1208–1211.
[doi:10.1214/aoms/1177705693](https://doi.org/10.1214/aoms/1177705693)

de Freitas, L. A. C. and Bonat, W. H. (2022). Hypothesis tests for
multiple responses regression models in R: the htmcglm package. arXiv
preprint arXiv:2208.01568.
[doi:10.48550/arXiv.2208.01568](https://doi.org/10.48550/arXiv.2208.01568)

## See also

[`mc_test_independence`](https://evandeilton.github.io/mcglmgraph/reference/mc_test_independence.md)
and
[`mc_test_adequacy`](https://evandeilton.github.io/mcglmgraph/reference/mc_test_adequacy.md),
which are this test with a particular contrast, and
[`mc_vcov_graph`](https://evandeilton.github.io/mcglmgraph/reference/mc_vcov_graph.md)
for the covariance.

Other inference on the graph layer:
[`mc_test_adequacy()`](https://evandeilton.github.io/mcglmgraph/reference/mc_test_adequacy.md),
[`mc_test_independence()`](https://evandeilton.github.io/mcglmgraph/reference/mc_test_independence.md),
[`mc_two_stage_available()`](https://evandeilton.github.io/mcglmgraph/reference/mc_two_stage_available.md),
[`mc_vcov_graph()`](https://evandeilton.github.io/mcglmgraph/reference/mc_vcov_graph.md)

## Examples

``` r
# \donttest{
set.seed(1); n <- 300
K <- diag(3) * 2; K[1, 2] <- K[2, 1] <- -0.7
Y <- matrix(rnorm(n * 3), n) %*% chol(solve(K))
d <- data.frame(y1 = Y[, 1], y2 = Y[, 2], y3 = Y[, 3], x = rnorm(n))
f <- mcglm_graph(list(y1 ~ x, y2 ~ x, y3 ~ x), data = d)
#> Automatic initial values selected.
mc_wald_graph(f, c("eta_13", "eta_23"))          # two edges at once
#> 
#>  Wald test of 2 linear restrictions on the graph layer (sandwich
#>  covariance)
#> 
#> data:  GGM layer with 6 parameters
#> X-squared = 1.0683, df = 2, p-value = 0.5862
#> sample estimates:
#> [1] -0.02798634 -0.04495199
#> 
mc_wald_graph(f, c(1, -1, 0, 0, 0, 0))           # eta_11 = eta_22
#> 
#>  Wald test of 1 linear restriction on the graph layer (sandwich
#>  covariance)
#> 
#> data:  GGM layer with 6 parameters
#> X-squared = 9.5848e-05, df = 1, p-value = 0.9922
#> sample estimates:
#> [1] -0.001224771
#> 
# }
```
