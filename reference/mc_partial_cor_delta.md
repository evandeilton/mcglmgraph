# Delta-method inference for the between-response partial correlations

Attaches standard errors and confidence intervals to the partial
correlations of
[`mc_partial_cor`](https://evandeilton.github.io/mcglmgraph/reference/mc_partial_cor.md)
by the delta method. The gradient of \\\rho\_{ij\mid\mathrm{rest}} =
-K\_{ij}/\sqrt{K\_{ii}K\_{jj}}\\ with respect to the covariance
parameters \\\theta_d\\ is \$\$\[\mathbf g\_{ij}\]\_d =
-\frac{\[\dot{\mathbf K}\_d\]\_{ij}}{\sqrt{K\_{ii}K\_{jj}}} -
\frac{\rho\_{ij\mid\mathrm{rest}}}{2}\left\\ \frac{\[\dot{\mathbf
K}\_d\]\_{ii}}{K\_{ii}} + \frac{\[\dot{\mathbf
K}\_d\]\_{jj}}{K\_{jj}}\right\\,\$\$ and the variance follows as
\\\mathbf g\_{ij}^\top
\widehat{\mathrm{Var}}(\hat{\boldsymbol\theta})\\\mathbf g\_{ij}\\ with
\\\dot{\mathbf K}\_d = \partial\mathbf K_b/\partial\theta_d\\. In the
GGM layer \\\dot{\mathbf K}\_d\\ is the canonical base \\\mathbf Z_d\\;
the implementation obtains it from \\\dot{\mathbf K}\_d = -\mathbf K_b
\dot{\boldsymbol\Sigma}\_d \mathbf K_b\\, which holds for any
parametrization, so the same expression also covers the DAG layer, where
the parameters are \\(\boldsymbol\gamma,\boldsymbol\delta)\\.

## Usage

``` r
mc_partial_cor_delta(
  fit,
  variability = c("model", "empirical", "two_stage"),
  cluster = NULL,
  small_sample = FALSE,
  level = 0.95,
  interval = c("wald", "fisher")
)
```

## Arguments

- fit:

  an object of class `"mcglm_graph"`.

- variability:

  `"model"` for the model-based covariance \\\mathbf J_n^{-1}\\,
  `"empirical"` for the sandwich form \\\mathbf H_n^{-1}(n\hat{\mathbf
  V})\mathbf H_n^{-1}\\, with \\\mathbf H_n\\ the observed curvature of
  the graph block, or `"two_stage"` for the sandwich corrected for the
  estimation of the residuals. When the correction takes the shortcut,
  the two empirical forms give the partial correlations exactly the same
  variance, since their gradients annihilate the shortcut's cross
  sensitivity algebraically; with the full form they agree to first
  order under the conditions of the section "What the correction
  reaches" of
  [`mc_vcov_graph`](https://evandeilton.github.io/mcglmgraph/reference/mc_vcov_graph.md).
  Whether that common variance is right depends on the same conditions:
  a residual correlation that moves with the covariates puts both off,
  and only the full form,
  [`mc_vcov_graph`](https://evandeilton.github.io/mcglmgraph/reference/mc_vcov_graph.md)`(type = "two_stage", method = "full")`,
  accounts for it. The model-based form agrees with the empirical ones
  only when, in addition, \\\mathbf V=\mathbf J\\, which fails outside
  the Gaussian fourth-moment regime.

- cluster:

  optional vector with one identifier per observational row, honored by
  the two empirical forms (`"empirical"` and `"two_stage"`).
  Contributions are summed within clusters before the variability is
  formed, which is what to use when the unit of dependence is not the
  row of the data matrix.

- small_sample:

  logical; multiply the clustered variability by \\G/(G-1)\\, with \\G\\
  the number of clusters (default `FALSE`).

- level:

  confidence level of the interval, a single number in \\(0, 1)\\
  (default 0.95).

- interval:

  `"wald"` (default), symmetric on the scale of the partial correlation,
  or `"fisher"`, symmetric on the scale of
  \\\operatorname{atanh}(\rho)\\, where the standard error is
  \\\mathrm{se}/(1-\rho^2)\\, and mapped back by \\\tanh\\, so that the
  bounds stay inside \\(-1, 1)\\. The standard error, the statistic and
  the p-value are the same under both.

## Value

A data frame with one row per pair of responses (\\q = R(R-1)/2\\ rows)
and columns `response1`, `response2`, `present` (logical support
indicator), `pcor` (the partial correlation), `se` (standard error),
`statistic` (Wald statistic `pcor / se`, `NA` when the standard error is
exactly zero), `p_value` (two-sided, normal reference) and `lower`,
`upper` (interval bounds). The rows follow the pairs \\(i, j)\\, \\i \<
j\\, in lexicographic order. The attribute `"gradient"` holds the matrix
of gradients \\\mathbf g\_{ij}\\, with one row per pair and one column
per graph parameter, and the attributes `"variability"`, `"level"` and
`"interval"` record the choices used.

## Scope of the inference

With `variability = "model"` or `"empirical"` the variance refers to the
second stage conditionally on the estimated residuals: it propagates the
uncertainty of the covariance parameters, not the uncertainty inherited
from the first-stage fit, which `"two_stage"` adds. With
`variability = "model"` it also assumes the working identity \\\mathbf V
= \mathbf J\\; outside the Gaussian fourth-moment regime the empirical
option is the appropriate one. Applied to a support chosen by a search
over graphs, the interval remains descriptive, since the selection is
not incorporated. The function refuses fits whose numerical diagnostics
explicitly mark graph-layer inference as uncertified.

## Structural zeros

In the GGM layer, a pair outside the specified precision support is zero
by construction: no parameter moves it, the gradient vanishes and the
reported standard error is zero up to rounding error, so that the
statistic of such a row is 0 or `NA` and carries no information. In the
DAG layer, structural zeros of the precision are instead determined by
the moral graph. In particular, two non-adjacent parents of a common
child can have a nonzero precision entry even though no arrow joins
them. Therefore `present` denotes adjacency in the specified GGM support
or in the moral graph of the specified DAG. Rows with `present = FALSE`
are structural constraints and should not be read as precise estimates
of a null association; a `TRUE` value says only that the entry is
structurally allowed, since parameter cancellation can still make it
zero.

## References

Lauritzen, S. L. (1996). *Graphical Models*. Oxford University Press.

## See also

[`mc_partial_cor`](https://evandeilton.github.io/mcglmgraph/reference/mc_partial_cor.md)
for the point estimates,
[`mc_vcov_graph`](https://evandeilton.github.io/mcglmgraph/reference/mc_vcov_graph.md)
for the covariance of the graph parameters and
[`mc_test_independence`](https://evandeilton.github.io/mcglmgraph/reference/mc_test_independence.md)
for the joint Wald test.

Other graph interpretation:
[`mc_get_graph()`](https://evandeilton.github.io/mcglmgraph/reference/mc_get_graph.md),
[`mc_partial_cor()`](https://evandeilton.github.io/mcglmgraph/reference/mc_partial_cor.md),
[`mc_plot_network()`](https://evandeilton.github.io/mcglmgraph/reference/mc_plot_network.md)

## Examples

``` r
set.seed(1); n <- 300
K <- diag(3) * 2; K[1, 2] <- K[2, 1] <- -0.7
Y <- matrix(rnorm(n * 3), n) %*% chol(solve(K))
d <- data.frame(y1 = Y[, 1], y2 = Y[, 2], y3 = Y[, 3], x = rnorm(n))
f <- mcglm_graph(list(y1~x, y2~x, y3~x), data = d)
#> Automatic initial values selected.
mc_partial_cor_delta(f)
#>   response1 response2 present       pcor         se statistic      p_value
#> 1        y1        y2    TRUE 0.34500759 0.05086281 6.7831008 1.176235e-11
#> 2        y1        y3    TRUE 0.02602162 0.05769593 0.4510130 6.519802e-01
#> 3        y2        y3    TRUE 0.04177389 0.05763428 0.7248098 4.685687e-01
#>         lower     upper
#> 1  0.24531830 0.4446969
#> 2 -0.08706033 0.1391036
#> 3 -0.07118722 0.1547350
mc_partial_cor_delta(f, variability = "empirical")
#>   response1 response2 present       pcor         se statistic      p_value
#> 1        y1        y2    TRUE 0.34500759 0.05402660 6.3858832 1.704110e-10
#> 2        y1        y3    TRUE 0.02602162 0.06203153 0.4194902 6.748579e-01
#> 3        y2        y3    TRUE 0.04177389 0.05467196 0.7640826 4.448180e-01
#>         lower     upper
#> 1  0.23911739 0.4508978
#> 2 -0.09555795 0.1476012
#> 3 -0.06538118 0.1489290
```
