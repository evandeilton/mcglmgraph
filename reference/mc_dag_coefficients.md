# Scale-free coefficients of a fitted directed layer

Reports, with standard errors and intervals, the quantities of the
directed layer that do not depend on the scale on which the residual
dependence is measured. Under a diagonal rescaling
\\\boldsymbol\Sigma_b\mapsto\mathbf A\boldsymbol\Sigma_b\mathbf A\\ the
raw coefficients move as \\\gamma\_{jk}\mapsto\gamma\_{jk}a_j/a_k\\ and
the innovation variances as \\\delta_j\mapsto a_j^2\delta_j\\, whereas
\$\$\gamma^\*\_{jk}=\gamma\_{jk}\left(\frac{\Sigma\_{b,kk}}{\Sigma\_{b,jj}}
\right)^{1/2},\qquad R^2\_{j\cdot
P_j}=1-\frac{\delta_j}{\Sigma\_{b,jj}}\$\$ are invariant. These are the
quantities to report: they play for the directed layer the role that
partial correlations play for the undirected one.

## Usage

``` r
mc_dag_coefficients(
  fit,
  variability = c("model", "empirical", "two_stage"),
  level = 0.95,
  interval = c("delta", "fixed_factor"),
  cluster = NULL,
  small_sample = FALSE,
  r2_interval = c("wald", "logit")
)

# S3 method for class 'mcg_dag_coef'
print(x, ...)
```

## Arguments

- fit:

  an object of class `"mcglm_graph"` fitted with the DAG layer.

- variability:

  `"model"` for \\\mathbf J_n^{-1}\\, `"empirical"` for the sandwich, or
  `"two_stage"` for the sandwich corrected for the estimation of the
  residuals.

- level:

  confidence level of the intervals.

- interval:

  `"delta"` (default) or `"fixed_factor"`, the interval for
  \\\gamma^\*\\ as described above.

- cluster, small_sample:

  passed to
  [`mc_vcov_graph`](https://evandeilton.github.io/mcglmgraph/reference/mc_vcov_graph.md)
  for the empirical forms.

- r2_interval:

  `"wald"` (default), symmetric on the scale of \\R^2\\, or `"logit"`,
  symmetric on the scale of \\\log\\R^2/(1-R^2)\\\\ and mapped back, so
  that the bounds stay inside \\(0, 1)\\. A response without parents has
  \\R^2=0\\ by construction, and its interval is that point.

- x:

  an object of class `"mcg_dag_coef"`.

- ...:

  ignored.

## Value

An object of class `"mcg_dag_coef"`: a list with `edges`, one row per
directed edge with the raw and standardized coefficient, the standard
error of each, the Wald statistic of the raw coefficient, its p-value
and the two intervals, and `nodes`, one row per response with the
innovation variance and its standard error, and \\R^2\_{j\cdot P_j}\\
with its standard error and interval. The gradients of \\\gamma^\*\\ and
of \\R^2\\ travel in the attributes `"gradient_std"` and
`"gradient_r2"`. The list also records `variability`, `level`,
`interval` and `r2_interval`. The `print` method shows the main columns
and returns `x` invisibly.

## What the coefficients are

\\\gamma\_{jk}\\ is the coefficient of response \\k\\ in the linear
projection of the residual component of response \\j\\ on those of its
parents, under the ordering of the fit, and \\\delta_j\\ is the variance
left by that projection. The factorization does not by itself define a
causal model. Reading a coefficient as the effect of an intervention
requires a linear structural model with independent exogenous
innovations, the correct ordering, no unobserved confounding and
equations that are invariant under intervention. The fit checks none of
these.

## When the first stage is such a rescaling

The first stage acts on the residual covariance, to first order, as a
diagonal rescaling under conditions (b) and (c) of the section "What the
correction reaches" of
[`mc_vcov_graph`](https://evandeilton.github.io/mcglmgraph/reference/mc_vcov_graph.md):
a correctly specified mean, and a conditional residual covariance that
does not depend on the covariates. Under them the direction in which the
first stage moves \\\boldsymbol\Sigma_b\\ is
\\\boldsymbol\Lambda\boldsymbol\Sigma_b+
\boldsymbol\Sigma_b\boldsymbol\Lambda\\ with \\\boldsymbol\Lambda\\
diagonal, whatever the parametrization, and the quantities above do not
see it. The two conditions are sufficient: with the constant variance
function (c) is not needed, and with a variance function that depends on
the mean a residual correlation that changes with the covariates moves
\\\boldsymbol\Sigma_b\\ in a direction the standardized coefficients do
see. The finer characterization of the section "When the shortcut is
right" of
[`mc_vcov_graph`](https://evandeilton.github.io/mcglmgraph/reference/mc_vcov_graph.md)
was established for the undirected layer; for this one the sufficient
conditions are the ones to rely on. The package checks none of them.

## What the standard errors cover

They come from
[`mc_vcov_graph`](https://evandeilton.github.io/mcglmgraph/reference/mc_vcov_graph.md)
and refer to the second stage conditionally on the estimated residuals,
unless `variability = "two_stage"`. The empirical forms use the observed
Hessian of the working criterion; the model form uses the Pearson
sensitivity and requires correct second-moment specification and its
Gaussian working variability identity. The Wald statistic of
\\\gamma\_{jk}=0\\ is invariant under the diagonal rescaling, because
the coefficient and its standard error are both multiplied by
\\a_j/a_k\\, so under (b) and (c) the test is adequate to first order
where the raw coefficient is not.

## Intervals for the scale-free quantities

Both \\\gamma^\*\_{jk}\\ and \\R^2\_{j\cdot P_j}\\ are functions of
\\(\boldsymbol\gamma,\boldsymbol\delta)\\, through the coefficient
itself and through \\\boldsymbol\Sigma_b=\mathbf T^{-1}\mathbf D\mathbf
T^{-\top}\\, so their variance follows by the delta method with the
analytical gradients
\$\$\frac{\partial\gamma^\*\_{jk}}{\partial\theta_d}=
\frac{\partial\gamma\_{jk}}{\partial\theta_d}\\f\_{jk}+
\frac{\gamma^\*\_{jk}}{2}\left(\frac{\dot\Sigma\_{b,d,kk}}{\Sigma\_{b,kk}}-
\frac{\dot\Sigma\_{b,d,jj}}{\Sigma\_{b,jj}}\right),\qquad \frac{\partial
R^2\_{j\cdot P_j}}{\partial\theta_d}=
-\frac{1}{\Sigma\_{b,jj}}\frac{\partial\delta_j}{\partial\theta_d}+
\frac{\delta_j\\\dot\Sigma\_{b,d,jj}}{\Sigma\_{b,jj}^2},\$\$ with
\\f\_{jk}=(\Sigma\_{b,kk}/\Sigma\_{b,jj})^{1/2}\\ and
\\\dot{\boldsymbol\Sigma}\_{b,d}=\partial\boldsymbol\Sigma_b/
\partial\theta_d\\. That is `interval = "delta"`, the default. With
`interval = "fixed_factor"` the interval for \\\gamma^\*\\ is the one
for \\\gamma\\ multiplied by \\f\_{jk}\\ as if the factor were known,
which leaves out the uncertainty of the variances in it. The interval
for \\R^2\\ is always by the delta method, since it has no fixed-factor
analogue.

## See also

[`mc_vcov_graph`](https://evandeilton.github.io/mcglmgraph/reference/mc_vcov_graph.md),
[`mc_dag_unit_diagonal`](https://evandeilton.github.io/mcglmgraph/reference/mc_dag_unit_diagonal.md)

Other directed layer:
[`mc_build_D_matrix()`](https://evandeilton.github.io/mcglmgraph/reference/mc_build_D_matrix.md),
[`mc_build_T_matrix()`](https://evandeilton.github.io/mcglmgraph/reference/mc_build_T_matrix.md),
[`mc_cholesky_derivatives()`](https://evandeilton.github.io/mcglmgraph/reference/mc_cholesky_derivatives.md),
[`mc_cholesky_to_covariance()`](https://evandeilton.github.io/mcglmgraph/reference/mc_cholesky_to_covariance.md),
[`mc_cholesky_to_precision()`](https://evandeilton.github.io/mcglmgraph/reference/mc_cholesky_to_precision.md),
[`mc_dag_closed_root()`](https://evandeilton.github.io/mcglmgraph/reference/mc_dag_closed_root.md),
[`mc_dag_sensitivity()`](https://evandeilton.github.io/mcglmgraph/reference/mc_dag_sensitivity.md),
[`mc_dag_unit_diagonal()`](https://evandeilton.github.io/mcglmgraph/reference/mc_dag_unit_diagonal.md)

## Examples

``` r
# \donttest{
set.seed(4); n <- 400; x <- rnorm(n)
y1 <- rnorm(n); y2 <- 0.6 * y1 + rnorm(n); y3 <- 0.4 * y2 + rnorm(n)
d <- data.frame(y1 = y1, y2 = y2, y3 = y3, x = x)
f <- mcglm_graph(list(y1 ~ x, y2 ~ x, y3 ~ x),
                 between_graph = mc_between_graph("dag"), data = d)
#> Automatic initial values selected.
mc_dag_coefficients(f)
#> Directed layer (model variability, level 0.95)
#> 
#> Edges (parent -> child):
#>   parent child    gamma gamma_std      se  se_std statistic   p_value
#> 1     y1    y2 0.538120  0.538120 0.04214 0.03552   12.7688 2.450e-37
#> 2     y1    y3 0.005664  0.005664 0.05431 0.05430    0.1043 9.169e-01
#> 3     y2    y3 0.399375  0.399375 0.05431 0.05114    7.3543 1.919e-13
#> 
#> Nodes:
#>   response n_parents  delta se_delta r_squared se_r_squared
#> 1       y1         0 0.9950  0.07036    0.0000      0.00000
#> 2       y2         1 0.7069  0.04998    0.2896      0.03823
#> 3       y3         2 0.8338  0.05896    0.1620      0.03373
# }
```
