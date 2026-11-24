# Introduction to mcglmgraph

## Overview

`mcglmgraph` extends Multivariate Covariance Generalized Linear Models
(MCGLM) by structuring the residual between-response covariance
$`\Sigma_b`$ with a graph. Three layers are available: the Gaussian
graphical model (GGM), which parametrizes the precision matrix
$`K_b = \Sigma_b^{-1}`$ so that its zeros are null partial correlations;
the directed acyclic graph (DAG), which uses the modified Cholesky (LDL)
decomposition under a fixed ordering of the responses; and the completed
partially directed acyclic graph (CPDAG), the equivalence class of DAGs
that partial-correlation tests can identify, learned by
structure-discovery routines. Numerically intensive operations run in
C++ via `RcppArmadillo`.

Estimation proceeds in two stages. The first stage fits the marginal
MCGLM (mean and dispersion blocks) and extracts the standardized Pearson
residuals. The second stage forms their sample between-response
covariance and fits the graph layer to this second-order object with the
Pearson estimating functions of the MCGLM, without a joint likelihood.
Graph edges therefore describe residual conditional association among
responses, not raw marginal association. In the GGM layer, a missing
edge means absence of second-order conditional association; it coincides
with conditional independence only under joint normality of the analyzed
residuals. In a Gaussian-copula model, that reading applies to a
precision matrix fitted on the latent Gaussian scale; the
Pearson-residual route does not perform that latent transformation.

``` r

library(mcglmgraph)
```

## A worked example: recovering a chain graph

We simulate four responses whose true second-order conditional
association graph is the chain $`1 - 2 - 3 - 4`$, add a covariate to two
of them, and fit the GGM layer.

``` r

set.seed(1)
n <- 300
K <- diag(4) * 2
K[1, 2] <- K[2, 1] <- -0.6
K[2, 3] <- K[3, 2] <- -0.6
K[3, 4] <- K[4, 3] <- -0.6
Y <- MASS::mvrnorm(n, rep(0, 4), solve(K))
x <- rnorm(n)
dat <- data.frame(y1 = Y[, 1] + 0.5 * x, y2 = Y[, 2],
                  y3 = Y[, 3], y4 = Y[, 4] - 0.3 * x, x = x)

fit <- mcglm_graph(list(y1 ~ x, y2 ~ x, y3 ~ x, y4 ~ x),
                   between_graph = mc_between_graph("ggm"), data = dat)
#> Automatic initial values selected.
fit
#> Multivariate Covariance GLM with graph-structured between-response covariance
#> Layer: GGM | responses: 4 | n: 300 | edges: 6
#> Between-response engine converged: TRUE (1 iterations)
#> 
#> Between-response covariance (Sigma_b):
#>         y1     y2      y3      y4
#> y1  0.9933 0.2606 -0.0033 -0.0142
#> y2  0.2606 0.9933  0.3174  0.0559
#> y3 -0.0033 0.3174  0.9933  0.2311
#> y4 -0.0142 0.0559  0.2311  0.9933
#> 
#> This is also a valid 'mcglm' object: use coef(), vcov(), gof(), ...
```

The partial correlations display the second-order conditional
association structure. The precision matrix used in the simulation gives
a partial correlation of 0.30 to each pair of adjacent responses and
zero to the other pairs; the estimates for adjacent pairs are close to
that value and those for non-adjacent pairs are close to zero.

``` r

round(mc_partial_cor(fit), 3)
#>        y1     y2     y3     y4
#> y1  1.000  0.278 -0.091 -0.009
#> y2  0.278  1.000  0.327 -0.016
#> y3 -0.091  0.327  1.000  0.225
#> y4 -0.009 -0.016  0.225  1.000
```

The regression coefficients are inherited from the marginal MCGLM fit.
The graph layer leaves the mean target fixed and re-expresses the
residual between-response covariance. The standard errors of the graph
parameters are conditional on the first-stage fit unless a two-stage
variance is requested (`mc_vcov_graph(type = "two_stage")`).

How much the first stage matters depends on the parameter. Suppose that
the Pearson residual of each response is $`(y - \mu)/\varsigma`$, with
mean $`\mu`$ and scale $`\varsigma`$ depending only on the first-stage
parameters and the covariates, that the mean is correctly specified, and
that the conditional covariance of the residuals given the covariates is
the same $`\Sigma_b`$ for every unit. In the GGM layer, with a support
that contains the nonzero off-diagonal entries of the true precision
matrix, estimating the first stage then acts on the residual covariance,
to first order, as a diagonal rescaling. Partial correlations and zeros
of the precision matrix are invariant to such a rescaling, so the
empirical (sandwich) variance of the partial correlations,
`mc_partial_cor_delta(variability = "empirical")`, and the two-stage
variance, `variability = "two_stage"`, agree to first order, and so does
the null distribution of the Wald test that an edge parameter is zero.
The correction remains for the parameters that depend on the scale: the
precision coefficients of the edges that are present and the diagonal
terms. The model-based variance agrees with the two empirical forms only
if, in addition, the working fourth moments are correct, as under
normality. When the residual correlation changes with the covariates, as
it can with count responses, the agreement is only approximate and the
two-stage variance is the form to report.

## Structure selection

In this sample, stepwise selection with BIC returns the three edges of
the chain used in the simulation.

``` r

sel <- mc_graph_stepwise(fit, criterion = "bic")
sel$edges
#>      [,1] [,2]
#> [1,]    1    2
#> [2,]    2    3
#> [3,]    3    4
```

The joint Wald test of the hypothesis that all edge parameters are zero
rejects it here. The fit is also an `mcglm` object, so the `mcglm`
methods (`coef`, `vcov`, `gof`, …) apply to it.

``` r

mc_test_independence(fit)
#> 
#>  Wald test of between-response conditional association
#> 
#> data:  GGM layer with 6 edges
#> X-squared = 60.215, df = 6, p-value = 4.071e-11
coef(fit)        # inherited mcglm method
#>      Estimates Parameters        Type Response
#> 1  -0.05565639     beta10        beta        1
#> 2   0.54683586     beta11        beta        1
#> 3   0.01179182     beta20        beta        2
#> 4  -0.02413884     beta21        beta        2
#> 5  -0.02278576     beta30        beta        3
#> 6   0.05732429     beta31        beta        3
#> 7  -0.02369288     beta40        beta        4
#> 8  -0.23584547     beta41        beta        4
#> 9   0.26238020      rho12 correlation       NA
#> 10 -0.00334519      rho13 correlation       NA
#> 11 -0.01427623      rho14 correlation       NA
#> 12  0.31956896      rho23 correlation       NA
#> 13  0.05626320      rho24 correlation       NA
#> 14  0.23267150      rho34 correlation       NA
#> 15  0.59120324      tau11         tau        1
#> 16  0.58985460      tau21         tau        2
#> 17  0.65344530      tau31         tau        3
#> 18  0.58735959      tau41         tau        4
```

## The DAG layer

The same data can be modeled with the directed acyclic graph layer,
which uses the modified Cholesky decomposition and guarantees positive
definiteness. In `mcglmgraph`, this is an ordered residual
parametrization. A causal reading requires additional subject-matter
assumptions, such as a defensible ordering, causal sufficiency and
faithfulness. When no ordering is justified, the CPDAG layer
([`mc_graph_pipeline()`](https://evandeilton.github.io/mcglmgraph/reference/mc_graph_pipeline.md),
[`mc_ges()`](https://evandeilton.github.io/mcglmgraph/reference/mc_ges.md))
learns the equivalence class instead.

``` r

fit_dag <- mcglm_graph(list(y1 ~ x, y2 ~ x, y3 ~ x, y4 ~ x),
                       between_graph = mc_between_graph("dag"), data = dat)
#> Automatic initial values selected.
mc_compute_criterion(fit_dag, c("aic", "bic"))
#>      aic      bic 
#> 3344.134 3381.172 
#> attr(,"criterion_status")
#> [1] "not_required"
#> attr(,"effective_dimension")
#> [1] NA
```

## References

Bonat, W. H. and Jorgensen, B. (2016). Multivariate covariance
generalized linear models. *Journal of the Royal Statistical Society:
Series C*, 65(5), 649–675. <doi:10.1111/rssc.12145>
