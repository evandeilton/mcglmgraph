# mcglmgraph

## Overview

`mcglmgraph` adds a graph to the covariance between responses in
Multivariate Covariance Generalized Linear Models (MCGLM), for Gaussian
and non-Gaussian responses. The method, MCGLM-Graph, is fitted in two
stages:

1.  The first stage fits the marginal MCGLM with the
    [`mcglm`](https://CRAN.R-project.org/package=mcglm) package and
    extracts its Pearson residuals.
2.  The second stage imposes a graph on the covariance of these
    residuals and estimates the graph parameters with the Pearson
    estimating functions of the MCGLM, without a joint likelihood.

Graph edges therefore describe residual association among the responses
after the marginal means are removed. A missing edge is read as absence
of second-order conditional association, which coincides with
conditional independence under joint normality of the residuals.

Three layers are available:

| Layer | Structure | Parametrization |
|----|----|----|
| GGM | Gaussian graphical model (undirected) | Precision matrix as a linear predictor; a zero is a null partial correlation |
| DAG | Directed acyclic graph | Modified Cholesky (LDL) decomposition under a fixed ordering of the responses |
| CPDAG | Completed partially directed acyclic graph | Markov equivalence class learned by partial-correlation tests or greedy equivalence search |

The analytical derivatives of the graph layer are implemented in C++
through `RcppArmadillo`.

## Installation

Install the development version from GitHub:

``` r

# install.packages("remotes")
remotes::install_github("evandeilton/mcglmgraph")
```

The package compiles C++17 code, so a working toolchain is required
(Rtools on Windows, Xcode command line tools on macOS).

## Example

Four responses whose residual association follows the chain 1 - 2 - 3 -
4, with a covariate acting on two of them.

``` r

library(mcglmgraph)

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

# Complete GGM layer on the residuals of four Gaussian regressions
fit <- mcglm_graph(list(y1 ~ x, y2 ~ x, y3 ~ x, y4 ~ x),
                   between_graph = mc_between_graph("ggm"), data = dat)

# Partial correlations with delta-method standard errors
mc_partial_cor_delta(fit)

# Select the support by backward search under an information criterion
sel <- mc_graph_stepwise(fit, direction = "backward", criterion = "bic")
sel
mc_plot_network(sel$fit)
```

## Main functions

| Task | Functions |
|----|----|
| Specify and fit | [`mc_between_graph()`](https://evandeilton.github.io/mcglmgraph/reference/mc_between_graph.md), [`mcglm_graph()`](https://evandeilton.github.io/mcglmgraph/reference/mcglm_graph.md) |
| Read the graph layer | [`mc_partial_cor()`](https://evandeilton.github.io/mcglmgraph/reference/mc_partial_cor.md), [`mc_partial_cor_delta()`](https://evandeilton.github.io/mcglmgraph/reference/mc_partial_cor_delta.md), [`mc_dag_coefficients()`](https://evandeilton.github.io/mcglmgraph/reference/mc_dag_coefficients.md), [`mc_get_graph()`](https://evandeilton.github.io/mcglmgraph/reference/mc_get_graph.md), [`mc_plot_network()`](https://evandeilton.github.io/mcglmgraph/reference/mc_plot_network.md) |
| Inference | [`mc_vcov_graph()`](https://evandeilton.github.io/mcglmgraph/reference/mc_vcov_graph.md), [`mc_wald_graph()`](https://evandeilton.github.io/mcglmgraph/reference/mc_wald_graph.md), [`mc_test_independence()`](https://evandeilton.github.io/mcglmgraph/reference/mc_test_independence.md), [`mc_test_adequacy()`](https://evandeilton.github.io/mcglmgraph/reference/mc_test_adequacy.md) |
| Select an undirected graph | [`mc_graph_stepwise()`](https://evandeilton.github.io/mcglmgraph/reference/mc_graph_stepwise.md), [`mc_graph_penalized()`](https://evandeilton.github.io/mcglmgraph/reference/mc_graph_penalized.md), [`mc_select_lambda()`](https://evandeilton.github.io/mcglmgraph/reference/mc_select_lambda.md) |
| Discover a directed structure | [`mc_graph_pipeline()`](https://evandeilton.github.io/mcglmgraph/reference/mc_graph_pipeline.md), [`mc_ges()`](https://evandeilton.github.io/mcglmgraph/reference/mc_ges.md), [`mc_pc_skeleton()`](https://evandeilton.github.io/mcglmgraph/reference/mc_pc_skeleton.md) |
| Stability and resampling | [`mc_graph_stability()`](https://evandeilton.github.io/mcglmgraph/reference/mc_graph_stability.md), [`mc_cpdag_stability()`](https://evandeilton.github.io/mcglmgraph/reference/mc_cpdag_stability.md), [`mc_graph_bootstrap()`](https://evandeilton.github.io/mcglmgraph/reference/mc_graph_bootstrap.md) |
| Compare with a known structure | [`mc_structure_metrics()`](https://evandeilton.github.io/mcglmgraph/reference/mc_structure_metrics.md), [`mc_cpdag_metrics()`](https://evandeilton.github.io/mcglmgraph/reference/mc_cpdag_metrics.md), [`mc_sid()`](https://evandeilton.github.io/mcglmgraph/reference/mc_sid.md) |

The vignettes cover an introduction, a comparison with a plain MCGLM fit
and the examples of the `mcglm` package revisited with a graph layer:

``` r

browseVignettes("mcglmgraph")
```

## References

Bonat, W. H. and Jorgensen, B. (2016). Multivariate covariance
generalized linear models. *Journal of the Royal Statistical Society:
Series C*, 65(5), 649-675. <https://doi.org/10.1111/rssc.12145>

Bonat, W. H. (2018). Multiple response variables regression models in R:
the mcglm package. *Journal of Statistical Software*, 84(4), 1-30.
<https://doi.org/10.18637/jss.v084.i04>

Murphy, K. M. and Topel, R. H. (1985). Estimation and inference in
two-step econometric models. *Journal of Business & Economic
Statistics*, 3(4), 370-379.
<https://doi.org/10.1080/07350015.1985.10509471>

Pourahmadi, M. (1999). Joint mean-covariance models with applications to
longitudinal data: unconstrained parameterisation. *Biometrika*, 86(3),
677-690. <https://doi.org/10.1093/biomet/86.3.677>

## Citation

``` r

citation("mcglmgraph")
```

## License

GPL-3. See
[LICENSE.md](https://evandeilton.github.io/mcglmgraph/LICENSE.md).
