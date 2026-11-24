# Gaussian working log-criterion of the between-response fit

Evaluates the Gaussian working log-criterion of the Pearson residuals
under the estimated between-response covariance, \$\$\ell^p_e =
-\tfrac{n}{2}\bigl\[R\log(2\pi) + \log\|\hat\Sigma_b\| +
\mathrm{tr}(\hat\Sigma_b^{-1}\hat S_e)\bigr\],\$\$ where \\n\\ is the
number of units, \\R\\ the number of responses, \\\hat\Sigma_b\\ the
fitted between-response covariance and \\\hat S_e\\ the sample
covariance (centered, divided by \\n\\) of the Pearson residuals of the
first stage. This scalar criterion generates the Pearson equations used
by the graph layer, and it is the quantity from which
[`mc_compute_criterion`](https://evandeilton.github.io/mcglmgraph/reference/mc_compute_criterion.md)
builds the information criteria. It is not presented as a composite
likelihood or as a fully specified likelihood for the marginal MCGLM
(multivariate covariance generalized linear model) responses.

## Usage

``` r
mc_graph_ploglik(fit)
```

## Arguments

- fit:

  an object of class `"mcglm_graph"`.

## Value

A single number, the Gaussian working log-criterion.

## See also

[`mc_compute_criterion`](https://evandeilton.github.io/mcglmgraph/reference/mc_compute_criterion.md),
[`mc_graph_cv`](https://evandeilton.github.io/mcglmgraph/reference/mc_graph_cv.md).

Other information criteria:
[`mc_compute_criterion()`](https://evandeilton.github.io/mcglmgraph/reference/mc_compute_criterion.md),
[`mc_graph_cv()`](https://evandeilton.github.io/mcglmgraph/reference/mc_graph_cv.md)

## Examples

``` r
# \donttest{
set.seed(1); n <- 150; x <- rnorm(n)
d <- data.frame(y1 = rnorm(n), y2 = rnorm(n), y3 = rnorm(n), x = x)
f <- mcglm_graph(list(y1~x, y2~x, y3~x), data = d)
#> Automatic initial values selected.
mc_graph_ploglik(f)
#> [1] -635.3592
# }
```
