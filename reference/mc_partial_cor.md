# Partial correlations between responses

Computes the between-response partial correlation matrix
\\\rho\_{ij\mid\mathrm{rest}} = -K\_{ij}/\sqrt{K\_{ii} K\_{jj}}\\ from
the estimated precision matrix \\K_b = \Sigma_b^{-1}\\ of the graph
layer (Lauritzen, 1996). It is invariant to the scale of the Pearson
residuals, and it is computed from `fit$graph$K_b` for both the GGM and
the DAG layer. A zero entry signals the absence of second-order
conditional association between two responses given the others, which
coincides with conditional independence only under joint normality. In
Gaussian-copula models, that interpretation applies only when the
precision is estimated on the latent Gaussian scale.

## Usage

``` r
mc_partial_cor(fit)
```

## Arguments

- fit:

  an object of class `"mcglm_graph"`.

## Value

An \\R \times R\\ matrix of partial correlations with unit diagonal,
with the response names as dimnames.

## References

Lauritzen, S. L. (1996). *Graphical Models*. Oxford University Press.

## See also

[`mc_partial_cor_delta`](https://evandeilton.github.io/mcglmgraph/reference/mc_partial_cor_delta.md)
for standard errors and intervals,
[`mc_get_graph`](https://evandeilton.github.io/mcglmgraph/reference/mc_get_graph.md),
[`mc_plot_network`](https://evandeilton.github.io/mcglmgraph/reference/mc_plot_network.md).

Other graph interpretation:
[`mc_get_graph()`](https://evandeilton.github.io/mcglmgraph/reference/mc_get_graph.md),
[`mc_partial_cor_delta()`](https://evandeilton.github.io/mcglmgraph/reference/mc_partial_cor_delta.md),
[`mc_plot_network()`](https://evandeilton.github.io/mcglmgraph/reference/mc_plot_network.md)

## Examples

``` r
set.seed(1); n <- 200; x <- rnorm(n)
d <- data.frame(y1 = rnorm(n), y2 = rnorm(n), y3 = rnorm(n), x = x)
f <- mcglm_graph(list(y1~x, y2~x, y3~x), data = d)
#> Automatic initial values selected.
mc_partial_cor(f)
#>              y1          y2           y3
#> y1  1.000000000 -0.02507805 -0.005680059
#> y2 -0.025078050  1.00000000 -0.040235728
#> y3 -0.005680059 -0.04023573  1.000000000
```
