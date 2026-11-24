# Edge score and Wald statistics

Computes, for every pair of responses, a chi-squared(1) statistic that
ranks edges of the Gaussian graphical model (GGM) layer without
refitting the whole model. The precision matrix is parametrized
linearly, \\K_b = \sum_d \eta_d Z_d\\, with one diagonal parameter per
response and one parameter \\\eta_e\\ per edge, so that an edge is
absent exactly when \\\eta_e = 0\\. Present edges are scored by the Wald
statistic \\\hat\eta_e^2 / \widehat{\mathrm{Var}}(\hat\eta_e)\\, with
the variance read from the model-based covariance \\J^{-1}\\ stored in
the fit. Absent edges are scored by the *efficient* Pearson score, that
is, the score of the candidate coordinate after projecting out the
nuisance parameters \\\lambda\\ of the current model, \$\$u_e = \psi_e -
J\_{e\lambda} J\_{\lambda\lambda}^{-1}\psi\_\lambda, \qquad S_e = u_e^2
/ \widehat{\mathrm{Var}}(u_e),\$\$ where \\\psi\\ is the second-order
Pearson estimating function and \\J\\ its sensitivity. With
`variability = "model"` the variance is taken under the working identity
\\V = J\\ (second Bartlett identity), in which case it equals the Schur
complement \\J\_{e\cdot\lambda} = J\_{ee} -
J\_{e\lambda}J\_{\lambda\lambda}^{-1} J\_{\lambda e}\\; with
`variability = "empirical"` it uses the empirical fourth-moment
variability of the Pearson contributions, which is the appropriate
choice when the second Bartlett identity fails.

## Usage

``` r
mc_edge_scores(
  fit,
  which = c("both", "wald", "score"),
  variability = c("model", "empirical")
)
```

## Arguments

- fit:

  an object of class `"mcglm_graph"` (GGM layer).

- which:

  one of `"both"` (default), `"wald"` (present edges only) or `"score"`
  (absent edges only).

- variability:

  `"model"` (default, working \\V = J\\) or `"empirical"` (fourth-moment
  variability) for the absent-edge scores. It does not affect the Wald
  statistics of present edges.

## Value

A data frame with one row per pair of responses (restricted by `which`)
and columns `i`, `j` (response indices with \\i\<j\\), `present`
(logical, edge in the support of the fit), `stat` (Wald statistic for
present edges, efficient score for absent ones), `p_value`
(chi-squared(1) upper tail of `stat`), `stat_unadjusted` (unadjusted
score \\\psi_e^2/J\_{ee}\\, only for absent edges, `NA` otherwise) and
`variance` (estimated variance of \\\hat\eta_e\\ for present edges, of
\\u_e\\ for absent ones). The attribute `"variability"` records the
choice used. An edge whose estimated variance is not finite and positive
has `stat` and `p_value` equal to `NA`, and
[`mc_graph_stepwise`](https://evandeilton.github.io/mcglmgraph/reference/mc_graph_stepwise.md)
does not take a present edge in that state as a removal candidate.

## Details

Dividing \\\psi_e^2\\ by the unadjusted \\J\_{ee}\\ ignores the
estimation of the remaining parameters and, because
\\J\_{e\cdot\lambda}\le J\_{ee}\\, understates the statistic. That
quantity is returned as `stat_unadjusted`, for ranking, but it does not
carry a chi-squared calibration.

The reported `p_value` is a working reference distribution: it is the
chi-squared(1) tail of `stat`. Its validity requires the regularity
conditions of the second-order estimating functions and, for
`variability = "model"`, the working identity \\V = J\\. Applied to
edges chosen by a search over graphs, it remains descriptive.

The function refuses fits whose numerical diagnostics explicitly mark
graph-layer inference as uncertified. Responses are indexed \\i\<j\\
over all \\R(R-1)/2\\ pairs.

## See also

[`mc_graph_stepwise`](https://evandeilton.github.io/mcglmgraph/reference/mc_graph_stepwise.md),
which uses these statistics to rank candidate edges.

Other structure selection:
[`mc_graph_penalized()`](https://evandeilton.github.io/mcglmgraph/reference/mc_graph_penalized.md),
[`mc_graph_stability()`](https://evandeilton.github.io/mcglmgraph/reference/mc_graph_stability.md),
[`mc_graph_stepwise()`](https://evandeilton.github.io/mcglmgraph/reference/mc_graph_stepwise.md),
[`mc_select_lambda()`](https://evandeilton.github.io/mcglmgraph/reference/mc_select_lambda.md)

## Examples

``` r
# \donttest{
set.seed(1); n <- 200; x <- rnorm(n)
d <- data.frame(y1 = rnorm(n), y2 = rnorm(n), y3 = rnorm(n), x = x)
f <- mcglm_graph(list(y1~x, y2~x, y3~x), data = d)
#> Automatic initial values selected.
mc_edge_scores(f)
#>   i j present        stat   p_value stat_unadjusted    variance
#> 1 1 2    TRUE 0.125702663 0.7229300              NA 0.005119502
#> 2 1 3    TRUE 0.006452406 0.9359773              NA 0.005113397
#> 3 2 3    TRUE 0.323259426 0.5696554              NA 0.005129649
# }
```
