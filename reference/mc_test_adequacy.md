# Adequacy of a support against the saturated graph

Tests whether the edges a support omits can all be zero at once. The
saturated graph is refitted on the same residuals, complete in the
undirected layer and complete along the same ordering in the directed
one, and the parameters of the omitted edges are tested jointly by
[`mc_wald_graph`](https://evandeilton.github.io/mcglmgraph/reference/mc_wald_graph.md)
under the covariance of that saturated fit (an adequacy check of a
selected structure). With \\E\_{\max}\\ the set of all pairs of
responses and \\\hat E\\ the edges of the support, the null hypothesis
is \\H_0:\theta_e=0\\ for every \\e\in E\_{\max}\setminus\hat E\\,
against the saturated graph, and the statistic is the Wald statistic of
[`mc_wald_graph`](https://evandeilton.github.io/mcglmgraph/reference/mc_wald_graph.md)
with the contrast that selects those \\\|E\_{\max}\|-\|\hat E\|\\
parameters.

## Usage

``` r
mc_test_adequacy(
  fit,
  type = c("sandwich", "model", "two_stage"),
  cluster = NULL,
  small_sample = FALSE,
  method = c("auto", "shortcut", "full")
)
```

## Arguments

- fit:

  an object of class `"mcglm_graph"` whose support omits at least one
  edge; the saturated graph is refitted from its residuals and must
  converge.

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

An object of class `"htest"` as in
[`mc_wald_graph`](https://evandeilton.github.io/mcglmgraph/reference/mc_wald_graph.md),
whose degrees of freedom are the number of omitted edges, with the extra
components `omitted` (a two-column matrix `i`, `j` of the omitted pairs
of responses), `n_saturated` (the number of edges of the saturated
graph) and `descriptive_after_selection`, always `TRUE`.

## A diagnostic, not a test at nominal level

If the omitted edges were fixed before looking at the data, the
statistic would be referred to a chi-squared distribution with as many
degrees of freedom as edges omitted. When the support was chosen by a
search on the same data, that reference is only the one of a fixed
support: the result is descriptive, with no claim of nominal level after
selection. A large p-value says that the zeros of the support are
compatible with the residual covariance; a small one says that at least
one omitted edge carries partial association. The returned object
records this in `descriptive_after_selection`.

## See also

[`mc_wald_graph`](https://evandeilton.github.io/mcglmgraph/reference/mc_wald_graph.md),
[`mc_test_independence`](https://evandeilton.github.io/mcglmgraph/reference/mc_test_independence.md)
and
[`mc_graph_stepwise`](https://evandeilton.github.io/mcglmgraph/reference/mc_graph_stepwise.md).

Other inference on the graph layer:
[`mc_test_independence()`](https://evandeilton.github.io/mcglmgraph/reference/mc_test_independence.md),
[`mc_two_stage_available()`](https://evandeilton.github.io/mcglmgraph/reference/mc_two_stage_available.md),
[`mc_vcov_graph()`](https://evandeilton.github.io/mcglmgraph/reference/mc_vcov_graph.md),
[`mc_wald_graph()`](https://evandeilton.github.io/mcglmgraph/reference/mc_wald_graph.md)

## Examples

``` r
# \donttest{
set.seed(1); n <- 400
K <- diag(4) * 2; K[1, 2] <- K[2, 1] <- -0.7; K[3, 4] <- K[4, 3] <- -0.7
Y <- matrix(rnorm(n * 4), n) %*% chol(solve(K))
d <- data.frame(y1 = Y[, 1], y2 = Y[, 2], y3 = Y[, 3], y4 = Y[, 4],
                x = rnorm(n))
f <- mcglm_graph(list(y1 ~ x, y2 ~ x, y3 ~ x, y4 ~ x), data = d)
#> Automatic initial values selected.
sel <- mc_graph_stepwise(f, criterion = "bic", trace = FALSE)
mc_test_adequacy(sel$fit)
#> 
#>  Wald test of adequacy against the saturated graph (sandwich covariance;
#>  descriptive after selection)
#> 
#> data:  4 of 6 edges omitted from a GGM layer
#> X-squared = 3.3534, df = 4, p-value = 0.5005
#> sample estimates:
#> [1] -0.007881532  0.017153068  0.035587530 -0.115230131
#> 
# }
```
