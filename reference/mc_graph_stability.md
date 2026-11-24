# Selection-frequency aggregation of the between-response graph structure

Repeats graph selection in `B` half-samples and records the empirical
selection frequency \\\hat\Pi_e\\ of every edge, the proportion of valid
subsamples in which the selector retains it. Edges with
\\\hat\Pi_e\ge\pi\_{\mathrm{thr}}\\ are returned as recurrent. The
frequency estimates the probability that the selector retains the edge
in a random half-sample; it is not the probability that the edge belongs
to the graph. With a selector tuned by an information criterion, the
default, the procedure is a diagnostic aggregation scheme and not the
error-controlled stability selection of Meinshausen and Buehlmann
(2010), because the number of edges each subsample selects is not
controlled: the recurrent set then carries no formal bound on its false
selections. Passing `q` gives the version for which such a bound is
reported, under the conditions stated below.

## Usage

``` r
mc_graph_stability(
  fit,
  B = 100L,
  pi_thr = 0.6,
  method = c("stepwise", "penalized"),
  criterion = "bic",
  refit = c("residual", "full"),
  workers = 1L,
  cluster = NULL,
  q = NULL,
  sampling = c("half", "complementary"),
  ...
)

# S3 method for class 'mcg_stability'
print(x, ...)
```

## Arguments

- fit:

  an object of class `"mcglm_graph"` (GGM layer).

- B:

  number of subsamples (default 100); with `sampling = "complementary"`,
  the number of pairs, so \\2B\\ subsamples are drawn.

- pi_thr:

  selection-frequency threshold in \\\[0,1\]\\ (default 0.6).

- method:

  per-subsample selector: `"stepwise"`
  ([`mc_graph_stepwise`](https://evandeilton.github.io/mcglmgraph/reference/mc_graph_stepwise.md),
  default) or `"penalized"`
  ([`mc_select_lambda`](https://evandeilton.github.io/mcglmgraph/reference/mc_select_lambda.md)).
  It is set to `"penalized"` when `q` is given.

- criterion:

  information criterion passed to the selector (`"bic"`, the default,
  `"aic"`, `"qic"`, `"gbic"`, its alias `"clbic"`, or `"ebic"`); not
  used when `q` is given.

- refit:

  `"residual"` (default) resamples the stored Pearson residuals with the
  first stage held fixed; `"full"` resamples rows of the data and refits
  both stages in every half-sample.

- workers:

  number of parallel workers for the resampling loop. The subsamples are
  drawn in the parent before any dispatch, so a given seed yields the
  same subsamples at any worker count and the result does not depend on
  it. Sequential by default, deliberately: user code often parallelizes
  at an outer level, and nesting a second level of workers underneath
  oversubscribes the machine. Forking is used where available; on
  Windows the call falls back to sequential evaluation. Inside a
  parallel worker the argument is ignored, with a warning.

- cluster:

  optional vector with one identifier per observational row; half of the
  clusters, rounded down, are drawn without replacement in every
  half-sample (at least two with `sampling = "half"`, at least one per
  half with `sampling = "complementary"`).

- q:

  optional integer between 1 and \\R(R-1)/2\\: the number of edges each
  subsample selects along the L1 path, which turns on the
  error-controlled procedure described above and requires
  `method = "penalized"` when `method` is given.

- sampling:

  `"half"` (default), `B` independent half-samples drawn without
  replacement, or `"complementary"`, `B` pairs of disjoint halves.

- ...:

  further arguments passed to the per-subsample selector
  ([`mc_graph_stepwise`](https://evandeilton.github.io/mcglmgraph/reference/mc_graph_stepwise.md)
  or
  [`mc_select_lambda`](https://evandeilton.github.io/mcglmgraph/reference/mc_select_lambda.md));
  ignored when `q` is given. For the `print` method, ignored.

- x:

  an object of class `"mcg_stability"`.

## Value

An object of class `"mcg_stability"` with the empirical
selection-frequency matrix among valid subsamples (`selection_prob`),
its conservative counterpart divided by all attempted subsamples
(`selection_prob_conservative`), the numbers of valid and failed
subsamples (`B_valid`, `B_failed`), and the recurrent edge set
(`stable_edges`, a two-column matrix of the pairs whose frequency among
valid subsamples is at least `pi_thr`). With `q` it also carries `q` and
the `error_bound`, `NA` when `pi_thr` is not above one half. The list
also records `B`, `pi_thr`, `method`, `criterion`, `refit`, `sampling`,
`clustered`, `n_subsamples` and the data frame `failures` (subsample
index and error message of each failed subsample). The matrices have the
response names as dimnames and hold the frequency of each pair in both
triangles. `print` returns its argument invisibly.

## Stability selection with error control

With `q`, every subsample selects the first `q` edges to enter along the
L1 path of
[`mc_graph_penalized`](https://evandeilton.github.io/mcglmgraph/reference/mc_graph_penalized.md),
run with warm starts over a geometric grid of 50 penalties from the
value that zeroes every edge in that subsample, \\\lambda\_{\max}\\ (or
\\10^{-3}\\ if larger), down to \\10^{-3}\\ times it; edges that enter
at the same grid point are ranked by the size of their coefficient
there. Grid points whose fit is not certified are skipped, and fewer
than `q` edges are selected when fewer enter along the grid, so \\q\\ is
an upper limit on the number of edges per subsample. The expected number
of falsely selected edges among those whose frequency reaches `pi_thr`
is then bounded by Theorem 1 of Meinshausen and Buehlmann (2010),
\$\$E(V)\le\frac{q^2}{(2\pi\_{\mathrm{thr}}-1)\\\|E\_{\max}\|},\qquad
\tfrac12\<\pi\_{\mathrm{thr}}\le1,\$\$ where \\V\\ is the number of
falsely selected edges, \\q\\ the number of edges each subsample selects
and \\\|E\_{\max}\|=R(R-1)/2\\ the number of pairs, reported as
`error_bound`. The bound rests on two conditions the package cannot
check: the selection indicators of the null edges are exchangeable, and
the selector is no worse than random guessing. It is a bound on an
expected count, not on a probability, and it refers to frequencies
computed over all the subsamples; when some subsamples fail,
`selection_prob_conservative` is the frequency that counts a failed
subsample as selecting no edge, whereas the recurrent set is read from
the frequencies among the valid ones. With `sampling = "complementary"`
the subsamples come in `B` pairs of disjoint halves, the complementary
pairs of Shah and Samworth (2013). What the bound controls is the edges
of the residual layer as the second stage sees them; with
`refit = "residual"` the first stage is held fixed.

## What each half-sample is

With `refit = "residual"` the half-sample is drawn from the stored
Pearson residuals and the first stage is held fixed, so the frequencies
describe the second stage alone. With `refit = "full"` the half-sample
is drawn from the rows of the data, the first-stage MCGLM (multivariate
covariance generalized linear model) is refitted on it, and the selector
runs on the residuals of that refit, so the frequencies carry the
variability of both stages. The second form accounts for the estimation
of the first stage and costs an MCGLM fit per half-sample; the first is
the default and is cheaper.

A half-sample holds \\\lfloor n/2\rfloor\\ of the \\n\\ rows, drawn
without replacement, and at least two. With `sampling = "complementary"`
one random permutation of the rows is split into two disjoint halves of
\\\lfloor n/2\rfloor\\ rows each (at least one), so that with \\n\\ odd
one row is left out of the pair.

With `cluster`, half of the clusters are drawn instead of half of the
rows, under the same rounding and the same minimum sizes counted in
clusters, and every row of a drawn cluster enters the half-sample. That
is the form to use when the unit of dependence is a cluster of rows, and
it is required with `refit = "full"` when the within-response structure
is not the identity: each drawn cluster then keeps its own block of the
structure, as in
[`mc_graph_bootstrap`](https://evandeilton.github.io/mcglmgraph/reference/mc_graph_bootstrap.md).

## References

Meinshausen, N. and Buehlmann, P. (2010). Stability selection. *Journal
of the Royal Statistical Society: Series B* 72(4), 417–473.
[doi:10.1111/j.1467-9868.2010.00740.x](https://doi.org/10.1111/j.1467-9868.2010.00740.x)

Shah, R. D. and Samworth, R. J. (2013). Variable selection with error
control: another look at stability selection. *Journal of the Royal
Statistical Society: Series B* 75(1), 55–80.
[doi:10.1111/j.1467-9868.2011.01034.x](https://doi.org/10.1111/j.1467-9868.2011.01034.x)

## See also

[`mc_graph_stepwise`](https://evandeilton.github.io/mcglmgraph/reference/mc_graph_stepwise.md),
[`mc_select_lambda`](https://evandeilton.github.io/mcglmgraph/reference/mc_select_lambda.md),
[`mc_graph_penalized`](https://evandeilton.github.io/mcglmgraph/reference/mc_graph_penalized.md).

Other structure selection:
[`mc_edge_scores()`](https://evandeilton.github.io/mcglmgraph/reference/mc_edge_scores.md),
[`mc_graph_penalized()`](https://evandeilton.github.io/mcglmgraph/reference/mc_graph_penalized.md),
[`mc_graph_stepwise()`](https://evandeilton.github.io/mcglmgraph/reference/mc_graph_stepwise.md),
[`mc_select_lambda()`](https://evandeilton.github.io/mcglmgraph/reference/mc_select_lambda.md)

## Examples

``` r
# \donttest{
set.seed(1); n <- 300
K <- diag(4)*2; K[1,2]<-K[2,1]<- -0.6; K[2,3]<-K[3,2]<- -0.6; K[3,4]<-K[4,3]<- -0.6
Y <- matrix(rnorm(n * 4), n) %*% chol(solve(K)); x <- rnorm(n)
d <- data.frame(y1=Y[,1], y2=Y[,2], y3=Y[,3], y4=Y[,4], x=x)
f <- mcglm_graph(list(y1~x,y2~x,y3~x,y4~x), data = d)
#> Automatic initial values selected.
mc_graph_stability(f, B = 40)
#> Selection-frequency aggregation (stepwise, criterion = bic, B = 40, threshold = 0.60)
#> Half-samples: residuals resampled, first stage held fixed
#> Valid subsamples: 40 | failed: 0
#> Recurrent edges: 3
#> Edges:
#>      [,1] [,2]
#> [1,]    1    2
#> [2,]    2    3
#> [3,]    3    4
#> 
#> Selection probabilities:
#>      y1   y2   y3   y4
#> y1 0.00 1.00 0.00 0.03
#> y2 1.00 0.00 1.00 0.05
#> y3 0.00 1.00 0.00 0.98
#> y4 0.03 0.05 0.98 0.00
# }
```
