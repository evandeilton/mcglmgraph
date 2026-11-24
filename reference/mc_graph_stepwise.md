# Stepwise selection of the dependence graph

Greedy search over supports of the GGM layer guided by an information
criterion (see
[`mc_compute_criterion`](https://evandeilton.github.io/mcglmgraph/reference/mc_compute_criterion.md)).
The backward search starts from the complete graph and, at each step,
removes the edge whose removal most reduces the criterion; the forward
search starts from the empty graph and adds the edge whose inclusion
most reduces it; `"both"` starts from the complete graph and considers
removals and additions at every step. Candidates are ranked by
[`mc_edge_scores`](https://evandeilton.github.io/mcglmgraph/reference/mc_edge_scores.md)
(smallest Wald statistic for removal, largest efficient score for
addition), and only the `k_confirm` top-ranked ones are refitted per
step, which is the score pruning that keeps the search affordable. The
search stops when no candidate lowers the criterion, or after
`max_steps` steps. Because the mean structure is held fixed, only the
between-response covariance is refitted at each step, which is justified
by the asymptotic orthogonality of the first-stage and graph-layer
estimating functions.

## Usage

``` r
mc_graph_stepwise(
  fit,
  direction = c("backward", "forward", "both"),
  criterion = c("bic", "aic", "qic", "gbic", "clbic", "ebic"),
  max_steps = NULL,
  k_confirm = 3L,
  trace = TRUE,
  ebic_gamma = 0.5,
  cluster = NULL,
  small_sample = FALSE
)
```

## Arguments

- fit:

  an object of class `"mcglm_graph"` (GGM layer).

- direction:

  `"backward"` (default), `"forward"` or `"both"`.

- criterion:

  `"bic"` (default), `"aic"`, `"qic"`, `"gbic"` or its alias `"clbic"`,
  or `"ebic"`.

- max_steps:

  maximum number of steps; `NULL` (default) means the number of possible
  edges, \\R(R-1)/2\\.

- k_confirm:

  number of top-ranked candidates refitted per step (default 3).

- trace:

  logical; record the accepted moves in the component `trace`.

- ebic_gamma:

  the parameter of the extended BIC, used only with
  `criterion = "ebic"`; see
  [`mc_compute_criterion`](https://evandeilton.github.io/mcglmgraph/reference/mc_compute_criterion.md).

- cluster:

  optional independent-cluster identifier, one per residual row, passed
  to
  [`mc_compute_criterion`](https://evandeilton.github.io/mcglmgraph/reference/mc_compute_criterion.md):
  Pearson contributions are summed within clusters in QIC and GBIC, and
  the number of clusters replaces the number of rows in the logarithmic
  penalty of BIC, GBIC and EBIC.

- small_sample:

  logical; apply the cluster factor \\G/(G-1)\\, with \\G\\ the number
  of clusters, to the empirical variability of the clustered criteria
  (default `FALSE`); ignored without `cluster`.

## Value

An object of class `"mcg_selection"`, a list with `method`
(`"stepwise"`), `criterion`, `graph` (the adjacency matrix of the
selected graph), `edges` (two-column matrix of selected edges),
`n_edges`, `trace` (data frame with `step`, `action`, `edge_i`, `edge_j`
and `criterion` for each accepted move, or `NULL`), `fit` (the
`"mcglm_graph"` refitted on the selected support) and
`candidate_failures` (data frame of the candidates that could not be
evaluated). The components `lambda`, `lambda_path`, `cpdag` and
`vstructures` are `NA` or `NULL` here.

## Details

The search is defined on the undirected layer: every candidate is a GGM
support, and the edge scores that rank the candidates are those of the
precision. A fit carrying the directed layer is refused rather than
converted, since the conversion would silently replace the structure the
caller fitted by a different one.

Candidates whose refit fails to converge or yields a non-finite
criterion are skipped and recorded in the component
`candidate_failures`.

## See also

[`mc_select_lambda`](https://evandeilton.github.io/mcglmgraph/reference/mc_select_lambda.md)
for the L1 route,
[`mc_graph_stability`](https://evandeilton.github.io/mcglmgraph/reference/mc_graph_stability.md)
for selection frequencies,
[`mcg_selection-methods`](https://evandeilton.github.io/mcglmgraph/reference/mcg_selection-methods.md).

Other structure selection:
[`mc_edge_scores()`](https://evandeilton.github.io/mcglmgraph/reference/mc_edge_scores.md),
[`mc_graph_penalized()`](https://evandeilton.github.io/mcglmgraph/reference/mc_graph_penalized.md),
[`mc_graph_stability()`](https://evandeilton.github.io/mcglmgraph/reference/mc_graph_stability.md),
[`mc_select_lambda()`](https://evandeilton.github.io/mcglmgraph/reference/mc_select_lambda.md)

## Examples

``` r
# \donttest{
set.seed(1); n <- 300
K <- diag(4) * 2; K[1,2]<-K[2,1]<- -0.6; K[2,3]<-K[3,2]<- -0.6; K[3,4]<-K[4,3]<- -0.6
Y <- matrix(rnorm(n * 4), n) %*% chol(solve(K)); x <- rnorm(n)
d <- data.frame(y1=Y[,1], y2=Y[,2], y3=Y[,3], y4=Y[,4], x=x)
f <- mcglm_graph(list(y1~x,y2~x,y3~x,y4~x), data = d)
#> Automatic initial values selected.
sel <- mc_graph_stepwise(f, criterion = "bic")
sel$n_edges
#> [1] 3
# }
```
