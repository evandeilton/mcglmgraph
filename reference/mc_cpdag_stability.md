# Resampling frequencies of the skeleton and of each orientation

Repeats the constraint-based discovery of
[`mc_graph_pipeline`](https://evandeilton.github.io/mcglmgraph/reference/mc_graph_pipeline.md)
over resamples and records, for every pair of responses, how often the
pair came out absent, undirected, or directed each way. An orientation
that the data barely support shows up here as a pair whose two
directions split the frequency, which a single run cannot reveal.

## Usage

``` r
mc_cpdag_stability(
  fit,
  B = 100L,
  alpha = 0.05,
  orient_method = c("pc", "hc", "tabu"),
  orient_rule = c("standard", "conservative", "majority"),
  required = NULL,
  forbidden = NULL,
  refit = c("residual", "full"),
  cluster = NULL,
  max_ord = NULL,
  df_residual = NULL,
  workers = 1L
)

# S3 method for class 'mcg_cpdag_stability'
print(x, ...)
```

## Arguments

- fit:

  an object of class `"mcglm_graph"`.

- B:

  number of resamples, an integer of at least two; default 100.

- alpha:

  significance level of the second-order separation test.

- orient_method:

  `"pc"` (native, default), `"hc"` or `"tabu"` (bnlearn). When bnlearn
  is not installed, the last two fall back on `"pc"` with a warning, and
  the result records it in `orientation_fallback`.

- orient_rule:

  rule for the unshielded triples of the native route: `"standard"`
  (default), `"conservative"` or `"majority"`. It has no effect when
  `orient_method` is a score-based route.

- required, forbidden:

  optional background knowledge, as two-column matrices of arcs \\(from,
  to)\\ given either by response name or by index. An arc in `required`
  must be present with that direction, an arc in `forbidden` must not.
  They are imposed before propagation, and Meek's fourth rule, idle
  without them, carries their consequences through the rest of the
  graph. Available on the native route only.

- refit:

  `"residual"` (default) or `"full"`, as above; `"full"` corresponds to
  `stage = "both"` of
  [`mc_graph_bootstrap`](https://evandeilton.github.io/mcglmgraph/reference/mc_graph_bootstrap.md).

- cluster:

  optional vector with one identifier per row, with no missing values
  and at least two distinct values; whole clusters are resampled.
  Required with `refit = "full"` when the within-response structure is
  not the identity, as in
  [`mc_graph_bootstrap`](https://evandeilton.github.io/mcglmgraph/reference/mc_graph_bootstrap.md).

- max_ord:

  maximum conditioning order; `NULL` uses \\R - 2\\.

- df_residual:

  optional residual degrees of freedom \\\nu_e\\: `NULL` (default) uses
  the raw-sample contract based on \\n\\; a number or the result of
  [`mc_residual_df`](https://evandeilton.github.io/mcglmgraph/reference/mc_residual_df.md)
  is used as given; `"auto"` computes it with
  [`mc_residual_df`](https://evandeilton.github.io/mcglmgraph/reference/mc_residual_df.md)
  and refuses a fit for which it is not exact.

- workers:

  number of parallel workers, default 1. The resamples are drawn in the
  parent, so a given seed yields the same ones at any worker count.

- x:

  an object of class `"mcg_cpdag_stability"`.

- ...:

  ignored.

## Value

An object of class `"mcg_cpdag_stability"`: a list with `states`, a data
frame with one row per pair of responses (`i`, `j`, `response1`,
`response2`) and the frequencies of the four states (`absent`,
`undirected`, `i_to_j`, `j_to_i`); the symmetric matrix `edge_prob` of
skeleton frequencies, one minus `absent`; `consensus`, the 0/1 adjacency
matrix that gives each pair its most frequent state, the first in the
order above on a tie (an entry `[i, j]` equal to one alone reads \\i\to
j\\, with `[j, i]` also one the edge is undirected); being formed pair
by pair, it need not be the CPDAG of a DAG; `valid_prob`, the proportion
of valid resamples whose output passed the validity gates of
[`mc_graph_pipeline`](https://evandeilton.github.io/mcglmgraph/reference/mc_graph_pipeline.md);
`mean_ambiguous`, the average number of ambiguous triples per valid
resample; the counts `B`, `B_valid` and `B_failed`; the data frame
`failures` (`replicate`, `error`); and the settings `refit`, `alpha`,
`orient_method`, `orient_rule`, `clustered` (logical) and
`response_names`. The `print` method shows the settings, the counts and
the table of frequencies, and returns its argument invisibly.

## What each resample is

With `refit = "residual"` the resample is drawn from the stored Pearson
residuals and the first stage is held fixed, so the frequencies describe
the discovery step alone. With `refit = "full"` rows of the data are
resampled and both stages are refitted, so they carry the variability of
the whole procedure. Units are drawn with replacement, keeping the
sample size, which is the bootstrap convention (Efron, 1979);
[`mc_graph_stability`](https://evandeilton.github.io/mcglmgraph/reference/mc_graph_stability.md)
instead uses half-samples, the convention of the stability-selection
literature (Meinshausen and Buhlmann, 2010). With \\B_v\\ valid
resamples, the frequency of a state is the proportion of valid resamples
in which the pair is in it, so the four frequencies of a pair sum to
one; a resample whose refit or pipeline raises an error is dropped and
recorded in `failures`, and at least one valid resample is required.

## What the frequencies estimate

Each frequency estimates, by resampling, the probability that the
discovery procedure, at the given `alpha` and orientation rule, returns
that state for the pair in a sample of the same size. It measures the
stability of the output of the procedure. It is not the probability that
the edge or the orientation belongs to the structure that generated the
data, and it carries no control of an error rate.

## References

Efron, B. (1979). Bootstrap methods: another look at the jackknife. *The
Annals of Statistics* 7(1), 1–26.
[doi:10.1214/aos/1176344552](https://doi.org/10.1214/aos/1176344552)

Meinshausen, N. and Buhlmann, P. (2010). Stability selection. *Journal
of the Royal Statistical Society, Series B* 72(4), 417–473.
[doi:10.1111/j.1467-9868.2010.00740.x](https://doi.org/10.1111/j.1467-9868.2010.00740.x)

## See also

[`mc_graph_pipeline`](https://evandeilton.github.io/mcglmgraph/reference/mc_graph_pipeline.md)
for a single run and
[`mc_graph_bootstrap`](https://evandeilton.github.io/mcglmgraph/reference/mc_graph_bootstrap.md)
for the bootstrap of the graph parameters.

Other resampling:
[`mc_graph_bootstrap()`](https://evandeilton.github.io/mcglmgraph/reference/mc_graph_bootstrap.md)

## Examples

``` r
# \donttest{
set.seed(1); n <- 400; x <- rnorm(n)
z1 <- rnorm(n); z2 <- rnorm(n); z3 <- 0.7 * z1 + 0.7 * z2 + rnorm(n)
d <- data.frame(y1 = z1, y2 = z2, y3 = z3, x = x)
f <- mcglm_graph(list(y1 ~ x, y2 ~ x, y3 ~ x), data = d)
#> Automatic initial values selected.
mc_cpdag_stability(f, B = 50)
#> Structure frequencies (pc, standard rule, alpha = 0.05, B = 50)
#> Resamples: residuals resampled, first stage held fixed | valid: 50 | failed: 0
#> Outputs passing the validity gates: 100%
#> 
#>   response1 response2 absent undirected i_to_j j_to_i
#> 1        y1        y2   0.94       0.06   0.00      0
#> 2        y1        y3   0.00       0.06   0.94      0
#> 3        y2        y3   0.00       0.06   0.94      0
# }
```
