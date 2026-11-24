# Bootstrap the graph layer with both stages refitted

Resamples units, or clusters of units, with replacement and refits the
first-stage MCGLM and the graph layer in every resample, returning the
bootstrap distribution of a statistic of the fit (by default the graph
parameters and the off-diagonal partial correlations). The support of
the graph is held fixed at the one the supplied fit carries, so the
result is an inferential check on a given structure and not a
post-selection estimator.

## Usage

``` r
mc_graph_bootstrap(
  fit,
  B = 200L,
  cluster = NULL,
  statistic = NULL,
  stage = c("both", "residual"),
  level = 0.95,
  workers = 1L
)

# S3 method for class 'mcg_bootstrap'
print(x, ...)
```

## Arguments

- fit:

  an object of class `"mcglm_graph"`.

- B:

  number of resamples, an integer of at least two; default 200.

- cluster:

  optional vector with one identifier per row, with no missing values
  and at least two distinct values; whole clusters are resampled, with
  replacement, keeping the total number of clusters. Required when the
  within-response structure is not the identity and `stage = "both"`.

- statistic:

  optional function of a refitted `"mcglm_graph"` object returning a
  named numeric vector, of the same length at every resample. The
  default returns the graph parameters and the off-diagonal partial
  correlations.

- stage:

  `"both"` refits the first stage and the graph layer; `"residual"`
  resamples the stored Pearson residuals and refits the graph layer
  alone.

- level:

  confidence level of the percentile intervals, a number strictly
  between 0 and 1; default 0.95.

- workers:

  number of parallel workers, default 1 (sequential). The resamples are
  drawn in the parent before any dispatch, so a given seed yields the
  same resamples at any worker count.

- x:

  an object of class `"mcg_bootstrap"`.

- ...:

  ignored.

## Value

An object of class `"mcg_bootstrap"`, a list with `estimate`, the named
numeric vector of the statistic on the original sample; `replicates`,
the `B_valid` by length-of-statistic matrix of valid resamples; `vcov`,
their covariance matrix, and `se`, their standard deviations; `bias`,
the mean of the replicates minus `estimate`; `ci`, the two-column matrix
(`lower`, `upper`) of percentile intervals; `level` and `stage`, as
given; the counts `B`, `B_valid` and `B_failed`; the data frame
`failures` (`replicate`, `error`); `clustered` (logical) and
`n_clusters` (`NA` without clusters). The `print` method shows the
estimates, standard errors, biases and intervals and returns its
argument invisibly.

## Details

With \\n\\ units (rows, or clusters when `cluster` is given), each
resample draws \\n\\ units with replacement (Efron, 1979), refits, and
evaluates the statistic \\T(\cdot)\\. With \\B_v\\ valid resamples
\\T^{\*}\_1,\dots,T^{\*}\_{B_v}\\, the bootstrap covariance is the
sample covariance of the replicates, `se` their standard deviations, the
bias estimate \\\bar T^{\*}-T(\hat{\mathrm F}\_n)\\, and the interval
for level \\1-a\\ the percentile interval with limits at the quantiles
\\a/2\\ and \\1-a/2\\ of the replicates. A resample is dropped, and
recorded with its message in `failures`, when its refit raises an error,
when its graph layer does not converge, or when its statistic is not a
finite vector of the length obtained on the original sample; the
summaries use the valid resamples only, and at least two are required.
Resampling independent clusters is the scheme that remains valid when
rows within a cluster are dependent. Because both stages are refitted,
the replicates carry the variability of the first stage, which the
model-based and sandwich covariances of
[`mc_vcov_graph`](https://evandeilton.github.io/mcglmgraph/reference/mc_vcov_graph.md)
take as known.

## What it is for

The variances the package reports elsewhere are conditional on the
estimated Pearson residuals unless the analytical correction of
[`mc_vcov_graph`](https://evandeilton.github.io/mcglmgraph/reference/mc_vcov_graph.md)
is asked for. Refitting propagates the first stage by construction, for
any family, which makes this the yardstick against which both routes of
that correction, the shortcut and the full form, are checked. With
`stage = "residual"` only the residuals are resampled and the first
stage is held fixed, which is the comparison that isolates how much the
first stage contributes.

## What it requires

The first-stage call is rebuilt from the fitted object. With an identity
within-response structure the rows are resampled, or whole clusters when
`cluster` is given. A richer matrix linear predictor describes
dependence among rows, so the rows it connects must be drawn together:
`cluster` is then required, every matrix of the structure must vanish
between distinct clusters, which is checked, and each drawn cluster
keeps its own block of the structure. A cluster drawn twice gives two
blocks that do not touch, as two independent units would; the submatrix
of the rows drawn would instead join the two copies. Without `cluster`
such a fit is refused.

## References

Efron, B. (1979). Bootstrap methods: another look at the jackknife. *The
Annals of Statistics* 7(1), 1–26.
[doi:10.1214/aos/1176344552](https://doi.org/10.1214/aos/1176344552)

## See also

[`mc_vcov_graph`](https://evandeilton.github.io/mcglmgraph/reference/mc_vcov_graph.md)
for the analytical forms and
[`mc_graph_stability`](https://evandeilton.github.io/mcglmgraph/reference/mc_graph_stability.md)
for selection frequencies.

Other resampling:
[`mc_cpdag_stability()`](https://evandeilton.github.io/mcglmgraph/reference/mc_cpdag_stability.md)

## Examples

``` r
# \donttest{
set.seed(1); n <- 200; x <- rnorm(n)
d <- data.frame(y1 = rnorm(n), y2 = rnorm(n), y3 = rnorm(n), x = x)
f <- mcglm_graph(list(y1 ~ x, y2 ~ x, y3 ~ x), data = d)
#> Automatic initial values selected.
bt <- mc_graph_bootstrap(f, B = 40)
bt$se
#>      eta_11      eta_22      eta_33      eta_12      eta_13      eta_23 
#> 0.009336896 0.012431239 0.010030779 0.082969202 0.049911598 0.069535071 
#>  pcor_y1_y2  pcor_y1_y3  pcor_y2_y3 
#> 0.080411716 0.048663508 0.067242523 
# }
```
