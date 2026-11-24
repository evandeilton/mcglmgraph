# Out-of-sample comparison of candidate supports

Splits the units into folds, fits the graph layer on the units outside
each fold and scores it on the units inside it, which compares candidate
supports without a penalty term and without scoring a support on the
sample that chose it. The assignment of units to folds is random, so
call [`set.seed()`](https://rdrr.io/r/base/Random.html) beforehand for a
reproducible split.

## Usage

``` r
mc_graph_cv(fit, supports = NULL, folds = 5L, workers = 1L)

# S3 method for class 'mcg_cv'
print(x, ...)
```

## Arguments

- fit:

  an object of class `"mcglm_graph"`.

- supports:

  a named list of candidate supports, each a two-column matrix of edges
  (rows \\(i,j)\\, response indices) or `NULL` for the complete graph.
  The default, `NULL`, compares the complete graph (named
  `"unstructured"`) with the support the fit carries (named `"fitted"`),
  when that one is sparse.

- folds:

  integer number of folds, at least two and at most the number of units
  (default 5).

- workers:

  number of parallel workers (default 1, sequential). The folds are
  formed in the parent before any dispatch, so a given seed yields the
  same split at any worker count.

- x:

  an object of class `"mcg_cv"`.

- ...:

  ignored.

## Value

An object of class `"mcg_cv"`: a list with `per_fold`, a data frame with
one row per support and fold (columns `support`, `fold`, `score`, the
held-out score, and `npar`, the number of free parameters of the
support, which is the same in every fold); `summary`, a data frame with
one row per support ordered by decreasing score (columns `support`,
`n_folds`, the number of folds with a usable fit, `npar`, `score`, the
mean over those folds, `se`, the standard deviation of the fold scores
divided by the square root of their number, and `delta`, the difference
from the best); `failures`, a data frame of the fits that failed
(`support`, `fold`, `error`); and `folds`, `n` and `supports` (the names
of the candidates). `print` returns its argument invisibly.

## What is scored

The Gaussian working criterion of the residual layer, the scalar the
information criteria of
[`mc_compute_criterion`](https://evandeilton.github.io/mcglmgraph/reference/mc_compute_criterion.md)
are built from, evaluated on the held-out residuals under the covariance
fitted on the rest and divided by the number of held-out units,
\$\$-\tfrac12\bigl\[R\log(2\pi)+\log\|\hat\Sigma\_{\mathrm{train}}\|
+\mathrm{tr}(\hat\Sigma\_{\mathrm{train}}^{-1}\hat
S\_{\mathrm{test}})\bigr\],\$\$ where \\R\\ is the number of responses,
\\\hat\Sigma\_{\mathrm{train}}\\ is the between-response covariance
fitted to the sample covariance of the training residuals, and \\\hat
S\_{\mathrm{test}}\\ is the covariance of the held-out residuals. Each
covariance is centered at the mean of its own part and divided by the
number of units in it. Larger is better. The scale is that of a log
density per unit, so differences between supports are directly
comparable and a difference of a given size means the same at any fold
size.

## What is held fixed

The first stage is not refitted: the split is of the stored Pearson
residuals, so the comparison is between structures for the residual
covariance and not between whole pipelines. That is the question the
sparse against unstructured comparison asks, and it keeps every support
scored on exactly the same residuals.

## See also

[`mc_compute_criterion`](https://evandeilton.github.io/mcglmgraph/reference/mc_compute_criterion.md)
for the in-sample criteria,
[`mc_graph_stepwise`](https://evandeilton.github.io/mcglmgraph/reference/mc_graph_stepwise.md)
and
[`mc_select_lambda`](https://evandeilton.github.io/mcglmgraph/reference/mc_select_lambda.md)
for supports to compare.

Other information criteria:
[`mc_compute_criterion()`](https://evandeilton.github.io/mcglmgraph/reference/mc_compute_criterion.md),
[`mc_graph_ploglik()`](https://evandeilton.github.io/mcglmgraph/reference/mc_graph_ploglik.md)

## Examples

``` r
# \donttest{
set.seed(1); n <- 400
K <- diag(4) * 2
K[1, 2] <- K[2, 1] <- -0.7; K[3, 4] <- K[4, 3] <- -0.7
Y <- matrix(rnorm(n * 4), n) %*% chol(solve(K))
d <- data.frame(y1 = Y[, 1], y2 = Y[, 2], y3 = Y[, 3], y4 = Y[, 4],
                x = rnorm(n))
f <- mcglm_graph(list(y1 ~ x, y2 ~ x, y3 ~ x, y4 ~ x), data = d)
#> Automatic initial values selected.
sel <- mc_graph_stepwise(f, criterion = "bic")
mc_graph_cv(f, supports = list(unstructured = NULL, selected = sel$edges))
#> Out-of-sample comparison of 2 supports (5 folds, n = 400)
#> Held-out working criterion per unit, larger is better.
#> 
#>        support n_folds npar  score      se    delta
#> 1     selected       5    6 -5.516 0.08428 0.000000
#> 2 unstructured       5   10 -5.519 0.08580 0.003232
# }
```
