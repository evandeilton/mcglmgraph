# PC-style skeleton with separating sets

Constraint-based skeleton discovery from a between-response covariance
matrix. Starting from the complete graph, an edge \\i - j\\ is deleted
as soon as a conditioning set \\S\\ of adjacent nodes makes the partial
correlation \\\rho\_{ij \cdot S}\\ non-significant, and that set is
stored as the separating set of the pair. The neighborhoods are frozen
within each conditioning order and both endpoint neighborhoods are
examined, as in the order-independent skeleton stage of the PC-stable
algorithm (Colombo and Maathuis, 2014), a variant of the PC algorithm
(Spirtes, Glymour and Scheines, 2000; Kalisch and Buehlmann, 2007).

## Usage

``` r
mc_pc_skeleton(Sigma, n, alpha = 0.05, max_ord = NULL, df_residual = NULL)
```

## Arguments

- Sigma:

  an \\R \times R\\ covariance matrix of the responses, e.g.
  `fit$graph$Shat`.

- n:

  number of experimental units, taken as independent.

- alpha:

  significance level of the Fisher z test of a zero partial correlation,
  a single number strictly between zero and one.

- max_ord:

  maximum size of the conditioning set; `NULL` uses \\R - 2\\.

- df_residual:

  optional residual degrees of freedom \\\nu_e\\ from the first-stage
  fit, as a number or as the result of
  [`mc_residual_df`](https://evandeilton.github.io/mcglmgraph/reference/mc_residual_df.md).
  When supplied, the Fisher factor uses \\\nu_e - \|S\| - 2\\; otherwise
  it uses \\n - \|S\| - 3\\.

## Value

A list with the symmetric integer `adjacency` matrix of the skeleton
(one for an edge, zero otherwise); the `sepset` list of length \\R^2\\,
whose elements `(i - 1) * R + j` and `(j - 1) * R + i` hold the integer
vector of nodes that separated \\i\\ from \\j\\, `integer(0)` for the
empty set, and `NULL` when the pair was never separated; the
`df_residual` used; and `diagnostics`, a list with the number of tests
(`n_tests`), the number of ridge-regularized blocks (`n_regularized`),
the number of tests that could not be formed (`n_unusable`, included in
`n_tests`), the largest ridge (`max_ridge`), the smallest Fisher degrees
of freedom among the valid tests (`min_df`) and the contract used
(`df_contract`). The search visits a pair from each of its two
responses, so a conditioning set common to both neighborhoods enters
these counts twice. The two counts are integers, and doubles when they
exceed the integer range.

## The test

The partial correlation \\\hat\rho\_{ij\cdot S}\\ is read off the
inverse of the \\(i, j, S)\\ block of `Sigma` and clamped to
\\\[-0.999999, 0.999999\]\\ before Fisher's transform
\\z=\operatorname{atanh}(\hat\rho\_{ij\cdot S})\sqrt{n-\|S\|-3}\\ (with
\\\nu_e-\|S\|-2\\ in place of \\n-\|S\|-3\\ when `df_residual` is
given), so that a perfectly collinear block gives a large but finite
statistic. A pair is declared conditionally uncorrelated, and its edge
removed, when the two-sided normal p-value of \\z\\ is at least `alpha`.
A block that is nearly singular is inverted after a ridge of \\10^{-8}\\
times its scale and the event is counted in `diagnostics`. A
conditioning set whose Fisher degrees of freedom are not positive, or
whose block gives no valid partial correlation even after the ridge,
gives no test: it is skipped, as in
[`mc_separating_sets`](https://evandeilton.github.io/mcglmgraph/reference/mc_separating_sets.md),
so it removes no edge, and the skipped tests are counted in `n_unusable`
and reported in one warning of class `"mcg_unusable_tests"`. Because
`Sigma` is usually computed from the estimated residuals of responses
that need not be Gaussian, the normal reference for \\z\\ is a working
asymptotic calibration and not an exact one. The factor takes `n` as the
number of independent units, the rows of the data for a fit with an
identity within-response structure. There is no `cluster` argument: with
dependent rows, `n` overstates the information and the test is
anti-conservative.

## References

Spirtes, P., Glymour, C. and Scheines, R. (2000). *Causation,
Prediction, and Search*, 2nd edition. MIT Press.

Kalisch, M. and Buehlmann, P. (2007). Estimating high-dimensional
directed acyclic graphs with the PC-algorithm. *Journal of Machine
Learning Research* 8, 613-636.

Colombo, D. and Maathuis, M. H. (2014). Order-independent
constraint-based causal structure learning. *Journal of Machine Learning
Research* 15, 3921-3962.

## See also

[`mc_separating_sets`](https://evandeilton.github.io/mcglmgraph/reference/mc_separating_sets.md),
[`mc_find_vstructures`](https://evandeilton.github.io/mcglmgraph/reference/mc_find_vstructures.md),
[`mc_graph_pipeline`](https://evandeilton.github.io/mcglmgraph/reference/mc_graph_pipeline.md)

Other structure discovery:
[`mc_find_vstructures()`](https://evandeilton.github.io/mcglmgraph/reference/mc_find_vstructures.md),
[`mc_ges()`](https://evandeilton.github.io/mcglmgraph/reference/mc_ges.md),
[`mc_graph_pipeline()`](https://evandeilton.github.io/mcglmgraph/reference/mc_graph_pipeline.md),
[`mc_orient_meek()`](https://evandeilton.github.io/mcglmgraph/reference/mc_orient_meek.md),
[`mc_residual_df()`](https://evandeilton.github.io/mcglmgraph/reference/mc_residual_df.md),
[`mc_separating_sets()`](https://evandeilton.github.io/mcglmgraph/reference/mc_separating_sets.md),
[`mc_skeleton_to_bnlearn()`](https://evandeilton.github.io/mcglmgraph/reference/mc_skeleton_to_bnlearn.md)

## Examples

``` r
set.seed(1); n <- 600
z1 <- rnorm(n); z2 <- rnorm(n); z3 <- 0.7 * z1 + 0.7 * z2 + rnorm(n)
S <- cov(cbind(z1, z2, z3))
mc_pc_skeleton(S, n)$adjacency
#>      [,1] [,2] [,3]
#> [1,]    0    0    1
#> [2,]    0    0    1
#> [3,]    1    1    0
```
