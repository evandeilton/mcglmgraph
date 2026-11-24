# Greedy equivalence search on the residual covariance

Searches over Markov equivalence classes of directed acyclic graphs
(DAGs) by the two-phase greedy equivalence search (GES) of Chickering
(2002): a forward phase that starts from the empty class (or from
`start`) and applies, at each step, the insertion that raises the score
most, and a backward phase that starts from what the forward phase left
and applies the best deletion. Each phase stops when no operator
improves the score. The score of a class, computed through any of its
members, is the penalized Gaussian working criterion
\$\$\sum\_{j=1}^{R}\left\\-\frac{n}{2}\log\hat\delta_j-
\frac{\lambda}{2}\\(\|P_j\|+1)\right\\,\$\$ where \\\hat\delta_j\\ is
the residual variance of the regression of response \\j\\ on its parents
\\P_j\\ in the residual covariance (the closed-form innovation variance
of the directed layer), \\\|P_j\|+1\\ counts the coefficients and the
innovation variance of the node, and \\\lambda\\ is the `penalty`.
Larger is better. Up to an additive constant the score is \\-1/2\\ times
the criterion \\-2\ell^p_e+\lambda q\\, with \\\ell^p_e\\ the maximized
working criterion and \\q=R+\|E\|\\ the number of parameters of the
member, so that maximizing it with \\\lambda=\log n\\ minimizes the
Bayesian information criterion (BIC) of the working criterion. A parent
set whose block of the residual covariance is numerically singular, or
that leaves a nonpositive residual variance, scores `-Inf` and is never
chosen.

## Usage

``` r
mc_ges(fit, penalty = NULL, start = NULL, max_set = 12L, max_steps = 200L)

# S3 method for class 'mcg_ges'
print(x, ...)
```

## Arguments

- fit:

  an object of class `"mcglm_graph"`, whose residual covariance the
  search runs on.

- penalty:

  single finite non-negative number, the multiplier \\\lambda\\ on the
  parameter count. `NULL` (default) uses `log(n)`, with `n` the number
  of units of the fit, which is the Bayesian information criterion; `2`
  gives Akaike's.

- start:

  optional square 0/1 adjacency matrix on the responses, in the package
  convention, to start the forward phase from; it should be a CPDAG. The
  default is the empty class.

- max_set:

  largest number of candidate neighbors enumerated by an operator (a
  node pair whose candidate set is larger is skipped), a guard against
  the exponential growth in dense neighborhoods. A single non-negative
  integer.

- max_steps:

  cap on the number of steps of each phase, a single positive integer.

- x:

  an object of class `"mcg_ges"`.

- ...:

  ignored.

## Value

An object of class `"mcg_ges"`: a list with the `cpdag` adjacency matrix
(labeled by response); its `score`; the `penalty` and `n` used; `steps`,
a named integer vector with the number of `forward` and `backward`
steps; and `trace`, a data frame with one row per operator applied
(`phase`, `step`, the nodes `x` and `y`, the set `set` of neighbors
involved, comma separated, and the score `gain`). `response_names` is
also recorded. The `print` method returns `x` invisibly.

## Why a search over classes is well posed here

The working criterion of the residual layer depends on the data through
the residual covariance alone, so two members of an equivalence class
fit the same covariance with the same number of parameters and score
alike. That is score equivalence, and without it a search over classes
would not be defined. It also means the output is a class and not a
graph: what the second-order separations identify is what the score can
distinguish. The score is also decomposable, a sum of terms that each
depend on one node and its parents, which is what lets an operator be
scored locally. The large-sample optimality of the search shown by
Chickering (2002) requires in addition a locally consistent score, a
property that is not claimed for the working criterion on Pearson
residuals.

## How it differs from the constraint-based route

[`mc_graph_pipeline`](https://evandeilton.github.io/mcglmgraph/reference/mc_graph_pipeline.md)
tests separations one at a time at a fixed level and then orients; this
one trades every edge against a penalty on the number of parameters and
never tests a hypothesis. The two answer the same question by different
means and need not agree, which is the point of having both. Neither is
run through an external library.

## References

Chickering, D. M. (2002). Optimal structure identification with greedy
search. *Journal of Machine Learning Research* 3, 507-554.

## See also

[`mc_graph_pipeline`](https://evandeilton.github.io/mcglmgraph/reference/mc_graph_pipeline.md)
for the constraint-based route,
[`mc_class_members`](https://evandeilton.github.io/mcglmgraph/reference/mc_class_members.md)
to enumerate the DAGs of the class found.

Other structure discovery:
[`mc_find_vstructures()`](https://evandeilton.github.io/mcglmgraph/reference/mc_find_vstructures.md),
[`mc_graph_pipeline()`](https://evandeilton.github.io/mcglmgraph/reference/mc_graph_pipeline.md),
[`mc_orient_meek()`](https://evandeilton.github.io/mcglmgraph/reference/mc_orient_meek.md),
[`mc_pc_skeleton()`](https://evandeilton.github.io/mcglmgraph/reference/mc_pc_skeleton.md),
[`mc_residual_df()`](https://evandeilton.github.io/mcglmgraph/reference/mc_residual_df.md),
[`mc_separating_sets()`](https://evandeilton.github.io/mcglmgraph/reference/mc_separating_sets.md),
[`mc_skeleton_to_bnlearn()`](https://evandeilton.github.io/mcglmgraph/reference/mc_skeleton_to_bnlearn.md)

## Examples

``` r
# \donttest{
set.seed(1); n <- 600; x <- rnorm(n)
z1 <- rnorm(n); z2 <- rnorm(n); z3 <- 0.8 * z1 + 0.8 * z2 + rnorm(n)
d <- data.frame(y1 = z1, y2 = z2, y3 = z3, x = x)
f <- mcglm_graph(list(y1 ~ x, y2 ~ x, y3 ~ x), data = d)
#> Automatic initial values selected.
mc_ges(f)
#> Greedy equivalence search (penalty = 6.4, n = 600)
#> Steps: 2 forward, 0 backward | score: 270.6351
#> Skeleton edges: 2
#> Adjacency:
#>    y1 y2 y3
#> y1  0  0  1
#> y2  0  0  1
#> y3  0  0  0
# }
```
