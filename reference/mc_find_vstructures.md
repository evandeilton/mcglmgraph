# Detect v-structures in the constraint-based skeleton

Runs the PC skeleton search on the estimated between-response covariance
and, for each unshielded triple \\(i, k, j)\\ (\\i - k - j\\ with \\i\\
and \\j\\ nonadjacent), decides whether \\k\\ is a collider \\i \to k
\leftarrow j\\, under one of three rules. A collider is identified by
the absence of \\k\\ from the sets that separate \\i\\ and \\j\\, which
is why the search conditions on subsets and not on all other responses
at once.

## Usage

``` r
mc_find_vstructures(
  fit,
  alpha = 0.05,
  max_ord = NULL,
  df_residual = NULL,
  orient_rule = c("standard", "conservative", "majority")
)
```

## Arguments

- fit:

  an object of class `"mcglm_graph"`.

- alpha:

  significance level of the second-order separation test.

- max_ord:

  maximum conditioning order passed to
  [`mc_pc_skeleton`](https://evandeilton.github.io/mcglmgraph/reference/mc_pc_skeleton.md)
  and to the enumeration of separating sets. A test the skeleton search
  cannot form is skipped with the warning described there.

- df_residual:

  optional residual degrees of freedom passed to
  [`mc_pc_skeleton`](https://evandeilton.github.io/mcglmgraph/reference/mc_pc_skeleton.md)
  and to the reported separating-set tests: `NULL` for the approximation
  by \\n\\, a number, or `"auto"`, which takes \\\nu_e\\ from
  [`mc_residual_df`](https://evandeilton.github.io/mcglmgraph/reference/mc_residual_df.md)
  and refuses a fit for which it is not exact.

- orient_rule:

  `"standard"` (default), `"conservative"` or `"majority"`.

## Value

A data frame with columns `i`, `k`, `j`, `is_collider`, `applied`,
`ambiguous`, `n_sets`, `prop_center`, `pcor` and `p_value`. The last two
refer to the test of \\i\\ and \\j\\ given the set the skeleton
recorded; `applied` is `FALSE` for a collider whose orientation was
withheld because another collider required the opposite direction on one
of its edges. Under `"standard"` the columns `n_sets` and `prop_center`
are `NA`, since no collection is formed.

## The three rules

`"standard"` consults the single set the skeleton recorded for the pair
and declares a collider when the center is absent from it. It is the
rule of Spirtes, Glymour and Scheines and it depends on which separating
set the enumeration happened to reach first. `"conservative"` enumerates
every set that separates the endpoints on the final skeleton and
declares a collider only when the center is in none of them, a
non-collider when it is in all of them, and leaves the triple ambiguous
otherwise (Ramsey, Zhang and Spirtes, 2006). `"majority"` decides by the
proportion of separating sets containing the center, with a tie left
ambiguous (Colombo and Maathuis, 2014). Combined with the
order-independent skeleton, the last two make the procedure independent
of the order of the responses.

## References

Spirtes, P., Glymour, C. and Scheines, R. (2000). *Causation,
Prediction, and Search*, 2nd edition. MIT Press.

Ramsey, J., Zhang, J. and Spirtes, P. (2006). Adjacency-faithfulness and
conservative causal inference. *Proceedings of the 22nd Conference on
Uncertainty in Artificial Intelligence*, 401-408.

Colombo, D. and Maathuis, M. H. (2014). Order-independent
constraint-based causal structure learning. *Journal of Machine Learning
Research* 15, 3921-3962.

## See also

[`mc_pc_skeleton`](https://evandeilton.github.io/mcglmgraph/reference/mc_pc_skeleton.md),
[`mc_separating_sets`](https://evandeilton.github.io/mcglmgraph/reference/mc_separating_sets.md),
[`mc_graph_pipeline`](https://evandeilton.github.io/mcglmgraph/reference/mc_graph_pipeline.md)

Other structure discovery:
[`mc_ges()`](https://evandeilton.github.io/mcglmgraph/reference/mc_ges.md),
[`mc_graph_pipeline()`](https://evandeilton.github.io/mcglmgraph/reference/mc_graph_pipeline.md),
[`mc_orient_meek()`](https://evandeilton.github.io/mcglmgraph/reference/mc_orient_meek.md),
[`mc_pc_skeleton()`](https://evandeilton.github.io/mcglmgraph/reference/mc_pc_skeleton.md),
[`mc_residual_df()`](https://evandeilton.github.io/mcglmgraph/reference/mc_residual_df.md),
[`mc_separating_sets()`](https://evandeilton.github.io/mcglmgraph/reference/mc_separating_sets.md),
[`mc_skeleton_to_bnlearn()`](https://evandeilton.github.io/mcglmgraph/reference/mc_skeleton_to_bnlearn.md)

## Examples

``` r
# \donttest{
set.seed(1); n <- 600
z1 <- rnorm(n); z2 <- rnorm(n); z3 <- 0.7*z1 + 0.7*z2 + rnorm(n)  # collider
d <- data.frame(y1 = z1, y2 = z2, y3 = z3, x = rnorm(n))
f <- mcglm_graph(list(y1~x, y2~x, y3~x), data = d)
#> Automatic initial values selected.
mc_find_vstructures(f)
#>   i k j is_collider applied ambiguous n_sets prop_center       pcor   p_value
#> 1 1 3 2        TRUE    TRUE     FALSE     NA          NA 0.03574128 0.3823036
# }
```
