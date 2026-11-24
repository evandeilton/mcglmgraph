# All separating sets of a pair of responses

Enumerates every conditioning set that makes a pair of responses
conditionally uncorrelated at the second order, given a skeleton. The
skeleton search of
[`mc_pc_skeleton`](https://evandeilton.github.io/mcglmgraph/reference/mc_pc_skeleton.md)
stops at the first such set and records only that one, which is all the
standard v-structure rule needs. The order-independent orientation rules
need the whole collection, and they need it on the final skeleton rather
than on the neighborhoods that happened to be live when the edge was
deleted, so it is rebuilt here from scratch: for the pair \\(i, j)\\
every subset of the neighbors of \\i\\ other than \\j\\, and every
subset of the neighbors of \\j\\ other than \\i\\, of size at most
`max_ord`, is tested, and the separating ones are returned. A set
reachable from both neighborhoods is reported once.

## Usage

``` r
mc_separating_sets(
  Sigma,
  n,
  adjacency,
  pairs = NULL,
  alpha = 0.05,
  max_ord = NULL,
  df_residual = NULL
)

# S3 method for class 'mcg_sepsets'
print(x, ...)
```

## Arguments

- Sigma:

  an \\R \times R\\ covariance matrix of the responses.

- n:

  number of experimental units, taken as independent; the test is the
  one of
  [`mc_pc_skeleton`](https://evandeilton.github.io/mcglmgraph/reference/mc_pc_skeleton.md),
  with the same clamp of the partial correlation and no `cluster`
  argument.

- adjacency:

  the symmetric 0/1 skeleton, for instance the `adjacency` component of
  [`mc_pc_skeleton`](https://evandeilton.github.io/mcglmgraph/reference/mc_pc_skeleton.md).

- pairs:

  optional integer matrix with two columns giving the pairs to examine;
  `NULL` takes every non-adjacent pair. The pairs must be non-adjacent
  in `adjacency`, since a separating set is not defined for an edge of
  the skeleton.

- alpha:

  significance level of the conditional-association test, on the same
  contract as
  [`mc_pc_skeleton`](https://evandeilton.github.io/mcglmgraph/reference/mc_pc_skeleton.md).

- max_ord:

  largest conditioning set enumerated; `NULL` uses \\R - 2\\.

- df_residual:

  optional residual degrees of freedom \\\nu_e\\ of the first-stage fit,
  a number or the result of
  [`mc_residual_df`](https://evandeilton.github.io/mcglmgraph/reference/mc_residual_df.md),
  passed to the same Fisher contract used by the skeleton.

- x:

  an object of class `"mcg_sepsets"`.

- ...:

  ignored.

## Value

An object of class `"mcg_sepsets"`: a list with `pairs` (one row per
pair examined, with the number of separating sets found and a
`truncated` flag), `sets` (one row per separating set, with its label,
size, partial correlation and p-value), the 0/1 matrix `members` of set
membership, and `diagnostics`. The `truncated` flag of a pair is `TRUE`
when its collection may be incomplete: a neighborhood of the pair is
larger than `max_ord`, or some conditioning set gave no valid test,
because its Fisher degrees of freedom are not positive or its block
gives no valid partial correlation. Such a set is skipped and counted in
`n_unusable`; unlike in
[`mc_pc_skeleton`](https://evandeilton.github.io/mcglmgraph/reference/mc_pc_skeleton.md),
it does not stop the enumeration. The object also records `alpha`,
`max_ord`, `df_residual` and `response_names`, and `diagnostics` with
the number of tests (`n_tests`), of ridge-regularized blocks
(`n_regularized`), of unusable blocks (`n_unusable`), the largest ridge
(`max_ridge`) and the Fisher contract (`df_contract`). The `print`
method returns `x` invisibly.

## Why the collection is informative on its own

Two sets can separate the same pair with very different strength, and
the skeleton keeps whichever the enumeration reached first. Reporting
the collection shows how much of the orientation rests on that arbitrary
choice.

## References

Ramsey, J., Zhang, J. and Spirtes, P. (2006). Adjacency-faithfulness and
conservative causal inference. *Proceedings of the 22nd Conference on
Uncertainty in Artificial Intelligence*, 401-408.

Colombo, D. and Maathuis, M. H. (2014). Order-independent
constraint-based causal structure learning. *Journal of Machine Learning
Research* 15, 3921-3962.

## See also

[`mc_graph_pipeline`](https://evandeilton.github.io/mcglmgraph/reference/mc_graph_pipeline.md),
whose `orient_rule` argument consumes this object.

Other structure discovery:
[`mc_find_vstructures()`](https://evandeilton.github.io/mcglmgraph/reference/mc_find_vstructures.md),
[`mc_ges()`](https://evandeilton.github.io/mcglmgraph/reference/mc_ges.md),
[`mc_graph_pipeline()`](https://evandeilton.github.io/mcglmgraph/reference/mc_graph_pipeline.md),
[`mc_orient_meek()`](https://evandeilton.github.io/mcglmgraph/reference/mc_orient_meek.md),
[`mc_pc_skeleton()`](https://evandeilton.github.io/mcglmgraph/reference/mc_pc_skeleton.md),
[`mc_residual_df()`](https://evandeilton.github.io/mcglmgraph/reference/mc_residual_df.md),
[`mc_skeleton_to_bnlearn()`](https://evandeilton.github.io/mcglmgraph/reference/mc_skeleton_to_bnlearn.md)

## Examples

``` r
set.seed(1); n <- 600
z1 <- rnorm(n); z2 <- rnorm(n); z3 <- 0.7 * z1 + 0.7 * z2 + rnorm(n)
S <- cov(cbind(z1, z2, z3))
sk <- mc_pc_skeleton(S, n)
mc_separating_sets(S, n, sk$adjacency)$sets
#>   pair i j response1 response2 set size       pcor   p_value
#> 1    1 1 2        z1        z2  {}    0 0.03853124 0.3462323
```
