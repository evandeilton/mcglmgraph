# Recovery of a partially directed structure

Compares an estimated completed partially directed acyclic graph (CPDAG)
with a known one, separating what the skeleton recovers from what the
orientation does.

## Usage

``` r
mc_cpdag_metrics(estimated, truth)
```

## Arguments

- estimated, truth:

  square adjacency matrices of the same dimension with zeros and ones
  and a zero diagonal, in the package convention: entries \\(i,j)\\ and
  \\(j,i)\\ both set mean an undirected edge, and only \\(i,j)\\ set
  means the arc \\i \to j\\.

## Value

A named numeric vector with `cpdag_shd`, `cpdag_skeleton_f1`,
`cpdag_skeleton_shd`, `cpdag_vstructure_recall`,
`cpdag_vstructure_precision`, `cpdag_orient_recall` and
`cpdag_overorient`.

## What each measure counts

Every pair of responses is in one of four states: absent, undirected, or
directed either way. `cpdag_shd` counts the pairs whose states differ.
The skeleton measures (`cpdag_skeleton_f1`, `cpdag_skeleton_shd`) are
those of
[`mc_structure_metrics`](https://evandeilton.github.io/mcglmgraph/reference/mc_structure_metrics.md)
and ignore direction. `cpdag_orient_recall` is the proportion of the
pairs the truth directs that the estimate carries in the same state, and
`cpdag_overorient` the proportion of the pairs the truth leaves
undirected while joining them that the estimate directs, which is the
rate at which a direction is invented. A v-structure is a triple \\i \to
k \leftarrow j\\ with \\i\\ and \\j\\ nonadjacent. V-structure recall is
the proportion of the true v-structures that the estimate declares, and
precision the proportion of the declared ones that are true; recall is
conditional on the truth containing one and precision on some
v-structure having been declared, and either is `NA` when its condition
is empty.

## See also

[`mc_true_cpdag`](https://evandeilton.github.io/mcglmgraph/reference/mc_true_cpdag.md)
to build the true CPDAG of a DAG,
[`mc_structure_metrics`](https://evandeilton.github.io/mcglmgraph/reference/mc_structure_metrics.md).

Other structure recovery:
[`mc_structure_metrics()`](https://evandeilton.github.io/mcglmgraph/reference/mc_structure_metrics.md),
[`mc_true_cpdag()`](https://evandeilton.github.io/mcglmgraph/reference/mc_true_cpdag.md)

## Examples

``` r
truth <- matrix(0L, 3, 3)
truth[1, 3] <- truth[2, 3] <- 1L          # 1 -> 3 <- 2
est <- truth
mc_cpdag_metrics(est, truth)
#>                  cpdag_shd          cpdag_skeleton_f1 
#>                          0                          1 
#>         cpdag_skeleton_shd    cpdag_vstructure_recall 
#>                          0                          1 
#> cpdag_vstructure_precision        cpdag_orient_recall 
#>                          1                          1 
#>           cpdag_overorient 
#>                         NA 
```
