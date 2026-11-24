# Convert a skeleton to the set of allowed bnlearn arcs

Lists both directions of every edge of a skeleton as named arcs, in the
`from`/`to` format of bnlearn. The complement of this set within all
ordered pairs of nodes is a blacklist that confines a score-based search
to the skeleton.

## Usage

``` r
mc_skeleton_to_bnlearn(edges, node_names)
```

## Arguments

- edges:

  two-column integer matrix of undirected skeleton edges, as indices
  into `node_names`; `NULL` or a matrix with no rows gives an empty
  result.

- node_names:

  character vector of node labels.

## Value

A data frame with character columns `from` and `to`, containing both
directions of every skeleton edge. Its complement can be supplied to
bnlearn as a blacklist.

## See also

Other structure discovery:
[`mc_find_vstructures()`](https://evandeilton.github.io/mcglmgraph/reference/mc_find_vstructures.md),
[`mc_ges()`](https://evandeilton.github.io/mcglmgraph/reference/mc_ges.md),
[`mc_graph_pipeline()`](https://evandeilton.github.io/mcglmgraph/reference/mc_graph_pipeline.md),
[`mc_orient_meek()`](https://evandeilton.github.io/mcglmgraph/reference/mc_orient_meek.md),
[`mc_pc_skeleton()`](https://evandeilton.github.io/mcglmgraph/reference/mc_pc_skeleton.md),
[`mc_residual_df()`](https://evandeilton.github.io/mcglmgraph/reference/mc_residual_df.md),
[`mc_separating_sets()`](https://evandeilton.github.io/mcglmgraph/reference/mc_separating_sets.md)

## Examples

``` r
mc_skeleton_to_bnlearn(rbind(c(1, 2), c(2, 3)), c("y1", "y2", "y3"))
#>   from to
#> 1   y1 y2
#> 2   y2 y3
#> 3   y2 y1
#> 4   y3 y2
```
