# Apply Meek's orientation rules

Orients every edge whose direction is implied by the edges already
oriented, by repeated application of Meek's four rules until no further
edge can be oriented.

## Usage

``` r
mc_orient_meek(amat, ambiguous = NULL, required = NULL, forbidden = NULL)
```

## Arguments

- amat:

  a square 0/1 adjacency matrix with zero diagonal, in the package
  convention: `amat[i, j] = amat[j, i] = 1` is an undirected edge \\i -
  j\\, `amat[i, j] = 1` with `amat[j, i] = 0` is an arc \\i \to j\\, and
  both zero is no edge.

- ambiguous:

  optional integer matrix with three columns \\(i, k, j)\\, one row per
  unshielded triple whose center could not be decided. Both orders of
  the endpoints are treated as the same triple.

- required:

  optional integer matrix with two columns \\(from, to)\\, arcs the
  orientation must contain.

- forbidden:

  optional integer matrix with two columns \\(from, to)\\, arcs the
  orientation must not contain.

## Value

The oriented adjacency matrix, with attributes `"meek_conflicts"` (edges
two rules wanted in opposite directions, left undirected),
`"meek_iterations"`, `"meek_blocked"` (the number of orientations
withheld by an ambiguous triple) and `"knowledge"`, a list recording the
arcs imposed and whether the result stayed acyclic and kept every
required arc.

## The four rules

\$\$\begin{array}{ll} \mathrm{R1} & c \to a,\\ a - b,\\ c\\
\mathrm{and}\\ b\\ \mathrm{nonadjacent} \Rightarrow a \to b\\
\mathrm{R2} & a \to c \to b,\\ a - b \Rightarrow a \to b\\ \mathrm{R3} &
a - b,\\ a - c,\\ a - d,\\ c \to b,\\ d \to b,\\ c\\ \mathrm{and}\\ d\\
\mathrm{nonadjacent} \Rightarrow a \to b\\ \mathrm{R4} & a - b,\\ a -
c,\\ a - d,\\ c \to d,\\ d \to b,\\ c\\ \mathrm{and}\\ b\\
\mathrm{nonadjacent} \Rightarrow a \to b \end{array}\$\$ R1 to R3 are
complete when the only oriented edges are those of the v-structures, so
the completed partially directed acyclic graph (CPDAG) of a Markov
equivalence class needs no more than those three (Meek, 1995). R4 is
idle in that case and becomes necessary as soon as edges are oriented
from outside, through `required` or `forbidden`, because the graph is
then no longer the CPDAG of any class and the first three rules stop
being complete. For an input that admits a consistent extension, the
four together return the maximally oriented partially directed graph: an
edge is left undirected exactly when the background knowledge and the
v-structures admit an acyclic extension with each of its two directions.

## Background knowledge

`required` lists arcs that must be present with the stated direction and
`forbidden` arcs that must not. An undirected edge whose reverse is
forbidden is oriented, since only one direction is left. Both are
applied before propagation and are checked for internal consistency: an
arc required in both directions, an arc both required and forbidden, an
arc required between responses the skeleton does not join, and an arc
required against an edge already oriented the other way are all refused.
An incompatibility that only surfaces during propagation, such as a
cycle, is reported in the returned attributes rather than raised,
because the partially oriented graph is still informative about where
the incompatibility lies.

## Ambiguous triples

Under the conservative and majority rules an unshielded triple can be
left undecided, and propagation must not settle by implication what the
separations left open. R1, R3 and R4 all orient an edge in order to
avoid creating a collider at its tail, so each is withheld when the
triple carrying that collider is ambiguous. R2 only forbids a cycle and
involves no unshielded triple, so it is applied unconditionally. Leaving
`ambiguous` at `NULL` reproduces the unrestricted rules exactly.

## References

Meek, C. (1995). Causal inference and causal explanation with background
knowledge. *Proceedings of the 11th Conference on Uncertainty in
Artificial Intelligence*, 403-410.

Verma, T. and Pearl, J. (1990). Equivalence and synthesis of causal
models. *Proceedings of the 6th Conference on Uncertainty in Artificial
Intelligence*, 220-227.

## See also

Other structure discovery:
[`mc_find_vstructures()`](https://evandeilton.github.io/mcglmgraph/reference/mc_find_vstructures.md),
[`mc_ges()`](https://evandeilton.github.io/mcglmgraph/reference/mc_ges.md),
[`mc_graph_pipeline()`](https://evandeilton.github.io/mcglmgraph/reference/mc_graph_pipeline.md),
[`mc_pc_skeleton()`](https://evandeilton.github.io/mcglmgraph/reference/mc_pc_skeleton.md),
[`mc_residual_df()`](https://evandeilton.github.io/mcglmgraph/reference/mc_residual_df.md),
[`mc_separating_sets()`](https://evandeilton.github.io/mcglmgraph/reference/mc_separating_sets.md),
[`mc_skeleton_to_bnlearn()`](https://evandeilton.github.io/mcglmgraph/reference/mc_skeleton_to_bnlearn.md)

## Examples

``` r
# 1 -> 3 <- 2 with 3 - 4: R1 orients 3 -> 4
a <- matrix(0L, 4, 4)
a[1, 3] <- a[2, 3] <- 1L
a[3, 4] <- a[4, 3] <- 1L
mc_orient_meek(a)
#>      [,1] [,2] [,3] [,4]
#> [1,]    0    0    1    0
#> [2,]    0    0    1    0
#> [3,]    0    0    0    1
#> [4,]    0    0    0    0
#> attr(,"meek_conflicts")
#> [1] from to  
#> <0 rows> (or 0-length row.names)
#> attr(,"meek_iterations")
#> [1] 2
#> attr(,"meek_blocked")
#> [1] 0
#> attr(,"knowledge")
#> attr(,"knowledge")$required
#>      [,1] [,2]
#> 
#> attr(,"knowledge")$forbidden
#>      [,1] [,2]
#> 
#> attr(,"knowledge")$n_imposed
#> [1] 0
#> 
#> attr(,"knowledge")$acyclic
#> [1] TRUE
#> 
#> attr(,"knowledge")$required_kept
#> [1] TRUE
#> 

# background knowledge on a fully undirected chain 1 - 2 - 3
b <- matrix(0L, 3, 3)
b[1, 2] <- b[2, 1] <- b[2, 3] <- b[3, 2] <- 1L
mc_orient_meek(b, required = rbind(c(1L, 2L)))   # 1 -> 2 forces 2 -> 3
#>      [,1] [,2] [,3]
#> [1,]    0    1    0
#> [2,]    0    0    1
#> [3,]    0    0    0
#> attr(,"meek_conflicts")
#> [1] from to  
#> <0 rows> (or 0-length row.names)
#> attr(,"meek_iterations")
#> [1] 2
#> attr(,"meek_blocked")
#> [1] 0
#> attr(,"knowledge")
#> attr(,"knowledge")$required
#>      [,1] [,2]
#> [1,]    1    2
#> 
#> attr(,"knowledge")$forbidden
#>      [,1] [,2]
#> 
#> attr(,"knowledge")$n_imposed
#> [1] 1
#> 
#> attr(,"knowledge")$acyclic
#> [1] TRUE
#> 
#> attr(,"knowledge")$required_kept
#> [1] TRUE
#> 
```
