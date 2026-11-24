# Members of a Markov equivalence class

Enumerates the directed acyclic graphs (DAGs) that a completed partially
directed acyclic graph (CPDAG) stands for: the acyclic orientations of
its undirected edges that keep every directed edge and introduce no
v-structure (immorality) the class does not have. By the
characterization of Markov equivalence by skeleton and v-structures
(Verma and Pearl, 1990; Andersson, Madigan and Perlman, 1997), these are
exactly the DAGs of the equivalence class. The members imply the same
separations, so observational data that carry only those constraints do
not distinguish among them. For a partially directed graph that is not a
CPDAG, such as one oriented further by background knowledge, the same
enumeration returns its consistent DAG extensions.

## Usage

``` r
mc_class_members(cpdag, max_members = 1000L)
```

## Arguments

- cpdag:

  a square 0/1 adjacency matrix in the package convention (an undirected
  edge has both entries one, an arc \\i \to j\\ only `[i, j]`), or an
  object of class `"mcg_cpdag"`.

- max_members:

  a single positive integer; the enumeration stops after this many
  members and reports the truncation.

## Value

A list of adjacency matrices, one per member, each a DAG in the same
convention and with the dimnames of `cpdag`. The list carries the
attribute `"truncated"`, `TRUE` when `max_members` was reached before
the enumeration ended. The list is empty when the graph admits no
consistent extension.

## How they are enumerated

One undirected edge is oriented at a time, both ways, Meek's rules
propagate the consequences, and a branch is abandoned as soon as it
creates a cycle or an immorality the class does not have. Because the
two branches fix opposite directions on the same edge, no member is
produced twice. The count grows quickly with the number of undirected
edges, so `max_members` caps it and the result says whether the cap was
reached.

## References

Verma, T. and Pearl, J. (1990). Equivalence and synthesis of causal
models. *Proceedings of the 6th Conference on Uncertainty in Artificial
Intelligence*, 220-227.

Andersson, S. A., Madigan, D. and Perlman, M. D. (1997). A
characterization of Markov equivalence classes for acyclic digraphs.
*The Annals of Statistics* 25(2), 505-541.
[doi:10.1214/aos/1031833662](https://doi.org/10.1214/aos/1031833662)

Meek, C. (1995). Causal inference and causal explanation with background
knowledge. *Proceedings of the 11th Conference on Uncertainty in
Artificial Intelligence*, 403-410.

## See also

[`mc_class_coefficients`](https://evandeilton.github.io/mcglmgraph/reference/mc_class_coefficients.md)
for what the members disagree about,
[`mc_orient_meek`](https://evandeilton.github.io/mcglmgraph/reference/mc_orient_meek.md)
for the propagation rules.

Other equivalence class:
[`mc_class_coefficients()`](https://evandeilton.github.io/mcglmgraph/reference/mc_class_coefficients.md),
[`mc_sid()`](https://evandeilton.github.io/mcglmgraph/reference/mc_sid.md)

## Examples

``` r
# the chain 1 - 2 - 3 has three members, the collider none of them
a <- matrix(0L, 3, 3)
a[1, 2] <- a[2, 1] <- a[2, 3] <- a[3, 2] <- 1L
length(mc_class_members(a))
#> [1] 3
```
