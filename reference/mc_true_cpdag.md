# CPDAG of a known directed structure

Builds the completed partially directed acyclic graph (CPDAG) of the
equivalence class of a directed acyclic graph (DAG) given through the
structural matrix of the modified Cholesky parametrization, where a
nonzero \\T\_{ji}\\ codes the arc \\i \to j\\. The equivalence class is
the set of DAGs with the same skeleton and the same v-structures (Verma
and Pearl, 1990; Andersson, Madigan and Perlman, 1997). The v-structures
are written first and Meek's rules complete the orientation, which gives
the object an estimated CPDAG should be compared against.

## Usage

``` r
mc_true_cpdag(T, tol = 0)
```

## Arguments

- T:

  a finite numeric square matrix whose entry \\(j, i)\\ is nonzero when
  \\i \to j\\; the structural matrix of the directed layer is one such,
  and the transpose of a directed adjacency is another.

- tol:

  non-negative number; entries not above it in absolute value count as
  zero (default 0).

## Value

The adjacency matrix of the CPDAG, in the package convention (see
[`mc_cpdag_metrics`](https://evandeilton.github.io/mcglmgraph/reference/mc_cpdag_metrics.md)),
with the dimnames of `T`.

## Arcs anywhere, not only below the diagonal

A structural matrix of the modified Cholesky parametrization is lower
triangular, so its arcs run from a lower index to a higher one and
reading only the lower triangle would suffice. A DAG that arrives from
elsewhere need not be numbered in a topological order, and its arcs can
sit anywhere, so every off-diagonal entry is read. For a triangular
input the two readings agree. A pair carrying an arc in both directions,
or a set of arcs with a cycle, is refused.

## References

Verma, T. and Pearl, J. (1990). Equivalence and synthesis of causal
models. *Proceedings of the Sixth Conference on Uncertainty in
Artificial Intelligence*, 220–227.

Andersson, S. A., Madigan, D. and Perlman, M. D. (1997). A
characterization of Markov equivalence classes for acyclic digraphs.
*The Annals of Statistics* 25(2), 505–541.
[doi:10.1214/aos/1031833662](https://doi.org/10.1214/aos/1031833662)

## See also

[`mc_cpdag_metrics`](https://evandeilton.github.io/mcglmgraph/reference/mc_cpdag_metrics.md).

Other structure recovery:
[`mc_cpdag_metrics()`](https://evandeilton.github.io/mcglmgraph/reference/mc_cpdag_metrics.md),
[`mc_structure_metrics()`](https://evandeilton.github.io/mcglmgraph/reference/mc_structure_metrics.md)

## Examples

``` r
T <- diag(3)
T[3, 1] <- -0.5; T[3, 2] <- -0.4        # 1 -> 3 <- 2, a collider
mc_true_cpdag(T)
#>      [,1] [,2] [,3]
#> [1,]    0    0    1
#> [2,]    0    0    1
#> [3,]    0    0    0
```
