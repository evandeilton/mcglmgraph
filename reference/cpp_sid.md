# Structural intervention distance between two DAGs (C++ kernel)

Counts the ordered pairs (i, j), i not j, for which the parents of i in
the estimated graph are not a valid adjustment set for the effect of i
on j in the true graph.

## Usage

``` r
cpp_sid(Gtrue, Hest)
```

## Arguments

- Gtrue:

  integer R x R adjacency of the true DAG, one in (a, b) meaning a -\>
  b.

- Hest:

  integer R x R adjacency of the estimated DAG, same convention.

## Value

list with the count and the R x R matrix of offending pairs.
