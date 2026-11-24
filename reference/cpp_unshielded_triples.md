# Enumerate unshielded triples of an undirected graph (C++ kernel)

Returns triples (i, k, j) such that i-k and k-j are edges but i-j is
not, the candidate v-structures of CPDAG discovery.

## Usage

``` r
cpp_unshielded_triples(adj)
```

## Arguments

- adj:

  integer adjacency matrix (R x R, 0/1, symmetric).

## Value

integer matrix with one row per triple, columns (i, k, j), 1-based.
