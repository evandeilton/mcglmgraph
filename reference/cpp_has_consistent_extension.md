# Does a consistent DAG extension exist? (C++ kernel)

The criterion of Dor and Tarsi (1992): repeatedly remove a node with no
outgoing arc whose undirected neighbors are adjacent to all its other
neighbors. A graph that empties out admits an extension preserving the
arcs and the v-structures.

## Usage

``` r
cpp_has_consistent_extension(amat)
```

## Arguments

- amat:

  integer adjacency in the package convention.

## Value

TRUE when an extension exists.
