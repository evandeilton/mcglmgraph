# V-structures of a partially directed graph (C++ kernel)

Triples i -\> k \<- j with i \< j and i, j nonadjacent, each arc
directed.

## Usage

``` r
cpp_vstructure_set(amat)
```

## Arguments

- amat:

  integer adjacency in the package convention.

## Value

integer matrix with columns (i, k, j), 1-based.
