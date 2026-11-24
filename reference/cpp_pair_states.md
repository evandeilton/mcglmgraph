# State of every pair of nodes (C++ kernel)

State of every pair of nodes (C++ kernel)

## Usage

``` r
cpp_pair_states(amat)
```

## Arguments

- amat:

  integer adjacency in the package convention.

## Value

integer vector over the pairs in upper-triangular order: 1 absent, 2
undirected, 3 i to j, 4 j to i.
