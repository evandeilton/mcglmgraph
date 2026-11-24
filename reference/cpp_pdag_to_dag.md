# A DAG consistent with a partially directed graph (C++ kernel)

The algorithm of Dor and Tarsi (1992): repeatedly take a node with no
outgoing arc whose undirected neighbors are adjacent to all its other
neighbors, orient every undirected edge at it inwards, and remove it.

## Usage

``` r
cpp_pdag_to_dag(amat)
```

## Arguments

- amat:

  the partially directed graph.

## Value

list with the DAG and a flag saying whether an extension exists.
