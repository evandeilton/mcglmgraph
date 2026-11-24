# Inverse of the block-diagonal sensitivity of the directed layer (C++ kernel)

The sensitivity of the directed layer is block diagonal with one block
per node, so its inverse is too: the block of node j is delta_j
Sigma\[P_j, P_j\]^-1 / n and the entry of its innovation variance is 2
delta_j^2 / n. Inverting block by block costs the sum of the cubes of
the parent-set sizes instead of the cube of the number of parameters.

## Usage

``` r
cpp_dag_sensitivity_inverse(Sigma, delta, n, edge_idx)
```

## Arguments

- Sigma:

  R x R fitted covariance of the layer.

- delta:

  vector of R innovation variances.

- n:

  number of experimental units.

- edge_idx:

  integer \|E\| x 2 matrix of edges, 1-based.

## Value

list with the inverse and n_singular, the number of nodes whose parent
block could not be inverted; a nonzero value means the caller must fall
back on the general route.
