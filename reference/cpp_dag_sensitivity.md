# Block-diagonal sensitivity of the directed layer (C++ kernel)

Assembles J_n on the accumulated scale from its closed form: the block
of the coefficients of node j is n Sigma\[P_j, P_j\] / delta_j, that of
its innovation variance is n / (2 delta_j^2), and every other block is
zero, including the one between coefficients of distinct nodes.

## Usage

``` r
cpp_dag_sensitivity(Sigma, delta, n, edge_idx)
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

the (\|E\| + R) x (\|E\| + R) sensitivity in the parameter order (gamma,
delta).
