# Pearson estimating function of the directed layer in closed form (C++ kernel)

Evaluates the estimating function on the accumulated scale directly from
the sample covariance, without forming the generic traces: psi of
gamma_j = (n/delta_j) (S\[P_j, j\] - S\[P_j, P_j\] gamma_j), psi of
delta_j = (n/2) (q_j - delta_j) / delta_j^2, with q_j = S\[j, j\] - 2
gamma_j' S\[P_j, j\] + gamma_j' S\[P_j, P_j\] gamma_j. It vanishes at
the root of cpp_dag_closed_root() and is what certifies it.

## Usage

``` r
cpp_dag_pearson_closed(Shat, gamma, delta, n, edge_idx)
```

## Arguments

- Shat:

  R x R sample covariance, on the ordering of the support.

- gamma:

  coefficients in the order of edge_idx.

- delta:

  vector of R innovation variances.

- n:

  number of experimental units.

- edge_idx:

  integer \|E\| x 2 matrix of edges, 1-based.

## Value

the estimating function in the order (gamma, delta).
