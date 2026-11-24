# Per-unit Pearson scores of the directed layer (C++ kernel)

The unit contribution of the coefficients of node j is eps_ij r\[i,
P_j\] / delta_j, with eps_ij = r\[i, j\] - gamma_j' r\[i, P_j\], and
that of its innovation variance is (eps_ij^2 - delta_j) / (2 delta_j^2).
Summed over units these are the estimating functions whose root is
cpp_dag_closed_root(), and their outer product is the variability of the
sandwich.

## Usage

``` r
cpp_dag_scores(Ec, gamma, delta, edge_idx)
```

## Arguments

- Ec:

  n x R matrix of centered Pearson residuals.

- gamma:

  coefficients in the order of edge_idx.

- delta:

  vector of R innovation variances.

- edge_idx:

  integer \|E\| x 2 matrix of edges, 1-based.

## Value

list with scores (n x (\|E\| + R), in the order (gamma, delta)) and
innovations (n x R).
