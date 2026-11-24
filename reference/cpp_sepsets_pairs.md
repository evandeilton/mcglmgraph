# All separating sets of given pairs of nodes (C++ kernel)

For each requested pair (i, j), enumerates every subset of the neighbors
of i other than j, and every subset of the neighbors of j other than i,
of size at most max_ord, and returns those that make the second-order
partial correlation non-significant at level alpha. Sets found in both
neighborhoods are returned once.

## Usage

``` r
cpp_sepsets_pairs(Sigma, n, adj, pairs, alpha, max_ord, df_residual = 0)
```

## Arguments

- Sigma:

  R x R covariance matrix of the responses.

- n:

  number of experimental units.

- adj:

  integer R x R symmetric 0/1 adjacency matrix of the skeleton.

- pairs:

  integer q x 2 matrix of pairs (i, j), 1-based, assumed non-adjacent in
  adj.

- alpha:

  significance level of the conditional-association test.

- max_ord:

  largest conditioning set size to enumerate.

- df_residual:

  residual degrees of freedom of the first stage; a value not greater
  than zero selects the raw contract n - \|S\| - 3.

## Value

list with set_pair (1-based row of pairs), set_members (n_sets x R,
0/1), set_size, set_pcor, set_p_value, and the diagnostics n_tests,
n_regularized, max_ridge, n_unusable and truncated (one flag per pair).
