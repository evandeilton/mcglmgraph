# Covariance from the LDL parametrization (C++ kernel)

Sigma_b = T^-1 D (T^-1)^T with D = diag(delta), symmetrized.

## Usage

``` r
cpp_cholesky_to_cov(gamma, delta, R, edge_idx)
```

## Arguments

- gamma:

  numeric vector of structural coefficients.

- delta:

  numeric vector of R positive innovation variances.

- R:

  integer number of responses.

- edge_idx:

  integer matrix \|E\| x 2 with edges (i, j), 1-based.

## Value

list with Sigma_b, T_inv and T.
