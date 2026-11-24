# Analytical LDL derivatives of Sigma_b (C++ kernel)

The derivative in the innovation variance delta_k is T^-1 E_kk (T^-1)^T,
and the derivative in the structural coefficient gamma_ji is T^-1 E_ji
Sigma_b + Sigma_b E_ij (T^-1)^T, which is symmetric because the second
term is the transpose of the first. E_ab has a one in position (a, b)
and zeros elsewhere.

## Usage

``` r
cpp_cholesky_derivatives(gamma, delta, R, edge_idx)
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

list with D_delta (R matrices) and D_gamma (\|E\| matrices).
