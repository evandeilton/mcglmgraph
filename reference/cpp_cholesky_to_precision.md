# Precision from the LDL parametrization (C++ kernel)

Precision from the LDL parametrization (C++ kernel)

## Usage

``` r
cpp_cholesky_to_precision(gamma, delta, R, edge_idx)
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

the R x R precision matrix K_b = T^T D^-1 T.
