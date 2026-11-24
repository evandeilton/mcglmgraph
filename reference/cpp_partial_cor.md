# Partial correlation matrix from a precision matrix (C++ kernel)

Entry (i, j) is -K_ij / sqrt(K_ii K_jj), with unit diagonal. The
diagonal of K is assumed positive; this is not checked.

## Usage

``` r
cpp_partial_cor(K)
```

## Arguments

- K:

  symmetric positive-definite precision matrix.

## Value

R x R matrix of partial correlations.
