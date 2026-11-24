# Invert a symmetric positive-definite precision matrix (C++ kernel)

Inverts the symmetric part of K by Cholesky-based inversion and stops
with an error when K is not positive definite.

## Usage

``` r
cpp_precision_to_cov(K)
```

## Arguments

- K:

  symmetric positive-definite matrix.

## Value

The covariance matrix Sigma_b = K^-1.
