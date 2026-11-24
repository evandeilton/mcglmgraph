# Analytical derivatives of Sigma_b with respect to eta (C++ kernel)

Since Sigma_b = K_b^-1 and dK_b/deta_d = Z_d, the derivative is
dSigma_b/deta_d = -Sigma_b Z_d Sigma_b; each matrix is symmetrized.

## Usage

``` r
cpp_precision_derivatives(Sigma_b, bases)
```

## Arguments

- Sigma_b:

  the R x R covariance matrix.

- bases:

  numeric cube R x R x D of base matrices Z_d.

## Value

list of D symmetric matrices, each -Sigma_b Z_d Sigma_b.
