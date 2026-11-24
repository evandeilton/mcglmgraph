# Pearson estimating function and sensitivity for Sigma_b (C++ kernel)

For n i.i.d. units with sample covariance Shat, computes the
second-order Pearson estimating function psi and the sensitivity J for
the covariance parameters, given the analytical derivatives Dlist of
Sigma: psi_d = -0.5 n \[tr(S^-1 D_d) - tr(S^-1 D_d S^-1 Shat)\], J_ab =
0.5 n tr(S^-1 D_a S^-1 D_b), with S = Sigma. Sigma must be symmetric
positive definite.

## Usage

``` r
cpp_pearson_between(Sigma, Shat, n, Dlist)
```

## Arguments

- Sigma:

  current R x R symmetric positive definite covariance matrix.

- Shat:

  R x R sample covariance of the Pearson residuals.

- n:

  number of units.

- Dlist:

  list of P derivative matrices dSigma/dtheta_d.

## Value

list with psi (length P) and J (P x P).
