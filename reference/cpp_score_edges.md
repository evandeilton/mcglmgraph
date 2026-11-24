# Score statistics for candidate (absent) edges (C++ kernel)

For each candidate edge e=(i,j) of the linear-precision parametrization,
the Pearson score is psi_e = n (Sigma - Shat)\[i,j\] and the sensitivity
is J_ee = n (Sigma_ii Sigma_jj + Sigma_ij^2); the statistic is S_e =
psi_e^2 / J_ee, with S_e = 0 when J_ee is not positive. It ignores the
estimation of the other parameters, so it ranks candidates but is not
calibrated as chi-square(1); the nuisance-adjusted statistic is formed
in R by mc_edge_scores().

## Usage

``` r
cpp_score_edges(Sigma, Shat, n, edges)
```

## Arguments

- Sigma:

  R x R fitted covariance matrix.

- Shat:

  R x R sample covariance of the Pearson residuals.

- n:

  number of units.

- edges:

  integer matrix q x 2 of candidate edges (i, j), 1-based.

## Value

numeric vector of score statistics.
