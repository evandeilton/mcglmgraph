# Cross-sensitivity factor of the graph layer (C++ kernel)

Forms the shortcut factor Q0(d, r) = (1/2) \[ K Ddot_d \](r, r), which
is the reduced form factor (1/2) \[ K Ddot_d K S_e \](r, r) evaluated at
the fitted covariance (S_e replaced by the fitted Sigma_b = K^-1),
before the division by the average squared residual m_r of each
response.

## Usage

``` r
cpp_two_stage_q(K, Dlist)
```

## Arguments

- K:

  the R x R precision matrix of the fitted layer.

- Dlist:

  list of P derivatives dSigma_b/dtheta_d, each R x R.

## Value

the P x R matrix Q0, whose entry (d, r) is the factor of parameter d and
response r before division by m_r.
