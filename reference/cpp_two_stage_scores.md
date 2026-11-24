# Corrected per-unit contributions of the graph layer (C++ kernel)

Returns c_i = psi_i - Q (r_i o r_i - m), with m the observed average of
the squared residuals of each response and Q(d, r) = (1/2) \[ K Ddot_d
\](r, r) / m_r. Summed over units these vanish, and their outer product
is the middle of the corrected sandwich.

## Usage

``` r
cpp_two_stage_scores(Psi, K, Dlist, Ec)
```

## Arguments

- Psi:

  n x P matrix of Pearson contributions.

- K:

  the R x R precision matrix of the fitted layer.

- Dlist:

  list of P derivatives dSigma_b/dtheta_d.

- Ec:

  n x R matrix of centered Pearson residuals; m is the column mean of
  its squares, which must be positive.

## Value

list with `scores` (n x P), the factor `Q` (P x R) applied, already
divided by m_r, and `center`, the vector of the R averages m.
