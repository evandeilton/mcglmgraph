# First stage of an MCGLM, unit by unit (C++ kernel)

Forms, for every unit of a partition of the rows, the sum over its rows
of the quasi-score and of the Pearson estimating functions of the first
stage; the sensitivity of the first stage from its definition, the block
of the covariance parameters on the regression parameters included; and
the derivative of the Pearson residuals in the first-stage parameters.

## Usage

``` r
cpp_first_stage_units(invC, C, Dblocks, DC, DCbeta, res, w, unit, n_units)
```

## Arguments

- invC:

  the \\nR \times nR\\ inverse of the first-stage covariance.

- C:

  the first-stage covariance.

- Dblocks:

  list of \\R\\ dense \\n \times p_r\\ derivatives of the means in the
  regression parameters of each response.

- DC:

  list of the derivatives of \\C\\ in the covariance parameters.

- DCbeta:

  list, one entry per regression parameter, of the derivatives of \\C\\
  in that parameter, or `NULL` where the variance does not depend on the
  mean.

- res:

  raw residuals \\y - \mu\\, stacked response by response.

- w:

  prior weights of the first stage, one per stacked row.

- unit:

  zero-based unit of each observational row, length \\n\\.

- n_units:

  number of units.

## Value

A list with `psi` (`n_units` by \\p_1\\), the sensitivity `S` (\\p_1
\times p_1\\, on the sign convention of mcglm,
\\E(\partial\psi/\partial\alpha^\top)\\), `G` (\\nR\\ by \\p_1\\, the
derivative of the Pearson residuals with the sign removed), `r` (the
Pearson residuals), and the numbers `n_beta` of regression parameters
and `n_cov` of covariance parameters, with \\p_1\\ their sum.
