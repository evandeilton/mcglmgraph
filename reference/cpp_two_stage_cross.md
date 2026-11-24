# Cross sensitivity between the graph layer and the first stage (C++ kernel)

Forms \\\mathbf B_n(d,\cdot)=\sum_i\mathbf
e_i^\top\boldsymbol\Sigma_b^{-1}
\dot{\boldsymbol\Sigma}\_d\boldsymbol\Sigma_b^{-1}\mathbf g_i\\, the
accumulated cross sensitivity \\-\partial\mathbf U_2/\partial
\boldsymbol\alpha^\top\\, from the centered residuals and the derivative
of the Pearson residuals returned by `cpp_first_stage_units`.

## Usage

``` r
cpp_two_stage_cross(Ec, Sinv, Dlist, G)
```

## Arguments

- Ec:

  \\n \times R\\ centered Pearson residuals.

- Sinv:

  the inverse of the fitted \\\boldsymbol\Sigma_b\\.

- Dlist:

  list of the \\P\\ derivatives of \\\boldsymbol\Sigma_b\\.

- G:

  \\nR \times p_1\\ derivative of the Pearson residuals, stacked
  response by response.

## Value

The \\P \times p_1\\ matrix \\\mathbf B_n\\.
