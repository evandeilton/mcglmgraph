# Totals of the first-stage estimating functions (C++ kernel)

The quasi-score and the Pearson functions of an MCGLM summed over all
rows, \\\mathbf D^\top\mathbf C^{-1}\mathbf W(\mathbf
y-\boldsymbol\mu)\\ and \\(\mathbf y-\boldsymbol\mu)^\top\mathbf
C^{-1}\dot{\mathbf C}\_l\mathbf W \mathbf C^{-1}(\mathbf
y-\boldsymbol\mu)-\operatorname{tr}(\mathbf C^{-1} \dot{\mathbf
C}\_l\mathbf W)\\, without the sensitivity: what a numerical derivative
of the first stage needs at every step.

## Usage

``` r
cpp_first_stage_totals(invC, Dblocks, DC, res, w)
```

## Arguments

- invC:

  the \\nR \times nR\\ inverse of the first-stage covariance.

- Dblocks:

  list of \\R\\ dense \\n \times p_r\\ derivatives of the means in the
  regression parameters of each response.

- DC:

  list of the derivatives of \\C\\ in the covariance parameters.

- res:

  raw residuals \\y - \mu\\, stacked response by response.

- w:

  prior weights of the first stage, one per stacked row.

## Value

The vector of the \\p_1\\ totals.
