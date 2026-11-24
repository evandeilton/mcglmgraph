# Second-order conditional-association test of one pair (C++ kernel)

The partial correlation of (i, j) given S, read off the inverse of the
(i, j, S) block, with the p-value of Fisher's z, z = atanh(pcor)
sqrt(df), where df = n - \|S\| - 3 or df_residual - \|S\| - 2. The
partial correlation is clamped to \[-0.999999, 0.999999\], and a nearly
singular block is inverted after a ridge, which is reported. It is the
test the skeleton and the enumeration of separating sets both apply.

## Usage

``` r
cpp_ci_pcor(Sigma, i, j, S, n, df_residual = 0)
```

## Arguments

- Sigma:

  R x R covariance matrix of the responses.

- i, j:

  the pair, 1-based.

- S:

  conditioning set, 1-based, possibly empty.

- n:

  number of experimental units.

- df_residual:

  residual degrees of freedom of the first stage; a value not greater
  than zero selects the raw contract n - \|S\| - 3.

## Value

list with pcor, p_value, df, ridge, regularized and usable, the last one
false when the block gives no valid test.
