# Build the between-response precision matrix (C++ kernel)

Computes K_b = sum_d eta_d Z_d and symmetrizes the result.

## Usage

``` r
cpp_build_precision(eta, bases)
```

## Arguments

- eta:

  numeric vector of length q_G, the number of graph parameters (diagonal
  then edge parameters).

- bases:

  numeric cube R x R x q_G of symmetric base matrices Z_d.

## Value

The R x R symmetric precision matrix K_b = sum_d eta_d Z_d.
