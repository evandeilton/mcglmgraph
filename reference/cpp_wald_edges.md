# Per-edge Wald statistics (C++ kernel)

Per-edge Wald statistics (C++ kernel)

## Usage

``` r
cpp_wald_edges(eta, Vcov, edge_pos)
```

## Arguments

- eta:

  full parameter vector.

- Vcov:

  estimated covariance matrix of eta (model-based, J^-1).

- edge_pos:

  0-based positions of edge parameters in eta.

## Value

vector of Wald statistics eta_e^2 / Vcov_ee; NA where the variance is
not finite and positive, since the statistic is then undefined.
