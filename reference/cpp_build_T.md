# Build the unit lower-triangular structural matrix T (C++ kernel)

Places -gamma_e at position (j, i) of the identity for each edge e = (i,
j), i \< j, so that T is unit lower triangular with T\[j, i\] =
-gamma_ji.

## Usage

``` r
cpp_build_T(gamma, R, edge_idx)
```

## Arguments

- gamma:

  numeric vector of structural coefficients, one per edge, in the row
  order of edge_idx.

- R:

  integer number of responses.

- edge_idx:

  integer matrix \|E\| x 2 with edges (i, j), i \< j, 1-based.

## Value

the R x R unit lower-triangular matrix T with T(j,i) = -gamma.
