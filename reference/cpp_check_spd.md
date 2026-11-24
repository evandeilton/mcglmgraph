# Test symmetric positive-definiteness via the smallest eigenvalue (C++ kernel)

The symmetric part of K is declared positive definite when its smallest
eigenvalue exceeds tol.

## Usage

``` r
cpp_check_spd(K, tol = 1e-10)
```

## Arguments

- K:

  square matrix; only its symmetric part is used.

- tol:

  numeric tolerance on the smallest eigenvalue (default 1e-10).

## Value

list with is_spd (logical) and min_eig (the smallest eigenvalue).
