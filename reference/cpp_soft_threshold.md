# Soft-thresholding operator (C++ kernel)

Soft-thresholding operator (C++ kernel)

## Usage

``` r
cpp_soft_threshold(x, thr)
```

## Arguments

- x:

  numeric vector.

- thr:

  non-negative threshold.

## Value

sign(x) \* max(\|x\| - thr, 0).
