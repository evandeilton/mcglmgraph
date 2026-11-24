# Confusion counts of two edge indicators (C++ kernel)

Confusion counts of two edge indicators (C++ kernel)

## Usage

``` r
cpp_confusion(est, tru)
```

## Arguments

- est:

  logical vector of estimated edges.

- tru:

  logical vector of true edges, same length.

## Value

integer vector (tp, fp, fn, tn) over the positions where both are not
missing.
