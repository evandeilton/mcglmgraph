# Apply Meek's four orientation rules (C++ kernel)

Apply Meek's four orientation rules (C++ kernel)

## Usage

``` r
cpp_meek_rules(amat, amb_keys, max_iter = 1000L)
```

## Arguments

- amat:

  integer R x R adjacency in the package convention: amat(a,b) and
  amat(b,a) both one is an undirected edge, only amat(a,b) one is a
  -\> b. Each of the four rules is applied in sweeps on the graph as it
  stood at the start of the sweep, and an edge two rules orient in
  opposite directions is left undirected and flagged in conflicts.

- amb_keys:

  keys a\*R\*R + b\*R + c of the ambiguous triples (a, b, c with center
  b), 0-based node indices, both endpoint orders present. Rules R1, R3
  and R4 are withheld where they would settle such a triple; R2 is not.

- max_iter:

  cap on the number of sweeps.

## Value

list with amat, conflicts (R x R), iterations and blocked.
