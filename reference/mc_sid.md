# Structural intervention distance

Counts the ordered pairs of responses for which the estimated graph
would send one to the wrong adjustment set (Peters and Buehlmann, 2015).
For the pair \\(i, j)\\ the estimate prescribes adjusting for the
parents of \\i\\ that it carries; the pair is counted when that set is
not a valid adjustment set for the effect of \\i\\ on \\j\\ in the true
graph. The distance is zero when the estimate gives a valid adjustment
set for every pair, and it is not symmetric in its two arguments.

## Usage

``` r
mc_sid(truth, estimated, max_members = 1000L)

# S3 method for class 'mcg_sid'
print(x, ...)
```

## Arguments

- truth:

  the true DAG, as a square 0/1 adjacency matrix where a one in position
  \\(a, b)\\ means \\a \to b\\. It must be acyclic and have no
  undirected edge.

- estimated:

  the estimated graph with the same number of nodes and the same
  convention: a DAG, a CPDAG (undirected edges carry a one in both
  positions), or an object of class `"mcg_cpdag"`.

- max_members:

  cap on the number of members enumerated for a partially directed
  estimate.

- x:

  an object of class `"mcg_sid"`.

- ...:

  ignored.

## Value

For a directed estimate, a single number, the count of offending pairs,
with the \\R \times R\\ 0/1 matrix of the offending ordered pairs
attached as the attribute `"pairs"` (a one in \\(i, j)\\ marks the pair
\\(i, j)\\). For a partially directed one, an object of class
`"mcg_sid"`, a list with `lower` and `upper`, the smallest and largest
distance over the members; `values`, the distance of each member;
`n_members`; `which_lower` and `which_upper`, the positions in the
member list (as returned by
[`mc_class_members`](https://evandeilton.github.io/mcglmgraph/reference/mc_class_members.md))
attaining each end; and `truncated`. The `print` method returns `x`
invisibly.

## Why it is not the Hamming distance

The structural Hamming distance counts the edges two graphs disagree
about and says nothing about what the disagreement costs. Two estimates
at the same Hamming distance from the truth can be worth very different
things, one still giving the right adjustment set for every pair and the
other for none. This distance asks that question instead.

## Which adjustment sets count as valid

Validity is decided by the adjustment criterion of Shpitser, VanderWeele
and Robins (2010), in the form given by Perkovic et al. (2018), which is
necessary and sufficient, and not by the back-door criterion of Pearl
(2009). The back-door criterion is sufficient without being necessary:
it refuses any set containing a descendant of \\i\\, while a descendant
that lies on no causal path from \\i\\ to \\j\\ does no harm, so using
it would count pairs the adjustment in fact recovers. Writing
\\\mathrm{cn}(i,j)\\ for the nodes other than \\i\\ on the directed
paths from \\i\\ to \\j\\ in the true graph, and \\\mathrm{forb}(i,j)\\
for their descendants (each node counted among its own) together with
\\i\\, the parent set \\S\\ is valid when it does not meet
\\\mathrm{forb}(i,j)\\ and when it d-separates \\i\\ from \\j\\ in the
graph with the arcs from \\i\\ into \\\mathrm{cn}(i,j)\\ removed. With
no causal path from \\i\\ to \\j\\ the second condition asks for
separation in the whole graph, which is what recovers a null effect.

One case needs care. When \\j\\ itself belongs to the adjustment set the
formula collapses to the marginal of \\j\\, since conditioning on it and
averaging over it returns it unchanged. That is right exactly when the
effect is null, and wrong otherwise. Reading it as always wrong would
make the distance of a graph to itself nonzero, because a parent of
\\i\\ always sits in its adjustment set and is never its descendant.

## A partially directed estimate gives an interval

A completed partially directed acyclic graph (CPDAG) does not prescribe
one adjustment set, so the distance is not one number but the range over
the members of the class. The members are enumerated with
[`mc_class_members`](https://evandeilton.github.io/mcglmgraph/reference/mc_class_members.md)
and the result reports the lower and upper end, the member attaining
each, and whether the enumeration had to be capped. When it was capped,
the interval covers the enumerated members only.

## What the distance measures

The reading in terms of interventions, in which a valid adjustment set
recovers the effect of one response on another, holds when the true
graph is a causal DAG: a structural model with independent exogenous
innovations, no unobserved confounding and no selection bias, whose
equations are invariant under intervention. Without those assumptions
the distance remains a purely graphical measure of disagreement between
the two structures, computed from the parent sets of the estimate. It
compares graphs and uses no data.

## References

Peters, J. and Buehlmann, P. (2015). Structural intervention distance
for evaluating causal graphs. *Neural Computation* 27(3), 771-799.
[doi:10.1162/NECO_a_00708](https://doi.org/10.1162/NECO_a_00708)

Shpitser, I., VanderWeele, T. and Robins, J. M. (2010). On the validity
of covariate adjustment for estimating causal effects. *Proceedings of
the 26th Conference on Uncertainty in Artificial Intelligence*, 527-536.

Perkovic, E., Textor, J., Kalisch, M. and Maathuis, M. H. (2018).
Complete graphical characterization and construction of adjustment sets
in Markov equivalence classes of ancestral graphs. *Journal of Machine
Learning Research* 18(220), 1-62.

Pearl, J. (2009). *Causality: Models, Reasoning and Inference*, 2nd
edition. Cambridge University Press.

## See also

[`mc_cpdag_metrics`](https://evandeilton.github.io/mcglmgraph/reference/mc_cpdag_metrics.md)
for the Hamming-type measures,
[`mc_class_members`](https://evandeilton.github.io/mcglmgraph/reference/mc_class_members.md)
for the members of a class.

Other equivalence class:
[`mc_class_coefficients()`](https://evandeilton.github.io/mcglmgraph/reference/mc_class_coefficients.md),
[`mc_class_members()`](https://evandeilton.github.io/mcglmgraph/reference/mc_class_members.md)

## Examples

``` r
# truth 1 -> 2 -> 3; an estimate that reverses the first arc
G <- matrix(0L, 3, 3); G[1, 2] <- G[2, 3] <- 1L
H <- matrix(0L, 3, 3); H[2, 1] <- H[2, 3] <- 1L
mc_sid(G, H)
#> [1] 3
#> attr(,"pairs")
#>      [,1] [,2] [,3]
#> [1,]    0    1    1
#> [2,]    1    0    0
#> [3,]    0    0    0
mc_sid(G, G)          # zero against itself
#> [1] 0
#> attr(,"pairs")
#>      [,1] [,2] [,3]
#> [1,]    0    0    0
#> [2,]    0    0    0
#> [3,]    0    0    0
```
