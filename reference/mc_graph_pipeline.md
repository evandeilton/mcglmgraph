# Constraint-based CPDAG discovery

Two-phase structure learning, self-contained in mcglmgraph, that
estimates the completed partially directed acyclic graph (CPDAG) of the
Markov equivalence class of the directed acyclic graph (DAG) underlying
the second-order dependence between the responses. Phase one recovers
the skeleton and the separating sets from the between-response
covariance via
[`mc_pc_skeleton`](https://evandeilton.github.io/mcglmgraph/reference/mc_pc_skeleton.md).
Phase two orients the edges up to the Markov equivalence class: the
native route (`orient_method = "pc"`) orients v-structures by the
separating-set rule and then applies Meek's rules; the optional route
(`"hc"` or `"tabu"`) delegates the second phase to bnlearn's score-based
search. A blacklist prevents arcs outside the native skeleton, but the
score search may remove native skeleton edges; this route is therefore
not a pure orientation of a fixed skeleton.

## Usage

``` r
mc_graph_pipeline(
  fit,
  alpha = 0.05,
  orient_method = c("pc", "hc", "tabu"),
  orient_rule = c("standard", "conservative", "majority"),
  required = NULL,
  forbidden = NULL,
  max_ord = NULL,
  df_residual = NULL
)

# S3 method for class 'mcg_cpdag'
print(x, ...)
```

## Arguments

- fit:

  an object of class `"mcglm_graph"`.

- alpha:

  significance level of the second-order separation test.

- orient_method:

  `"pc"` (native, default), `"hc"` or `"tabu"` (bnlearn). When bnlearn
  is not installed, the last two fall back on `"pc"` with a warning, and
  the result records it in `orientation_fallback`.

- orient_rule:

  rule for the unshielded triples of the native route: `"standard"`
  (default), `"conservative"` or `"majority"`. It has no effect when
  `orient_method` is a score-based route.

- required, forbidden:

  optional background knowledge, as two-column matrices of arcs \\(from,
  to)\\ given either by response name or by index. An arc in `required`
  must be present with that direction, an arc in `forbidden` must not.
  They are imposed before propagation, and Meek's fourth rule, idle
  without them, carries their consequences through the rest of the
  graph. Available on the native route only.

- max_ord:

  maximum conditioning order; `NULL` uses \\R - 2\\.

- df_residual:

  optional residual degrees of freedom \\\nu_e\\: `NULL` (default) uses
  the raw-sample contract based on \\n\\; a number or the result of
  [`mc_residual_df`](https://evandeilton.github.io/mcglmgraph/reference/mc_residual_df.md)
  is used as given; `"auto"` computes it with
  [`mc_residual_df`](https://evandeilton.github.io/mcglmgraph/reference/mc_residual_df.md)
  and refuses a fit for which it is not exact.

- x:

  an object of class `"mcg_cpdag"`.

- ...:

  ignored.

## Value

An object of class `"mcg_cpdag"`, a list with the 0/1 matrix `skeleton`
of the oriented structure and its `edges` (two-column integer matrix);
the constraint-stage `pc_skeleton` and `sepset`, as returned by
[`mc_pc_skeleton`](https://evandeilton.github.io/mcglmgraph/reference/mc_pc_skeleton.md);
the `cpdag` adjacency in the convention of
[`mc_orient_meek`](https://evandeilton.github.io/mcglmgraph/reference/mc_orient_meek.md),
labeled by response; the `vstructures` table (`NULL` on a score-based
route); the `conflicts` table of edges required in opposite directions;
the named logical vector `gates`; the logical `valid_cpdag`; the `label`
to be used when reporting the output; `knowledge` and `n_imposed`, the
background knowledge with the checks on it and the number of arcs it
oriented; `orient_method`, `orient_rule`, `orient_method_requested` and
`orientation_fallback`; `sepsets`, the collections enumerated under a
rule other than `"standard"` (`NULL` otherwise); `ambiguous`, the
triples left undecided, and their number `n_ambiguous`;
`meek_iterations` and `meek_blocked`, the number of sweeps of Meek's
rules and of orientations withheld by an ambiguous triple;
`separation_diagnostics`, the diagnostics of the skeleton search; and
`alpha`, `df_residual` and `response_names`. The `print` method returns
`x` invisibly.

## Which separating sets the orientation consults

`orient_rule` selects between the rule that reads the single set the
skeleton recorded (`"standard"`) and the two rules that enumerate every
set separating the endpoints of each unshielded triple on the final
skeleton, `"conservative"` and `"majority"`. The enumeration is native
and uses the same conditional-association test and the same
degrees-of-freedom contract as the skeleton, so a difference between two
runs is attributable to the rule and not to a change of test. A triple
that the rule leaves ambiguous receives no orientation, and Meek's R1,
R3 and R4 are withheld where they would settle it by propagation; the
count is reported in `n_ambiguous`.

## What background knowledge changes

Without it the output is the CPDAG of an equivalence class, the object
the second-order separations identify. With it the output is no longer
the CPDAG of any class: when no triple is left ambiguous it is the
maximally oriented partially directed graph of Meek (1995), the union of
the members of the class that satisfy the knowledge. The `label` records
which of the two it is, and `knowledge` carries the arcs that were
imposed. Knowledge is refused, rather than silently dropped, when it
contradicts itself, when it requires an arc between responses the
skeleton does not join, and when it contradicts an edge oriented by a
v-structure. A contradiction with an edge that only Meek's rules would
orient is not refused beforehand: it surfaces during propagation as a
conflict or a cycle and the output then fails the checks below.

## What the output may be called

The output is submitted to seven checks, reported in `gates`: no pair of
colliders or of rules required opposite directions on an edge
(`no_orientation_conflict`); no separation test needed a ridge to invert
its block (`unregularized_separation`); the orientation stayed within
the learned skeleton (`orientation_within_skeleton`) and kept all of it
(`skeleton_preserved`); the directed part has no cycle (`acyclic`);
every collider that was written survives Meek propagation
(`vstructures_preserved`); and a consistent DAG extension exists, by the
algorithm of Dor and Tarsi (1992) (`consistent_extension`). When all
pass, `valid_cpdag` is `TRUE` and `label` is `"CPDAG"`, or
`"maximally oriented graph under background knowledge"` when the
knowledge oriented some edge. When any fails, `valid_cpdag` is `FALSE`
and `label` becomes `"second-order partially directed structure"`: the
object still records what the separations supported, but it is not
certified as the CPDAG of a Markov equivalence class. The first two
checks record facts of the run and not properties of the returned graph,
so a refused label does not establish that the output differs from the
CPDAG of a class. A triple left ambiguous by the conservative or
majority rule stays unoriented, so a certified output can be less
oriented than the CPDAG of the class of its extension.

## Assumptions behind the reading

Passing the checks is a statement about internal consistency only.
Reading the output as the CPDAG of the second-order structure assumes
that the between-response covariance is generated by some DAG and that
every zero partial correlation corresponds to a d-separation in it
(second-order Markov and faithfulness). Reading it as the Markov
equivalence class of the distribution of the responses assumes, in
addition, that a zero partial correlation of the residuals is equivalent
to conditional independence of the responses given the conditioning
responses and the covariates, which holds for jointly Gaussian residuals
whose conditional covariance does not depend on the covariates and fails
in general otherwise. Recovery of the class is asymptotic, for a fixed
number of responses and a level `alpha` that decreases with the sample
size. A causal reading of the arcs requires further causal sufficiency,
absence of selection bias and substantive justification. The procedure
checks none of these assumptions.

## References

Spirtes, P., Glymour, C. and Scheines, R. (2000). *Causation,
Prediction, and Search*, 2nd edition. MIT Press.

Meek, C. (1995). Causal inference and causal explanation with background
knowledge. *Proceedings of the 11th Conference on Uncertainty in
Artificial Intelligence*, 403-410.

Dor, D. and Tarsi, M. (1992). A simple algorithm to construct a
consistent extension of a partially oriented graph. Technical Report
R-185, Cognitive Systems Laboratory, University of California, Los
Angeles.

Andersson, S. A., Madigan, D. and Perlman, M. D. (1997). A
characterization of Markov equivalence classes for acyclic digraphs.
*The Annals of Statistics* 25(2), 505-541.
[doi:10.1214/aos/1031833662](https://doi.org/10.1214/aos/1031833662)

## See also

[`mc_pc_skeleton`](https://evandeilton.github.io/mcglmgraph/reference/mc_pc_skeleton.md),
[`mc_separating_sets`](https://evandeilton.github.io/mcglmgraph/reference/mc_separating_sets.md),
[`mc_orient_meek`](https://evandeilton.github.io/mcglmgraph/reference/mc_orient_meek.md),
[`mc_class_members`](https://evandeilton.github.io/mcglmgraph/reference/mc_class_members.md),
[`mc_ges`](https://evandeilton.github.io/mcglmgraph/reference/mc_ges.md)
for the score-based alternative.

Other structure discovery:
[`mc_find_vstructures()`](https://evandeilton.github.io/mcglmgraph/reference/mc_find_vstructures.md),
[`mc_ges()`](https://evandeilton.github.io/mcglmgraph/reference/mc_ges.md),
[`mc_orient_meek()`](https://evandeilton.github.io/mcglmgraph/reference/mc_orient_meek.md),
[`mc_pc_skeleton()`](https://evandeilton.github.io/mcglmgraph/reference/mc_pc_skeleton.md),
[`mc_residual_df()`](https://evandeilton.github.io/mcglmgraph/reference/mc_residual_df.md),
[`mc_separating_sets()`](https://evandeilton.github.io/mcglmgraph/reference/mc_separating_sets.md),
[`mc_skeleton_to_bnlearn()`](https://evandeilton.github.io/mcglmgraph/reference/mc_skeleton_to_bnlearn.md)

## Examples

``` r
# \donttest{
set.seed(1); n <- 600
z1 <- rnorm(n); z2 <- rnorm(n); z3 <- 0.7*z1 + 0.7*z2 + rnorm(n)
d <- data.frame(y1 = z1, y2 = z2, y3 = z3, x = rnorm(n))
f <- mcglm_graph(list(y1~x, y2~x, y3~x), data = d)
#> Automatic initial values selected.
mc_graph_pipeline(f)$cpdag
#>    y1 y2 y3
#> y1  0  0  1
#> y2  0  0  1
#> y3  0  0  0
# }
```
