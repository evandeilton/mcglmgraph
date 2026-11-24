# Is the two-stage correction available for a fit

Reports which of the two routes that propagate the first-stage
uncertainty into the graph layer apply to a fit, and if one does not,
why. The shortcut is the reduced form evaluated at the fitted covariance
(see the section "Propagating the first stage" of
[`mc_vcov_graph`](https://evandeilton.github.io/mcglmgraph/reference/mc_vcov_graph.md)).
The reduced form needs every response to carry a single dispersion
parameter that multiplies its variance function,
\\v\_{ri}=\tau\_{r}V_r(\mu\_{ri};p_r)\\, an identity within-response
structure and constant prior weights, because together they make the
first-stage Pearson equations fix the average of the squared residuals
of each response and let the correction be written, to first order, with
second-stage quantities alone. It also needs the first stage to estimate
the correlation between responses without restriction, which mcglm does
and the check therefore does not read. The check below is that
structural one: on a DAG with an immorality the shortcut can qualify and
still be, for the parameters of the collider and of its descendants,
only a first-order approximation to the reduced form, which is one of
the two cases in which `method = "auto"` takes the full form when it is
available. The other is a first stage with a variance function other
than `"constant"`: the shortcut then rests on conditions on the mean and
on the residual covariance that the fit does not show (see the section
"What the check cannot verify"). The full form needs only that the first
stage be a sum over units, which is what the rows are when the
within-response structure is the identity, and what the clusters must be
otherwise.

## Usage

``` r
mc_two_stage_available(fit, cluster = NULL)
```

## Arguments

- fit:

  an object of class `"mcglm_graph"`.

- cluster:

  optional vector with one identifier per observational row, as in
  [`mc_vcov_graph`](https://evandeilton.github.io/mcglmgraph/reference/mc_vcov_graph.md);
  default `NULL`. It changes only whether the full form is available for
  a within-response structure that is not the identity, and with it the
  route `method = "auto"` takes on a DAG with an immorality or on a fit
  with a variance function other than `"constant"`.

## Value

A list with `shortcut_ok` (logical) and `why_shortcut` (character
vector, empty when the shortcut applies), the verdict and the reasons
for the shortcut; `full_ok` and `why_full`, the same for the full form;
and `method`, a string, `"shortcut"` or `"full"`, the route
`method = "auto"` takes in `mc_vcov_graph(type = "two_stage")` and in
the functions that consume it, or `NA_character_` when neither applies.
The shortcut can apply (`shortcut_ok = TRUE`) while `method` is
`"full"`: `"auto"` takes the full form, when it is available, for a fit
with a variance function other than `"constant"` and for a DAG with an
immorality. The components `ok` and `why` repeat `shortcut_ok` and
`why_shortcut`.

## Which fits qualify for the shortcut

The variance functions `"constant"`, `"tweedie"`, `"binomialP"` and
`"binomialPQ"` are multiplicative in the dispersion parameter and
qualify. `"poisson_tweedie"` is additive, \\v=\mu+\tau\mu^{p}\\, and
qualifies only at \\p=1\\, where it collapses to \\(1+\tau)\mu\\; the
power must be fixed there rather than estimated, since an estimated
power is a first-stage parameter of its own. A response with more than
one covariance component has more than one dispersion parameter and does
not qualify. Prior weights must be constant: with weights that vary, the
first stage fixes a weighted average of the squared residuals rather
than their average.

A single dispersion parameter is necessary and not sufficient. The
identity comes from the first-stage equations being the score of a
working Gaussian model that factorizes over units, which needs the
within-response structure to be a multiple of the identity. One
covariance component that is not diagonal keeps the parameter count at
one and breaks the factorization all the same, so it is refused too.

## Which fits qualify for the full form

Any family, any number of covariance components and any power, estimated
or fixed, provided the fit carries its first-stage components and the
graph layer sits on its Pearson residuals. With an identity
within-response structure the units are the rows, or the clusters when
`cluster` is given. With any other structure `cluster` is required and
has to nest it: every matrix of the matrix linear predictor must vanish
between rows of distinct clusters, which is checked. A random-effect
block per subject nested in a cluster per hunter, say, qualifies with
the hunter as the cluster.

## What the check cannot verify

The check reads the structure of the first-stage fit and nothing else:
the variance functions, the number of dispersion parameters, whether the
within-response structure is the identity or nests in the clusters,
whether a Poisson-Tweedie power was fixed at one and whether the prior
weights are constant. For the shortcut, that structure makes the first
stage fix the average of the squared residuals, which is the algebraic
half of it. The statistical half, that the first stage moves the
residual covariance along a diagonal rescaling, rests on statements
about the data that cannot be read off the fit:

- a correctly specified mean, condition (b) of the section "What the
  correction reaches" of
  [`mc_vcov_graph`](https://evandeilton.github.io/mcglmgraph/reference/mc_vcov_graph.md).
  It is not needed when all responses have the same design, with an
  intercept, and the derivative of the mean in the linear predictor
  divided by the standard deviation is the same for every unit, as with
  the constant variance function and the identity link, or with the
  Tweedie variance function of power fixed at two and the log link;

- a conditional residual covariance that does not depend on the
  covariates, condition (c) there. It is not needed, whatever the link,
  when every response has the constant variance function; a variance
  function that depends on the mean needs it at least in the weaker form
  of conditional second moments of the Pearson residuals that change
  with the covariates only through a factor common to the responses;

- independence between the units, that is between the rows once the
  within-response structure is the identity, or between the clusters.

The first two concern the shortcut only; the full form, with its
observed sensitivity, does not need them. The third concerns every form.
A fit can pass the check and still violate them. The spatial example is
the plain case: a soil survey fitted with a trend surface in the mean
and an identity within-response structure qualifies, while its
neighboring sites remain dependent, since coordinates in the mean do not
whiten the residuals. There the rows are not independent units, and
neither the corrected form nor the sandwich it corrects has the variance
it claims.

Since the first two cannot be verified, and the constant variance
function is the case that does not need the second, `method = "auto"`
takes the full form, when it is available, for a fit in which some
response has another variance function.

## References

Murphy, K. M. and Topel, R. H. (1985). Estimation and inference in
two-step econometric models. *Journal of Business and Economic
Statistics* 3(4), 370–379.
[doi:10.1080/07350015.1985.10509471](https://doi.org/10.1080/07350015.1985.10509471)

Newey, W. K. and McFadden, D. (1994). Large sample estimation and
hypothesis testing. In Engle, R. F. and McFadden, D. L. (eds), *Handbook
of Econometrics*, vol. 4, 2111–2245. Elsevier.
[doi:10.1016/S1573-4412(05)80005-4](https://doi.org/10.1016/S1573-4412%2805%2980005-4)

## See also

[`mc_vcov_graph`](https://evandeilton.github.io/mcglmgraph/reference/mc_vcov_graph.md),
whose `type = "two_stage"` consumes this check, and
[`mc_graph_bootstrap`](https://evandeilton.github.io/mcglmgraph/reference/mc_graph_bootstrap.md),
which propagates the first stage by refitting.

Other inference on the graph layer:
[`mc_test_adequacy()`](https://evandeilton.github.io/mcglmgraph/reference/mc_test_adequacy.md),
[`mc_test_independence()`](https://evandeilton.github.io/mcglmgraph/reference/mc_test_independence.md),
[`mc_vcov_graph()`](https://evandeilton.github.io/mcglmgraph/reference/mc_vcov_graph.md),
[`mc_wald_graph()`](https://evandeilton.github.io/mcglmgraph/reference/mc_wald_graph.md)

## Examples

``` r
# \donttest{
set.seed(1); n <- 200; x <- rnorm(n)
d <- data.frame(y1 = rnorm(n), y2 = rnorm(n), y3 = rnorm(n), x = x)
f <- mcglm_graph(list(y1 ~ x, y2 ~ x, y3 ~ x), data = d)
#> Automatic initial values selected.
mc_two_stage_available(f)
#> $ok
#> [1] TRUE
#> 
#> $why
#> character(0)
#> 
#> $shortcut_ok
#> [1] TRUE
#> 
#> $why_shortcut
#> character(0)
#> 
#> $full_ok
#> [1] TRUE
#> 
#> $why_full
#> character(0)
#> 
#> $method
#> [1] "shortcut"
#> 
# }
```
