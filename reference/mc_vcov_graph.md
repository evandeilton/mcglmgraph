# Covariance matrix of the graph-layer parameters

Returns the estimated covariance matrix of the graph-layer parameter
vector \\\hat{\boldsymbol\theta}\\ of a fitted `mcglm_graph` object in
the three forms the package supports: model-based, empirical sandwich,
and sandwich corrected for the estimation of the first stage.

## Usage

``` r
mc_vcov_graph(
  fit,
  type = c("model", "sandwich", "two_stage"),
  cluster = NULL,
  small_sample = FALSE,
  method = c("auto", "shortcut", "full"),
  sensitivity = c("observed", "expected")
)
```

## Arguments

- fit:

  an object of class `"mcglm_graph"`.

- type:

  `"model"` for \\\mathbf J_n^{-1}\\, `"sandwich"` for the empirical
  sandwich, or `"two_stage"` for the sandwich corrected for the
  estimation of the residuals.

- cluster:

  optional vector with one identifier per observational row. With
  `type = "sandwich"` or `"two_stage"`, the contributions are summed
  within clusters before the variability is formed, which is the form to
  use when the unit of dependence is not the row of the data matrix. It
  must have no missing value and at least two distinct values; the
  default `NULL` takes the rows as the units. It is refused with
  `type = "model"`.

- small_sample:

  logical, default `FALSE`; multiply the clustered variability by
  \\G/(G-1)\\, with \\G\\ the number of clusters. It has no effect
  without `cluster` and is refused with `type = "model"`.

- method:

  for `type = "two_stage"`, the route of the correction: `"shortcut"`,
  `"full"`, or `"auto"` (default). `"auto"` takes the full form when it
  is available and either a response has a variance function other than
  `"constant"` or the fitted DAG has an immorality; in every other case
  it takes the shortcut when it applies and the full form when it does
  not. The full form rebuilds the first stage and is slower. It has no
  meaning for the other two types.

- sensitivity:

  for the full form, `"observed"` (default) or `"expected"`, the
  sensitivity of the first stage described in the section "The full
  form". The shortcut and the other two types do not use it:
  `"expected"` is refused with them, and gives a warning when
  `method = "auto"` takes the shortcut.

## Value

A \\q\_{\mathcal G} \times q\_{\mathcal G}\\ symmetric matrix, with
\\q\_{\mathcal G}\\ the number of graph parameters and their names in
its dimnames, carrying the attributes `"type"`, `"variability"`,
`"edge_pos"` (positions of the edge parameters), `"n_clusters"`,
`"small_sample"`, which is `TRUE` exactly when the factor \\G/(G-1)\\
was applied, and the inversion diagnostics `"inversion_status"` and
`"inversion_ridge"` of the sensitivity solve (with `type = "model"` also
`"inversion_rcond"`, the reciprocal condition number met by that solve).
`"n_clusters"` is `NA` without `cluster`. With `type = "two_stage"` it
also carries `"two_stage_method"`, the route taken, and either the
factor `"Q"` the shortcut applies, \\\[\mathbf K_b\dot{\boldsymbol
\Sigma}\_{b,d}\]\_{rr}/(2m_r)\\, with the averages \\m_r\\ in
`"center"`, or, for the full form, the cross sensitivity `"B"`, the
first-stage sensitivity `"A"`, both on the accumulated scale, which
`"sensitivity"` was used and its reciprocal condition number
`"sensitivity_rcond"`. The function refuses fits whose numerical
diagnostics mark graph-layer inference as uncertified.

## Details

The graph parameters solve second-order (Pearson) estimating functions
\\\mathbf
U_n(\boldsymbol\theta)=\sum_i\boldsymbol\psi_i(\boldsymbol\theta)\\
built on the Pearson residuals of the first-stage MCGLM (Bonat and
Jorgensen, 2016). Write \\\mathbf
J=-E(\partial\boldsymbol\psi_i/\partial \boldsymbol\theta^\top)\\ for
the sensitivity and \\\mathbf V=\lim_n
n^{-1}\sum_i\mathrm{Var}(\boldsymbol\psi_i)\\ for the variability, both
per unit. The Godambe information is \\\mathbf G=\mathbf J^\top\mathbf
V^{-1}\mathbf J\\. With the number of responses fixed and the number of
units \\n\to\infty\\, under regularity conditions that include
estimating functions unbiased at \\\boldsymbol\theta_0\\ for every unit
or cluster, with independent contributions over units or clusters, a
root \\\boldsymbol\theta_0\\ that is interior and locally unique, a
nonsingular \\\mathbf J\\, a finite and positive definite \\\mathbf V\\
with a Lindeberg condition, and uniform control of the fourth moments of
the residuals, \$\$\sqrt
n(\hat{\boldsymbol\theta}-\boldsymbol\theta_0)\xrightarrow{d} \mathcal
N\left(\mathbf 0,\\\mathbf G^{-1}\right),\qquad \mathbf G^{-1}=\mathbf
J^{-1}\mathbf V\mathbf J^{-\top}\$\$ (Godambe, 1960; Jorgensen and
Knudsen, 2004). The matrix returned estimates \\n^{-1}\mathbf G^{-1}\\,
the covariance of \\\hat{\boldsymbol\theta}\\ itself. With independent
clusters the same statement holds for the cluster sums of the
contributions, with the number of clusters in the place of \\n\\. The
robustness of the sandwich is to the fourth-order structure around the
root identified by the estimating functions; it does not correct bias
from a misspecified mean or second-order covariance, and it treats the
graph as fixed.

With `type = "model"` the covariance is the model-based \\\mathbf
J_n^{-1}\\, valid under the working identity \\\mathbf V = \mathbf J\\
when the relevant second moment is correctly specified and Gaussian
fourth-moment identities hold. With `type = "sandwich"`, also called the
naive sandwich because it takes the Pearson residuals as observed, it is
\\\mathbf H_n^{-1}\mathbf V_n\mathbf H_n^{-1}\\, with \\\mathbf H_n\\
the curvature of the section "Scale convention" and \\\mathbf
V_n=n\widehat{\mathbf V}\\ the sum of the outer products of the Pearson
contributions, an empirical fourth-moment variability that does not
assume that identity. In the DAG layer its block for the structural
coefficients of a node is the heteroskedasticity-consistent estimator of
White (1980) for the least-squares regression of the centered residual
of that node on those of its parents. The unit contribution to the
estimating function of parameter \\\theta_d\\ is
\$\$\psi\_{d,i}=-\tfrac12\left\\\operatorname{tr}(\boldsymbol\Sigma_b^{-1}
\dot{\boldsymbol\Sigma}\_{b,d})-\mathbf r_i^\top\boldsymbol\Sigma_b^{-1}
\dot{\boldsymbol\Sigma}\_{b,d}\boldsymbol\Sigma_b^{-1}\mathbf
r_i\right\\,\$\$ and \\\widehat{\mathbf
V}=n^{-1}\sum_i\boldsymbol\psi_i\boldsymbol\psi_i^\top\\, where
\\\mathbf r_i\\ is the vector of centered Pearson residuals of unit
\\i\\. This empirical variability requires scores centered in
expectation for each independent unit or cluster, or a vanishing
aggregate of their squared means. A zero total score at the fitted root
does not establish that condition under heterogeneous second-moment
misspecification. With `cluster`, the scores are summed by cluster
before taking outer products; independence of those clusters remains an
assumption about the data.

## Scale convention

The model-based sensitivity stored in `fit$graph$J` is on the
accumulated scale, \\\mathbf J_n=n\mathbf J\\, while the variability
above is per unit. For a sparse DAG, the empirical sandwich and
two-stage correction use the observed negative Hessian \\\mathbf H_n\\
of the working criterion. Its parent blocks use the residual sample
covariance; it can differ from \\\mathbf J_n\\ if the DAG misses
second-moment structure. For the linear-precision GGM they coincide. The
empirical matrix returned here is \\\mathbf H_n^{-1}(n\widehat{\mathbf
V})\mathbf H_n^{-1}\\, which is the covariance of
\\\hat{\boldsymbol\theta}\\ and not of \\\sqrt
n\\\hat{\boldsymbol\theta}\\, so its square roots are directly
comparable with the estimates returned by `coef(fit, type = "graph")`.

## Propagating the first stage

The first two forms refer to the second stage conditionally on the
estimated Pearson residuals: they carry the sampling variability of the
estimating functions of the graph layer and not the uncertainty
inherited from the first-stage fit. `type = "two_stage"` propagates it,
with the number of responses fixed, \\n\to\infty\\ and the regularity
conditions of the Details holding for the stacked system. Stacking the
two sets of estimating functions gives a block triangular sensitivity,
and block inversion gives the influence function of the second stage,
\$\$\mathrm{IF}\_i=\mathbf H^{-1}\left(\boldsymbol\psi_i-\mathbf B
\mathbf A^{-1}\boldsymbol\psi\_{\alpha,i}\right),\qquad \mathbf
B=-E\\\left(\frac{\partial\boldsymbol\psi_i}
{\partial\boldsymbol\alpha^\top}\right).\$\$ Here \\\mathbf
H=-E(\partial\boldsymbol\psi_i/\partial \boldsymbol\theta^\top)\\ is the
curvature of the graph block under the distribution of the data,
estimated by the observed curvature \\\mathbf H_n/n\\, and \\\mathbf
A=-E(\partial
\boldsymbol\psi\_{\alpha,i}/\partial\boldsymbol\alpha^\top)\\ is the
sensitivity of the first stage, whose parameters are
\\\boldsymbol\alpha\\; both are taken nonsingular. The covariance of
\\\hat{\boldsymbol\theta}\\ is \\n^{-1}\mathbf H^{-1}\widetilde{\mathbf
V}\mathbf H^{-\top}\\ with \$\$\widetilde{\mathbf V}=\mathbf
V\_{\theta\theta} -\mathbf B\mathbf A^{-1}\mathbf V\_{\alpha\theta}
-\mathbf V\_{\theta\alpha}\mathbf A^{-\top}\mathbf B^\top +\mathbf
B\mathbf A^{-1}\mathbf V\_{\alpha\alpha}\mathbf A^{-\top} \mathbf
B^\top,\$\$ where \\\mathbf V\_{\alpha\alpha}\\, \\\mathbf
V\_{\alpha\theta}\\ and \\\mathbf V\_{\theta\theta}\\ are the blocks of
the variability of the stacked contributions. The first two types keep
only \\\mathbf V\_{\theta\theta}\\; the difference between the corrected
and the naive covariance has no determined sign. This is the correction
of Murphy and Topel (1985) and, as a stacked system, the two-step result
of Newey and McFadden (1994, Theorem 6.1).

When every response carries a single dispersion parameter that
multiplies its variance function (a Poisson-Tweedie variance only with
its power fixed at one), the within-response structure is the identity
and the prior weights are constant, the first-stage term can be written,
to first order, with second-stage quantities only. This is the *reduced
form*, and these three conditions are those of the shortcut below. It
also uses that the first stage estimates the correlation between
responses without restriction, as mcglm does, and that the first stage
acts on the residual covariance as a diagonal rescaling, for which
conditions (a) to (c) of the section "What the correction reaches" are
sufficient: \$\$\mathbf B\mathbf
A^{-1}\boldsymbol\psi\_{\alpha,i}=\mathbf Q (\mathbf r_i\circ\mathbf
r_i-\mathbf m),\qquad Q\_{dr}=\frac{1}{2m_r}\bigl\[\mathbf
K_b\dot{\boldsymbol\Sigma}\_{b,d} \mathbf K_b\widehat{\mathbf
S}\_e\bigr\]\_{rr},\qquad
\dot{\boldsymbol\Sigma}\_{b,d}=\partial\boldsymbol\Sigma_b/\partial\theta_d,\$\$
with \\\circ\\ the elementwise product, \\\mathbf K_b=
\boldsymbol\Sigma_b^{-1}\\, \\\widehat{\mathbf S}\_e\\ the sample
covariance of the Pearson residuals and \\m_r\\ the average of the
squared centered Pearson residuals of response \\r\\, which, with
constant prior weights and residuals that average zero, is the value at
which the first-stage Pearson equations fix it. Since \\r\_{ri}^2\\
scales as \\1/\tau_r\\, unit \\i\\ moves \\\log\hat\tau_r\\ by
\\(r\_{ri}^2-m_r)/(n m_r)\\, which is where the factor \\1/m_r\\ comes
from. Without the bias correction of mcglm
(`control_algorithm = list(correct = FALSE)`) \\m_r=1\\. With it, its
default, and unit weights, \\m_r=(n-p_r)/n\\ when the responses share
the design (or there is a single response), \\p_r\\ the number of
regression parameters: exactly for the constant variance function with
the identity link and for power two with the log link, and closely for
the other multiplicative variance functions. A constant weight \\w_0\\
turns \\p_r/n\\ into \\p_r/(nw_0^2)\\, and with designs that differ
\\m_r\\ is the observed average, slightly below one. Centering on
\\m_r\\ makes the correction term sum to zero in the sample, and a
parameter the first stage fixes exactly, such as the innovation variance
of a root node of a DAG with Gaussian intercept-only margins, gets a
zero standard error. The *shortcut*, which is what `method = "shortcut"`
computes, is the reduced form evaluated at the fitted covariance,
\\Q\_{dr}=\[\mathbf K_b\dot{\boldsymbol\Sigma}\_{b,d}\]\_{rr}/(2m_r)\\.
The two coincide in every row of \\\mathbf Q\\ whose entries the fitted
and the sample covariances agree on: at every support of the GGM layer,
and for the parameters of every node of a DAG without immoralities. For
the parameters of a collider with non-adjacent parents, and of its
descendants, the shortcut is a first-order approximation, across every
response column of the parent block, the parents included, that also
needs the DAG to be correctly specified; for the parameters of every
other node it equals the reduced form. The *full form* evaluates
\\\mathbf B\\ at the data from the first-stage contributions and needs
none of these conditions.

## When the shortcut is right

The statements of this section are of first order, with the number of
responses fixed, \\n\to\infty\\ and expectations read as limits of
averages over independent units. The reduced form rests on two facts of
different kinds. One is algebraic: under the three conditions of the
previous section, and with the correlation between responses of the
first stage unrestricted, the squared Pearson residual of each response,
minus its average over the units, is a linear combination of the
centered first-stage contributions of the same unit, whatever the data.
That is what
[`mc_two_stage_available`](https://evandeilton.github.io/mcglmgraph/reference/mc_two_stage_available.md)
checks. The other is statistical. Let \\\mathbf
M(\boldsymbol\alpha)=E\\\mathbf r_i(\boldsymbol\alpha) \mathbf
r_i(\boldsymbol\alpha)^\top\\\\ be the second moment of the Pearson
residuals as a function of the first-stage parameters and \\\mathbf
M_0\\ its value at their limit, which has unit diagonal. Let the Pearson
residuals have mean zero at that limit, so that centering them, as the
graph layer does, changes neither \\\mathbf M_0\\ nor its derivatives.
The reduced form holds when every first-stage parameter \\\alpha_m\\
moves \\\mathbf M\\ as a diagonal rescaling, that is, when for some
diagonal matrix \\\boldsymbol\Lambda_m\\ \$\$\frac{\partial\mathbf
M}{\partial\alpha_m}= -(\boldsymbol\Lambda_m\mathbf M_0+\mathbf
M_0\boldsymbol\Lambda_m)\$\$ on the diagonal and in the entries that the
matrices \\\mathbf K_b \dot{\boldsymbol\Sigma}\_{b,d}\mathbf K_b\\ read:
the edges in the GGM layer, the block of each node and its parents in
the DAG layer. Conditions (b) and (c) of the section "What the
correction reaches", the latter with \\\mathbf M_0\\ in the place of
\\\boldsymbol\Sigma_b\\, give it, with
\\\boldsymbol\Lambda_m=\operatorname{Diag}\\E(\partial\log\varsigma\_{ri}/
\partial\alpha_m)\\\\. They are sufficient and can be weakened in two
ways:

- \(c\) can be replaced by \\E(\mathbf r_i\mathbf r_i^\top\mid \mathbf
  x_i)=c(\mathbf x_i)\mathbf M_0\\, with \\c\\ a positive function
  common to the responses, and \\\boldsymbol\Lambda_m\\ then carries
  \\E\\c(\mathbf x_i)\\\partial\log\varsigma\_{ri}/\partial\alpha_m\\\\.
  When every response has the constant variance function (c) is not
  needed in any form, whatever the link, because
  \\\partial\log\varsigma\_{ri}/\partial\alpha_m\\ is then the same for
  every unit;

- \(b\) is used only to give the Pearson residuals mean zero and to make
  \\E(r\_{si}\\\varsigma\_{ri}^{-1}\partial\mu\_{ri}/\partial
  \boldsymbol\beta_r)=\mathbf 0\\ for every pair of responses \\(r,s)\\.
  When all responses have the same design, with an intercept, and the
  variance function and the link make that vector the row of the design
  times a factor that is the same for every unit, as the constant
  variance function with the identity link and the Tweedie variance
  function of power fixed at two with the log link do, the quasi-score
  equations of the first stage make the Pearson residuals orthogonal to
  the design in the sample, and (b) is not needed. In the first of these
  two cases the reduced form then holds with no condition on the mean or
  on the conditional second moments. In the second, (c), or its weaker
  form, is required of the Pearson residuals about the limit of the
  fitted mean; when it holds for those about the true mean, a
  misspecified mean violates it only at second order in the relative
  difference between the two means.

Outside these cases the shortcut can omit a term of first order, whose
sign is not determined. That is so for a mean that misses a covariate
when the designs of the responses differ and, with a variance function
that depends on the mean, for a residual correlation that moves with the
covariates; the variances of the partial correlations can then be
affected as well. The full form needs none of these conditions: with the
observed sensitivity it estimates the influence function of Newey and
McFadden (1994, Theorem 6.1) whenever the contributions of each unit,
row or cluster, to the estimating functions of both stages have zero
mean at their limits and the units are independent. Where both routes
apply and the conditions hold they agree to first order; where the
conditions fail, the full form is the one to use, and `method = "full"`
asks for it. A bootstrap that resamples independent units or clusters
and refits both stages
([`mc_graph_bootstrap`](https://evandeilton.github.io/mcglmgraph/reference/mc_graph_bootstrap.md))
propagates the first stage by construction and is the reference against
which both routes can be checked.

## The full form

Outside the shortcut the three pieces are formed from the first stage
itself. The fitted MCGLM is rebuilt at its estimate with the routines of
mcglm, and split into contributions \\\boldsymbol\psi\_{\alpha,u}\\ of
its quasi-score and Pearson functions over units \\u\\; \\\mathbf A\\ is
its sensitivity; and \\\mathbf B\\ is formed by the chain rule on the
Pearson residual \\r\_{ri}=(y\_{ri}-\mu\_{ri})/\varsigma\_{ri}\\,
\\\varsigma\_{ri}^2=C\_{(ri),(ri)}\\, \$\$-\frac{\partial
r\_{ri}}{\partial\boldsymbol\alpha}=
\frac{1}{\varsigma\_{ri}}\frac{\partial\mu\_{ri}}{\partial\boldsymbol\alpha}
+r\_{ri}\frac{\partial\log\varsigma\_{ri}}{\partial\boldsymbol\alpha},\$\$
with no condition on the mean or on the residual covariance. The
correlations between responses of the first stage do not enter
\\\varsigma\_{ri}\\, so their block of \\\mathbf B\\ is zero.

The sensitivity \\\mathbf A\\ of the first stage is where the full form
departs from what mcglm reports. By default, `sensitivity = "observed"`,
it is the Jacobian at the estimate of the first-stage estimating
equations as mcglm solved them, its bias correction of the Pearson
equations included when the fit used it (`fit$con$correct`), formed by
central differences of their totals, each one rebuilt with the routines
of mcglm. Divided by the number of units it estimates the matrix the
influence function asks for whether or not the mean and the variance
function of the first stage are correctly specified. With
`sensitivity = "expected"` \\\mathbf A\\ takes its model-based form,
with blocks \\\mathbf D^\top\mathbf C^{-1}\mathbf W\mathbf D\\ for the
regression parameters, \\\operatorname{tr}(\mathbf C^{-1}\dot{\mathbf
C}\_\lambda\mathbf W \mathbf C^{-1}\dot{\mathbf C}\_{\lambda'})\\ for
the covariance parameters \\\lambda\\ and the same trace with
\\\partial\mathbf C/\partial \beta\\ in the place of \\\dot{\mathbf
C}\_{\lambda'}\\ for the covariance parameters on the regression
parameters, where \\\mathbf D\\ is the derivative of the stacked means,
\\\mathbf C\\ the first-stage covariance and \\\mathbf W\\ the diagonal
matrix of prior weights. It is cheaper and estimates the same limit when
the mean and the covariance of the first stage are correctly specified;
the corrected covariances of the two forms then differ by a term of
order at most \\n^{-1/2}\\ relative to the naive sandwich. For a
parameter whose corrected variance vanishes to first order, as that of
the innovation variance of a root node of a DAG does when the first
stage fixes the average of the squared residuals of that response and
their mean is zero at the limit, the variance reported under either form
is of lower order than the naive one, and the two forms can differ by as
much as their own size. When the mean or the variance function of the
first stage is misspecified the expected form in general no longer
estimates \\\mathbf A\\, the difference between the two corrected
covariances can be of first order, and the observed sensitivity is the
one to use. The observed form inherits the conditioning of the first
stage: when that stage is weakly identified, as a Poisson-Tweedie fit
whose estimated power and dispersion trade off against each other can
be, its Jacobian is close to singular and the propagated variance is
numerically unstable. The attribute `"sensitivity_rcond"` reports the
reciprocal condition number of the sensitivity used, and comparing the
two forms on the same fit is a direct check. Either way the block of the
covariance parameters on the regression parameters is kept: the
`joint_inv_sensitivity` of mcglm has it at zero, it is not zero when the
variance depends on the mean, and without it the regression columns of
\\\mathbf B\mathbf A^{-1}\\ do not vanish under the conditions of the
reduced form, where they must. The unit contributions of the first stage
are centered: the bias correction is one term for the whole sample
rather than a sum over units, so the contributions without it sum to
minus it at the estimate, and a term common to all units does not enter
the variability.

The first stage is a sum over units only when its covariance is block
diagonal over them. With an identity within-response structure the units
are the rows, or the clusters when `cluster` is given. With any other
structure `cluster` is required and has to nest it: every matrix of the
matrix linear predictor must vanish between distinct clusters, which is
checked.

`method = "auto"` reads the route from what the fit shows. It takes the
full form when the full form is available and at least one response has
a variance function other than `"constant"`. With a variance that
depends on the mean the shortcut needs condition (c) of the section
"What the correction reaches", at least in its weaker form, and that is
a statement about the data that cannot be read off the fit. It also
takes the full form, when it is available, on a fitted DAG with an
immorality: the shortcut evaluates the reduced form at the fitted
covariance, which is exact for the parameters of every node of a DAG
without immoralities but, for the parameters of a collider with
non-adjacent parents and of its descendants, only a first-order
approximation under correct specification. In every other case `"auto"`
takes the shortcut when it applies and the full form when it does not.
The shortcut, when it applies, is therefore the route in two cases:
every response has the constant variance function and the layer is a GGM
or a DAG without immoralities, or the full form is not available. On a
DAG with an immorality the second case raises a warning. The variance
function is read by its name in the first-stage fit, one per response,
whatever its power: a Tweedie variance with power fixed at zero counts
as depending on the mean, and so does a fit whose names cannot be read.
`method = "shortcut"` and `method = "full"` force a route, and an
explicit `"shortcut"` is honored without a warning. On the fits both
routes cover they agree to first order when the conditions of the
section "When the shortcut is right" hold. The attribute
`"two_stage_method"` of the result records the route taken, and
[`mc_two_stage_available`](https://evandeilton.github.io/mcglmgraph/reference/mc_two_stage_available.md)
reports it beforehand; a fit neither route covers raises an error rather
than returning the uncorrected matrix under the corrected name.

## What the correction reaches

How much of the propagation reaches a given quantity is settled, for the
undirected layer, by a diagonal-rescaling argument, under three
conditions on the first stage, which are sufficient and, as the section
"When the shortcut is right" says, not always necessary:

- (a):

  the Pearson residual has the form
  \\r\_{ri}=(y\_{ri}-\mu\_{ri})/\varsigma\_{ri}\\, with the mean and the
  scale depending only on the first-stage parameters and the covariates,
  which is what mcglm returns;

- (b):

  the mean is correctly specified, \\E(r\_{ri}\mid\mathbf x_i)=0\\;

- (c):

  the conditional residual covariance does not depend on the covariates,
  \\E(\mathbf r_i\mathbf r_i^\top\mid\mathbf x_i)=
  \boldsymbol\Sigma_b\\.

Under the three, the propagation term is the coordinate vector of
\\\mathbf K_b\boldsymbol\Lambda_i+\boldsymbol\Lambda_i\mathbf K_b\\,
with \\\boldsymbol\Lambda_i\\ diagonal: the first stage acts on the
residual covariance, to first order, as a diagonal rescaling. The
propagation vanishes on every edge whose precision entry is zero and is
annihilated by the gradient of every partial correlation, so for those
quantities the sandwich and the corrected form agree to first order. The
diagonal coefficients and the fitted \\\boldsymbol\Sigma_b\\ are not
invariant and do need the correction. The DAG layer has the same
reading, for any fixed DAG, with the second moment of the residuals in
the place of \\\boldsymbol\Sigma_b\\ in (c): the propagation vanishes on
every structural coefficient whose target value is zero and on every
function of the parameters that is invariant to a diagonal rescaling of
the residuals, such as the partial correlations, while the nonzero
structural coefficients and the innovation variances need the
correction.

For the shortcut the agreement on the partial correlations is exact and
holds whatever the data: the gradients of the partial correlations
annihilate \\\mathbf H^{-1}\mathbf Q\\ algebraically, since \\2\mathbf
H^{-1}\mathbf Q\_{\cdot r}\\ is the direction of a rescaling, so
`type = "sandwich"` and the shortcut give them the same variance. That
says nothing about whether the variance is right. When the first stage
does not act as a rescaling, as can happen with a variance function that
depends on the mean and a residual correlation that moves with the
covariates, both can miss the same first-order term on the partial
correlations, and only the full form accounts for it.

The model-based form \\\mathbf J_n^{-1}\\ agrees with the empirical
forms only when, in addition, \\\mathbf H=\mathbf J\\ and \\\mathbf
V=\mathbf J\\; the latter is the correctly specified Gaussian
fourth-moment regime. Outside it the difference is not a matter of the
first stage, and (a) to (c) do not remove it. The package checks none of
(a) to (c): they are statements about the process that generated the
data. Applied to a support chosen by a search over graphs, every form
remains descriptive, since the selection is not incorporated.

## References

Bonat, W. H. and Jorgensen, B. (2016). Multivariate covariance
generalized linear models. *Journal of the Royal Statistical Society,
Series C* 65(5), 649–675.
[doi:10.1111/rssc.12145](https://doi.org/10.1111/rssc.12145)

Godambe, V. P. (1960). An optimum property of regular maximum likelihood
estimation. *Annals of Mathematical Statistics* 31(4), 1208–1211.
[doi:10.1214/aoms/1177705693](https://doi.org/10.1214/aoms/1177705693)

Jorgensen, B. and Knudsen, S. J. (2004). Parameter orthogonality and
bias adjustment for estimating functions. *Scandinavian Journal of
Statistics* 31(1), 93–114.
[doi:10.1111/j.1467-9469.2004.00375.x](https://doi.org/10.1111/j.1467-9469.2004.00375.x)

Murphy, K. M. and Topel, R. H. (1985). Estimation and inference in
two-step econometric models. *Journal of Business and Economic
Statistics* 3(4), 370–379.
[doi:10.1080/07350015.1985.10509471](https://doi.org/10.1080/07350015.1985.10509471)

Newey, W. K. and McFadden, D. (1994). Large sample estimation and
hypothesis testing. In Engle, R. F. and McFadden, D. L. (eds), *Handbook
of Econometrics*, vol. 4, 2111–2245. Elsevier.
[doi:10.1016/S1573-4412(05)80005-4](https://doi.org/10.1016/S1573-4412%2805%2980005-4)

White, H. (1980). A heteroskedasticity-consistent covariance matrix
estimator and a direct test for heteroskedasticity. *Econometrica*
48(4), 817–838. [doi:10.2307/1912934](https://doi.org/10.2307/1912934)

## See also

[`mc_test_independence`](https://evandeilton.github.io/mcglmgraph/reference/mc_test_independence.md)
for the joint Wald test, which consumes this matrix,
[`mc_wald_graph`](https://evandeilton.github.io/mcglmgraph/reference/mc_wald_graph.md)
for general linear hypotheses,
[`mc_two_stage_available`](https://evandeilton.github.io/mcglmgraph/reference/mc_two_stage_available.md)
for the route the correction takes,
[`mc_graph_bootstrap`](https://evandeilton.github.io/mcglmgraph/reference/mc_graph_bootstrap.md)
for resampling, and
[`mc_partial_cor_delta`](https://evandeilton.github.io/mcglmgraph/reference/mc_partial_cor_delta.md)
for the delta-method inference on the partial correlations.

Other inference on the graph layer:
[`mc_test_adequacy()`](https://evandeilton.github.io/mcglmgraph/reference/mc_test_adequacy.md),
[`mc_test_independence()`](https://evandeilton.github.io/mcglmgraph/reference/mc_test_independence.md),
[`mc_two_stage_available()`](https://evandeilton.github.io/mcglmgraph/reference/mc_two_stage_available.md),
[`mc_wald_graph()`](https://evandeilton.github.io/mcglmgraph/reference/mc_wald_graph.md)

## Examples

``` r
# \donttest{
set.seed(1); n <- 300
K <- diag(3) * 2; K[1, 2] <- K[2, 1] <- -0.7
Y <- matrix(rnorm(n * 3), n) %*% chol(solve(K))
d <- data.frame(y1 = Y[, 1], y2 = Y[, 2], y3 = Y[, 3], x = rnorm(n))
f <- mcglm_graph(list(y1 ~ x, y2 ~ x, y3 ~ x), data = d)
#> Automatic initial values selected.
V <- mc_vcov_graph(f, type = "sandwich")
sqrt(diag(V))[attr(V, "edge_pos")]      # standard errors of the edges
#>     eta_12     eta_13     eta_23 
#> 0.07693955 0.06668838 0.05918675 
# }
```
