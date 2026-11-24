# mcglmgraph: Graphical Multivariate Covariance Generalized Linear Models

Graph-structured covariance between responses for Multivariate
Covariance Generalized Linear Models (MCGLM; Bonat and Jorgensen, 2016),
fitted by two-stage residual analysis. The method is called MCGLM-Graph.

## Two stages

The first stage fits the marginal MCGLM, that is, the mean and
dispersion structure of each response, through
[`mcglm::mcglm()`](https://rdrr.io/pkg/mcglm/man/mcglm.html), and
extracts the standardized Pearson residuals. The second stage imposes a
graph on the between-response covariance \\\Sigma_b\\ of these residuals
and estimates the graph parameters with the second-order (Pearson)
estimating functions of the MCGLM, with analytical derivatives
implemented in C++ through RcppArmadillo. No joint likelihood is used.
The regression coefficients are those of the first stage. Graph edges
therefore describe residual association among responses after the
marginal means are removed. Variances of the graph parameters are
conditional on the first stage unless the two-stage variance, which
propagates the first-stage uncertainty (Murphy and Topel, 1985), is
requested
([`mc_vcov_graph`](https://evandeilton.github.io/mcglmgraph/reference/mc_vcov_graph.md)).

## Three layers

- GGM:

  Gaussian graphical model. The precision matrix follows a linear
  predictor on known symmetric basis matrices, \\K_b = \Sigma_b^{-1} =
  \sum_d \eta_d Z_d\\. A zero off-diagonal entry is a null partial
  correlation (Lauritzen, 1996), that is, absence of second-order
  conditional association between two responses given the others. This
  coincides with conditional independence only under joint normality of
  the analyzed residuals.

- DAG:

  Directed acyclic graph. The modified Cholesky (LDL) decomposition
  (Pourahmadi, 1999) \\\Sigma_b = T^{-1} D (T^{-1})^\top\\ under a fixed
  ordering of the responses, with \\T\\ unit lower triangular and \\D\\
  diagonal, guarantees positive definiteness. A causal reading requires
  a substantively justified ordering and further assumptions.

- CPDAG:

  Completed partially directed acyclic graph, the Markov equivalence
  class of DAGs. The skeleton is learned by partial-correlation tests
  and the identifiable edges are oriented by v-structures and Meek's
  rules (Spirtes, Glymour and Scheines, 2000; Meek, 1995);
  alternatively, the class is searched by greedy equivalence search
  (GES; Chickering, 2002).

## Functions by task

- Specify and fit:

  [`mc_between_graph`](https://evandeilton.github.io/mcglmgraph/reference/mc_between_graph.md),
  [`mcglm_graph`](https://evandeilton.github.io/mcglmgraph/reference/mcglm_graph.md)
  (alias
  [`mcglmgraph()`](https://evandeilton.github.io/mcglmgraph/reference/mcglm_graph.md)).

- Read the graph layer:

  [`mc_partial_cor`](https://evandeilton.github.io/mcglmgraph/reference/mc_partial_cor.md),
  [`mc_partial_cor_delta`](https://evandeilton.github.io/mcglmgraph/reference/mc_partial_cor_delta.md),
  [`mc_dag_coefficients`](https://evandeilton.github.io/mcglmgraph/reference/mc_dag_coefficients.md),
  [`mc_get_graph`](https://evandeilton.github.io/mcglmgraph/reference/mc_get_graph.md),
  [`mc_plot_network`](https://evandeilton.github.io/mcglmgraph/reference/mc_plot_network.md).

- Inference on the graph layer:

  [`mc_vcov_graph`](https://evandeilton.github.io/mcglmgraph/reference/mc_vcov_graph.md),
  [`mc_two_stage_available`](https://evandeilton.github.io/mcglmgraph/reference/mc_two_stage_available.md),
  [`mc_wald_graph`](https://evandeilton.github.io/mcglmgraph/reference/mc_wald_graph.md),
  [`mc_test_independence`](https://evandeilton.github.io/mcglmgraph/reference/mc_test_independence.md),
  [`mc_test_adequacy`](https://evandeilton.github.io/mcglmgraph/reference/mc_test_adequacy.md),
  [`mc_edge_scores`](https://evandeilton.github.io/mcglmgraph/reference/mc_edge_scores.md).

- Compare supports:

  [`mc_compute_criterion`](https://evandeilton.github.io/mcglmgraph/reference/mc_compute_criterion.md),
  [`mc_graph_ploglik`](https://evandeilton.github.io/mcglmgraph/reference/mc_graph_ploglik.md),
  [`mc_graph_cv`](https://evandeilton.github.io/mcglmgraph/reference/mc_graph_cv.md).

- Select an undirected graph:

  [`mc_graph_stepwise`](https://evandeilton.github.io/mcglmgraph/reference/mc_graph_stepwise.md),
  [`mc_select_lambda`](https://evandeilton.github.io/mcglmgraph/reference/mc_select_lambda.md),
  [`mc_graph_penalized`](https://evandeilton.github.io/mcglmgraph/reference/mc_graph_penalized.md).

- Discover a directed structure:

  [`mc_graph_pipeline`](https://evandeilton.github.io/mcglmgraph/reference/mc_graph_pipeline.md),
  [`mc_ges`](https://evandeilton.github.io/mcglmgraph/reference/mc_ges.md),
  [`mc_pc_skeleton`](https://evandeilton.github.io/mcglmgraph/reference/mc_pc_skeleton.md),
  [`mc_separating_sets`](https://evandeilton.github.io/mcglmgraph/reference/mc_separating_sets.md),
  [`mc_find_vstructures`](https://evandeilton.github.io/mcglmgraph/reference/mc_find_vstructures.md),
  [`mc_orient_meek`](https://evandeilton.github.io/mcglmgraph/reference/mc_orient_meek.md),
  [`mc_residual_df`](https://evandeilton.github.io/mcglmgraph/reference/mc_residual_df.md),
  [`mc_skeleton_to_bnlearn`](https://evandeilton.github.io/mcglmgraph/reference/mc_skeleton_to_bnlearn.md).

- Equivalence classes and directed fits:

  [`mc_class_members`](https://evandeilton.github.io/mcglmgraph/reference/mc_class_members.md),
  [`mc_class_coefficients`](https://evandeilton.github.io/mcglmgraph/reference/mc_class_coefficients.md),
  [`mc_dag_closed_root`](https://evandeilton.github.io/mcglmgraph/reference/mc_dag_closed_root.md).

- Stability and resampling:

  [`mc_graph_stability`](https://evandeilton.github.io/mcglmgraph/reference/mc_graph_stability.md),
  [`mc_cpdag_stability`](https://evandeilton.github.io/mcglmgraph/reference/mc_cpdag_stability.md),
  [`mc_graph_bootstrap`](https://evandeilton.github.io/mcglmgraph/reference/mc_graph_bootstrap.md).

- Check against a known truth:

  [`mc_structure_metrics`](https://evandeilton.github.io/mcglmgraph/reference/mc_structure_metrics.md),
  [`mc_cpdag_metrics`](https://evandeilton.github.io/mcglmgraph/reference/mc_cpdag_metrics.md),
  [`mc_true_cpdag`](https://evandeilton.github.io/mcglmgraph/reference/mc_true_cpdag.md),
  [`mc_sid`](https://evandeilton.github.io/mcglmgraph/reference/mc_sid.md).

- Covariance algebra:

  [`mc_precision_bases`](https://evandeilton.github.io/mcglmgraph/reference/mc_precision_bases.md),
  [`mc_build_precision`](https://evandeilton.github.io/mcglmgraph/reference/mc_build_precision.md),
  [`mc_precision_to_covariance`](https://evandeilton.github.io/mcglmgraph/reference/mc_precision_to_covariance.md),
  [`mc_precision_derivatives`](https://evandeilton.github.io/mcglmgraph/reference/mc_precision_derivatives.md),
  [`mc_spd_line_search`](https://evandeilton.github.io/mcglmgraph/reference/mc_spd_line_search.md),
  [`mc_check_spd`](https://evandeilton.github.io/mcglmgraph/reference/mc_check_spd.md),
  [`mc_build_T_matrix`](https://evandeilton.github.io/mcglmgraph/reference/mc_build_T_matrix.md),
  [`mc_build_D_matrix`](https://evandeilton.github.io/mcglmgraph/reference/mc_build_D_matrix.md),
  [`mc_cholesky_to_covariance`](https://evandeilton.github.io/mcglmgraph/reference/mc_cholesky_to_covariance.md),
  [`mc_cholesky_to_precision`](https://evandeilton.github.io/mcglmgraph/reference/mc_cholesky_to_precision.md),
  [`mc_cholesky_derivatives`](https://evandeilton.github.io/mcglmgraph/reference/mc_cholesky_derivatives.md),
  [`mc_dag_sensitivity`](https://evandeilton.github.io/mcglmgraph/reference/mc_dag_sensitivity.md),
  [`mc_dag_unit_diagonal`](https://evandeilton.github.io/mcglmgraph/reference/mc_dag_unit_diagonal.md),
  [`mc_build_sigma_between_graph`](https://evandeilton.github.io/mcglmgraph/reference/mc_build_sigma_between_graph.md).

Selection results have `print`, `summary` and `plot` methods
([`mcg_selection-methods`](https://evandeilton.github.io/mcglmgraph/reference/mcg_selection-methods.md)).

## References

Bonat, W. H. and Jorgensen, B. (2016). Multivariate covariance
generalized linear models. *Journal of the Royal Statistical Society:
Series C*, 65(5), 649–675.
[doi:10.1111/rssc.12145](https://doi.org/10.1111/rssc.12145)

Bonat, W. H. (2018). Multiple response variables regression models in R:
the mcglm package. *Journal of Statistical Software*, 84(4), 1–30.
[doi:10.18637/jss.v084.i04](https://doi.org/10.18637/jss.v084.i04)

Chickering, D. M. (2002). Optimal structure identification with greedy
search. *Journal of Machine Learning Research*, 3, 507–554.

Lauritzen, S. L. (1996). *Graphical Models*. Clarendon Press, Oxford.

Meek, C. (1995). Causal inference and causal explanation with background
knowledge. In *Proceedings of the Eleventh Conference on Uncertainty in
Artificial Intelligence*, 403–410.

Murphy, K. M. and Topel, R. H. (1985). Estimation and inference in
two-step econometric models. *Journal of Business & Economic
Statistics*, 3(4), 370–379.
[doi:10.1080/07350015.1985.10509471](https://doi.org/10.1080/07350015.1985.10509471)

Pourahmadi, M. (1999). Joint mean-covariance models with applications to
longitudinal data: unconstrained parameterisation. *Biometrika*, 86(3),
677–690.
[doi:10.1093/biomet/86.3.677](https://doi.org/10.1093/biomet/86.3.677)

Spirtes, P., Glymour, C. and Scheines, R. (2000). *Causation,
Prediction, and Search*, 2nd edition. MIT Press, Cambridge,
Massachusetts.

## Author

**Maintainer**: Jose Evandeilton Lopes <evandeilton@gmail.com>

Authors:

- Jose Evandeilton Lopes <evandeilton@gmail.com>

- Wagner Hugo Bonat <wbonat@ufpr.br> (Author of the MCGLM framework that
  this package extends) \[copyright holder\]
