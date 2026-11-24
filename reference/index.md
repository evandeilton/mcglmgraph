# Package index

## Specify and fit

- [`mc_between_graph()`](https://evandeilton.github.io/mcglmgraph/reference/mc_between_graph.md)
  : Specify the between-response graph layer
- [`mcglm_graph()`](https://evandeilton.github.io/mcglmgraph/reference/mcglm_graph.md)
  [`mcglmgraph()`](https://evandeilton.github.io/mcglmgraph/reference/mcglm_graph.md)
  [`print(`*`<mcglm_graph>`*`)`](https://evandeilton.github.io/mcglmgraph/reference/mcglm_graph.md)
  [`summary(`*`<mcglm_graph>`*`)`](https://evandeilton.github.io/mcglmgraph/reference/mcglm_graph.md)
  [`print(`*`<summary.mcglm_graph>`*`)`](https://evandeilton.github.io/mcglmgraph/reference/mcglm_graph.md)
  [`coef(`*`<mcglm_graph>`*`)`](https://evandeilton.github.io/mcglmgraph/reference/mcglm_graph.md)
  [`plot(`*`<mcglm_graph>`*`)`](https://evandeilton.github.io/mcglmgraph/reference/mcglm_graph.md)
  : Fit a graph-structured multivariate covariance generalized linear
  model

## Read the graph layer

- [`mc_partial_cor()`](https://evandeilton.github.io/mcglmgraph/reference/mc_partial_cor.md)
  : Partial correlations between responses
- [`mc_partial_cor_delta()`](https://evandeilton.github.io/mcglmgraph/reference/mc_partial_cor_delta.md)
  : Delta-method inference for the between-response partial correlations
- [`mc_dag_coefficients()`](https://evandeilton.github.io/mcglmgraph/reference/mc_dag_coefficients.md)
  [`print(`*`<mcg_dag_coef>`*`)`](https://evandeilton.github.io/mcglmgraph/reference/mc_dag_coefficients.md)
  : Scale-free coefficients of a fitted directed layer
- [`mc_get_graph()`](https://evandeilton.github.io/mcglmgraph/reference/mc_get_graph.md)
  : Extract the estimated dependence graph
- [`mc_plot_network()`](https://evandeilton.github.io/mcglmgraph/reference/mc_plot_network.md)
  : Plot the estimated dependence network

## Inference on the graph layer

- [`mc_vcov_graph()`](https://evandeilton.github.io/mcglmgraph/reference/mc_vcov_graph.md)
  : Covariance matrix of the graph-layer parameters
- [`mc_two_stage_available()`](https://evandeilton.github.io/mcglmgraph/reference/mc_two_stage_available.md)
  : Is the two-stage correction available for a fit
- [`mc_wald_graph()`](https://evandeilton.github.io/mcglmgraph/reference/mc_wald_graph.md)
  : Wald test of a linear hypothesis on the graph layer
- [`mc_test_independence()`](https://evandeilton.github.io/mcglmgraph/reference/mc_test_independence.md)
  : Wald test of between-response conditional association
- [`mc_test_adequacy()`](https://evandeilton.github.io/mcglmgraph/reference/mc_test_adequacy.md)
  : Adequacy of a support against the saturated graph
- [`mc_edge_scores()`](https://evandeilton.github.io/mcglmgraph/reference/mc_edge_scores.md)
  : Edge score and Wald statistics

## Compare supports

- [`mc_compute_criterion()`](https://evandeilton.github.io/mcglmgraph/reference/mc_compute_criterion.md)
  : Information criteria for a graph-structured fit
- [`mc_graph_ploglik()`](https://evandeilton.github.io/mcglmgraph/reference/mc_graph_ploglik.md)
  : Gaussian working log-criterion of the between-response fit
- [`mc_graph_cv()`](https://evandeilton.github.io/mcglmgraph/reference/mc_graph_cv.md)
  [`print(`*`<mcg_cv>`*`)`](https://evandeilton.github.io/mcglmgraph/reference/mc_graph_cv.md)
  : Out-of-sample comparison of candidate supports

## Select an undirected graph

- [`mc_graph_stepwise()`](https://evandeilton.github.io/mcglmgraph/reference/mc_graph_stepwise.md)
  : Stepwise selection of the dependence graph
- [`mc_select_lambda()`](https://evandeilton.github.io/mcglmgraph/reference/mc_select_lambda.md)
  : Select the penalty parameter along an L1 path
- [`mc_graph_penalized()`](https://evandeilton.github.io/mcglmgraph/reference/mc_graph_penalized.md)
  : L1-penalized estimation of the precision graph
- [`print(`*`<mcg_selection>`*`)`](https://evandeilton.github.io/mcglmgraph/reference/mcg_selection-methods.md)
  [`summary(`*`<mcg_selection>`*`)`](https://evandeilton.github.io/mcglmgraph/reference/mcg_selection-methods.md)
  [`plot(`*`<mcg_selection>`*`)`](https://evandeilton.github.io/mcglmgraph/reference/mcg_selection-methods.md)
  : Methods for graph structure selection results

## Discover a directed structure

- [`mc_graph_pipeline()`](https://evandeilton.github.io/mcglmgraph/reference/mc_graph_pipeline.md)
  [`print(`*`<mcg_cpdag>`*`)`](https://evandeilton.github.io/mcglmgraph/reference/mc_graph_pipeline.md)
  : Constraint-based CPDAG discovery
- [`mc_ges()`](https://evandeilton.github.io/mcglmgraph/reference/mc_ges.md)
  [`print(`*`<mcg_ges>`*`)`](https://evandeilton.github.io/mcglmgraph/reference/mc_ges.md)
  : Greedy equivalence search on the residual covariance
- [`mc_pc_skeleton()`](https://evandeilton.github.io/mcglmgraph/reference/mc_pc_skeleton.md)
  : PC-style skeleton with separating sets
- [`mc_separating_sets()`](https://evandeilton.github.io/mcglmgraph/reference/mc_separating_sets.md)
  [`print(`*`<mcg_sepsets>`*`)`](https://evandeilton.github.io/mcglmgraph/reference/mc_separating_sets.md)
  : All separating sets of a pair of responses
- [`mc_find_vstructures()`](https://evandeilton.github.io/mcglmgraph/reference/mc_find_vstructures.md)
  : Detect v-structures in the constraint-based skeleton
- [`mc_orient_meek()`](https://evandeilton.github.io/mcglmgraph/reference/mc_orient_meek.md)
  : Apply Meek's orientation rules
- [`mc_residual_df()`](https://evandeilton.github.io/mcglmgraph/reference/mc_residual_df.md)
  [`print(`*`<mcg_residual_df>`*`)`](https://evandeilton.github.io/mcglmgraph/reference/mc_residual_df.md)
  : Residual degrees of freedom of the first stage
- [`mc_skeleton_to_bnlearn()`](https://evandeilton.github.io/mcglmgraph/reference/mc_skeleton_to_bnlearn.md)
  : Convert a skeleton to the set of allowed bnlearn arcs

## Equivalence classes and directed fits

- [`mc_class_members()`](https://evandeilton.github.io/mcglmgraph/reference/mc_class_members.md)
  : Members of a Markov equivalence class
- [`mc_class_coefficients()`](https://evandeilton.github.io/mcglmgraph/reference/mc_class_coefficients.md)
  [`print(`*`<mcg_class_coef>`*`)`](https://evandeilton.github.io/mcglmgraph/reference/mc_class_coefficients.md)
  : What the members of a class disagree about
- [`mc_dag_closed_root()`](https://evandeilton.github.io/mcglmgraph/reference/mc_dag_closed_root.md)
  : Closed-form root of the directed layer

## Stability and resampling

- [`mc_graph_stability()`](https://evandeilton.github.io/mcglmgraph/reference/mc_graph_stability.md)
  [`print(`*`<mcg_stability>`*`)`](https://evandeilton.github.io/mcglmgraph/reference/mc_graph_stability.md)
  : Selection-frequency aggregation of the between-response graph
  structure
- [`mc_cpdag_stability()`](https://evandeilton.github.io/mcglmgraph/reference/mc_cpdag_stability.md)
  [`print(`*`<mcg_cpdag_stability>`*`)`](https://evandeilton.github.io/mcglmgraph/reference/mc_cpdag_stability.md)
  : Resampling frequencies of the skeleton and of each orientation
- [`mc_graph_bootstrap()`](https://evandeilton.github.io/mcglmgraph/reference/mc_graph_bootstrap.md)
  [`print(`*`<mcg_bootstrap>`*`)`](https://evandeilton.github.io/mcglmgraph/reference/mc_graph_bootstrap.md)
  : Bootstrap the graph layer with both stages refitted

## Check against a known truth

- [`mc_structure_metrics()`](https://evandeilton.github.io/mcglmgraph/reference/mc_structure_metrics.md)
  : Recovery of an undirected support
- [`mc_cpdag_metrics()`](https://evandeilton.github.io/mcglmgraph/reference/mc_cpdag_metrics.md)
  : Recovery of a partially directed structure
- [`mc_true_cpdag()`](https://evandeilton.github.io/mcglmgraph/reference/mc_true_cpdag.md)
  : CPDAG of a known directed structure
- [`mc_sid()`](https://evandeilton.github.io/mcglmgraph/reference/mc_sid.md)
  [`print(`*`<mcg_sid>`*`)`](https://evandeilton.github.io/mcglmgraph/reference/mc_sid.md)
  : Structural intervention distance

## Covariance algebra

Building blocks of the two parametrizations.

- [`mc_precision_bases()`](https://evandeilton.github.io/mcglmgraph/reference/mc_precision_bases.md)
  : Build the base matrices of the GGM precision predictor
- [`mc_build_precision()`](https://evandeilton.github.io/mcglmgraph/reference/mc_build_precision.md)
  : Assemble the between-response precision matrix
- [`mc_precision_to_covariance()`](https://evandeilton.github.io/mcglmgraph/reference/mc_precision_to_covariance.md)
  : Convert a precision matrix into the between-response covariance
- [`mc_precision_derivatives()`](https://evandeilton.github.io/mcglmgraph/reference/mc_precision_derivatives.md)
  : Analytical derivatives of the covariance with respect to eta
- [`mc_spd_line_search()`](https://evandeilton.github.io/mcglmgraph/reference/mc_spd_line_search.md)
  : Backtracking line search preserving positive-definiteness
- [`mc_check_spd()`](https://evandeilton.github.io/mcglmgraph/reference/mc_check_spd.md)
  : Check positive-definiteness of a matrix
- [`mc_build_T_matrix()`](https://evandeilton.github.io/mcglmgraph/reference/mc_build_T_matrix.md)
  : Build the structural matrix T of the LDL parametrization
- [`mc_build_D_matrix()`](https://evandeilton.github.io/mcglmgraph/reference/mc_build_D_matrix.md)
  : Build the diagonal innovation-variance matrix D
- [`mc_cholesky_to_covariance()`](https://evandeilton.github.io/mcglmgraph/reference/mc_cholesky_to_covariance.md)
  : Covariance from the modified Cholesky (LDL) parametrization
- [`mc_cholesky_to_precision()`](https://evandeilton.github.io/mcglmgraph/reference/mc_cholesky_to_precision.md)
  : Precision from the modified Cholesky (LDL) parametrization
- [`mc_cholesky_derivatives()`](https://evandeilton.github.io/mcglmgraph/reference/mc_cholesky_derivatives.md)
  : Analytical LDL derivatives of the covariance
- [`mc_dag_sensitivity()`](https://evandeilton.github.io/mcglmgraph/reference/mc_dag_sensitivity.md)
  : Block-diagonal sensitivity of the directed layer
- [`mc_dag_unit_diagonal()`](https://evandeilton.github.io/mcglmgraph/reference/mc_dag_unit_diagonal.md)
  : Innovation variances that force a unit diagonal
- [`mc_build_sigma_between_graph()`](https://evandeilton.github.io/mcglmgraph/reference/mc_build_sigma_between_graph.md)
  : Build the between-response covariance and its derivatives under a
  graph

## Package

- [`mcglmgraph-package`](https://evandeilton.github.io/mcglmgraph/reference/mcglmgraph-package.md)
  : mcglmgraph: Graphical Multivariate Covariance Generalized Linear
  Models
