#' Fit a (Zero-Inflated) Bayesian Causal Forest for count regression.
#'
#' \code{countbcf} is the sampler behind the CountBCF and Zero-Inflated CountBCF
#' models of Souto (2026a, 2026b). It estimates Bayesian Causal Forests for
#' count and zero-inflated count outcomes by combining the log-linear BART
#' backbone of Murray (2021) with the (mu, tau) decomposition of Hahn,
#' Murray and Carvalho (2020). The C++ MCMC machinery, GIG leaf prior, and
#' zero-inflation bookkeeping are inherited from the \code{countbart} source
#' code of Wikle and Zigler (2023).
#'
#' @details
#' For a count outcome \eqn{y_i}, binary treatment indicator
#' \eqn{z_i \in \{0,1\}}, and covariates \eqn{x_i}:
#'
#' Non-zero-inflated count models (\code{"poisson"}, \code{"nb"}):
#' \deqn{\log E[y_i \mid x_i, z_i] = \mu_f(x_i, \hat\pi_i) + z_i\, \tau_f(x_i),}
#' where \eqn{\mu_f} is a "vanilla" prognostic forest and \eqn{\tau_f} is a
#' moderating forest entered linearly in \eqn{z}.
#'
#' Zero-inflated count models (\code{"zipoisson"}, \code{"zinb"}): the rate
#' component above is augmented with two log-odds components for the
#' structural-zero indicator \eqn{S_i}:
#' \deqn{\mathrm{logit}\, P(S_i = 1 \mid x_i, z_i) =
#'   \bigl(\mu_{f_0}(x_i, \hat\pi_i) + z_i\, \tau_{f_0}(x_i)\bigr)
#'   - \bigl(\mu_{f_1}(x_i, \hat\pi_i) + z_i\, \tau_{f_1}(x_i)\bigr),}
#' with \eqn{Y_i = 0} when \eqn{S_i = 1} and
#' \eqn{Y_i \sim \mathrm{Poisson}(\lambda_i)} (or
#' \eqn{\mathrm{NegBin}(\lambda_i, \kappa)}) otherwise. The total is therefore
#' 2 forests for non-ZI models and 6 for ZI models.
#'
#' The leaf prior is the GIG mixture from Murray (2021), with
#' concentration \code{a0} (count component) and \code{z_conc}
#' (zero-inflation components). The \eqn{\mu} and \eqn{\tau} forests within a
#' component share this leaf concentration; the stronger shrinkage of the
#' \eqn{\tau} (treatment-effect) forests is delivered by the tree-count ratio
#' and the stub-biased tree-shape prior rather than by a separate
#' concentration. With the leaf scale calibrated for \code{ntree_control}
#' trees but applied to the smaller \code{ntree_moderate}-tree \eqn{\tau}
#' forest, the implied prior on the treatment-effect function has standard
#' deviation \code{a0 * sqrt(ntree_moderate / ntree_control)}, so heterogeneous
#' treatment effects are more strongly shrunk than the prognostic component.
#'
#' Data are internally sorted by \code{y} so that zeros come first; the
#' permutation is returned as \code{order_vec} and is applied to every
#' observation-indexed input, including the rows of \code{randeff_design}.
#' The reported
#' \code{yhat}, \code{yhat_log}, and \code{tau_f_log} rows are in the
#' sorted order. The per-forest coefficient matrices (\code{mu_f_post},
#' \code{tau_f_post}, etc.) are unsorted (in the original input order)
#' and are typically what users want for downstream CATE / ATE inference.
#'
#' @references
#' Souto, H. G. (2026a). \emph{CountBCF: Bayesian Causal Forests for
#'   Count Outcomes.} Working paper.
#'
#' Souto, H. G. (2026b). \emph{Zero-Inflated CountBCF: Bayesian Causal
#'   Forests for Zero-Inflated Count Outcomes.} Working paper.
#'
#' Hahn, P. R., Murray, J. S., and Carvalho, C. M. (2020). Bayesian
#'   regression tree models for causal inference: regularization,
#'   confounding, and heterogeneous effects. \emph{Bayesian Analysis},
#'   15(3), 965--1056. \doi{10.1214/19-BA1195}.
#'
#' Murray, J. S. (2021). Log-linear Bayesian additive regression trees for
#'   multinomial logistic and count regression models. \emph{Journal of
#'   the American Statistical Association}, 116(534), 756--769.
#'   \doi{10.1080/01621459.2020.1813587}.
#'
#' Wikle, N. B., and Zigler, C. M. (2023). Causal health impacts of power
#'   plant emission controls under modeled and uncertain physical process
#'   interference. \emph{Annals of Applied Statistics}.
#'   Companion package: \url{https://github.com/nbwikle/estimating-interference}.
#'
#' @param y Response variable (vector of non-negative integers).
#' @param z Treatment assignment (vector of 0s and 1s).
#' @param x_control Design matrix for the prognostic function \eqn{\mu^{(0)}(x)}.
#' @param x_moderate Design matrix for the treatment effect function \eqn{\tau^{(0)}(x)}. Defaults to x_control.
#' @param x_zero Design matrix for the zero-inflation odds component \eqn{f^{(0)}}. Defaults to x_control.
#' @param x_pos Design matrix for the not-zero-inflation odds component \eqn{f^{(1)}}. Defaults to x_control.
#' @param pihat Length n estimate of \eqn{E(Z|X)}, included in mu-forest covariates. Defaults to 0.5.
#' @param offset Offset value (scalar or vector). If NULL, an offset of 1 is used.
#' @param nburn Number of burn-in MCMC iterations.
#' @param nsim Number of MCMC iterations to save after burn-in.
#' @param nthin Save every nthin'th MCMC iterate.
#' @param update_interval Print status every update_interval MCMC iterations.
#' @param ntree_control Number of trees in mu^{(0)} (count component).
#' @param ntree_moderate Number of trees in tau^{(0)} (count component).
#' @param nztree_control Number of trees in mu^{(0)} for zero-inflation odds (used in ZI models).
#' @param nztree_moderate Number of trees in tau^{(0)} for zero-inflation odds (used in ZI models).
#' @param a0 Concentration parameter for count-component leaf prior (shared by
#'   the mu and tau forests). NA gives a default estimate.
#' @param z_conc Concentration parameter for zero-inflation odds leaf prior
#'   (shared by the mu and tau forests).
#' @param base_control Base for tree prior on mu^{(0)} trees.
#' @param power_control Power for tree prior on mu^{(0)} trees.
#' @param base_moderate Base for tree prior on tau^{(0)} trees.
#' @param power_moderate Power for tree prior on tau^{(0)} trees.
#' @param kappa_a Shape parameter alpha for kappa (NB dispersion) beta-prime prior.
#' @param kappa_b Shape parameter beta for kappa (NB dispersion) beta-prime prior.
#' @param kappa_prop_sd SD of kappa M-H proposal.
#' @param count_model One of "poisson", "nb", "zipoisson", "zinb".
#' @param include_pihat One of "control", "moderate", "both", "none". Whether to
#'   include pihat among covariates of mu (control) or tau (moderate) forests.
#' @param randeff_design,randeff_variance_component_design,randeff_scales,randeff_df
#'   Random-effects setup; the 1x1 matrix defaults disable random effects. When
#'   used, \code{randeff_design} must have one row per observation (in the same
#'   order as \code{y}) and one column per random-effect coefficient, and
#'   \code{randeff_variance_component_design} must have one row per column of
#'   \code{randeff_design}. Note that the random effects are sampled but do
#'   \emph{not} enter the count likelihood (the fitted log-mean is
#'   \code{offset + mu_f + z * tau_f} only), so they are carried for
#'   compatibility with the Gaussian BCF samplers rather than being usable
#'   group-level effects; a warning is issued when they are enabled. Put group
#'   indicators in \code{x_control} to model group structure.
#' @param return_trees If TRUE, the serialized \code{tree_samples} for each
#'   forest are returned, enabling counterfactual prediction at new
#'   \code{X} or \code{z} via \code{\link{get_forest_fit}}. Defaults to FALSE
#'   because the tree trace is memory-heavy.
#' @param debug If TRUE, returns \code{const_args} and the raw C++ fit object
#'   for debugging.
#' @param x_control_predict,x_moderate_predict Optional covariate matrices at
#'   which to evaluate the count-component prognostic and moderating forests
#'   during posterior sampling. Both matrices must have the same positive
#'   number of rows and match the columns of their training matrices. These
#'   predictions require no stored tree trace and do not affect model fitting.
#'   The zero-inflation forests, if present, are not predicted by these arguments.
#' @param pihat_predict Propensity estimates for prediction rows, required when
#'   the fitted model includes the propensity among either forest's covariates.
#'
#' @return An object of class \code{"countbcf_fit"}: a list with posterior
#' draws of the log mean and the per-forest fits. Notable elements:
#' \describe{
#'   \item{\code{yhat_log}, \code{yhat}}{\eqn{n \times n_{\mathrm{sim}}} posterior
#'     draws of \eqn{\log E[Y_i \mid X_i, Z_i = z_i^{\mathrm{obs}}]} and its exponential.}
#'   \item{\code{order_vec}}{Permutation used to sort \code{y} (zeros first).
#'     Rows of \code{yhat*} follow this order.}
#'   \item{\code{mu_f_post}, \code{tau_f_post}}{
#'     \eqn{n_{\mathrm{sim}} \times n} per-iter raw coefficients
#'     of the count component's prognostic / moderating forest in the
#'     \emph{original} input order. \code{tau_f_post} stores
#'     \eqn{\tau_f(x_i)} (not multiplied by \eqn{z_i}).}
#'   \item{\code{mu_f0_post}, \code{tau_f0_post}, \code{mu_f1_post}, \code{tau_f1_post}}{
#'     Analogous quantities for the ZI log-odds components (ZI models only).}
#'   \item{\code{kappa}, \code{kappa_acpt}}{NB dispersion posterior and M-H
#'     acceptance rate (NB / ZINB only).}
#'   \item{\code{control_fit$tree_samples}, \code{moderate_fit$tree_samples}, ...}{
#'     Serialized tree traces for each of the 2 (non-ZI) or 6 (ZI)
#'     forests, when \code{return_trees = TRUE}.}
#' }
#'
#' @section Recovering CATE and ATE:
#' For a Poisson/NB model the per-iter, per-unit CATE on the log-rate scale
#' is \code{fit$tau_f_post} and the CATE on the rate scale is
#' \code{exp(fit$mu_f_post + fit$tau_f_post) - exp(fit$mu_f_post)}.
#'
#' For a ZIP/ZINB model the response-scale CATE is obtained by combining the
#' two potential outcomes:
#' \preformatted{
#'   sigmoid <- function(z) 1 / (1 + exp(-z))
#'   log_lambda_0 <- fit$mu_f_post
#'   log_lambda_1 <- fit$mu_f_post + fit$tau_f_post
#'   zi_logit_0   <- fit$mu_f0_post                   - fit$mu_f1_post
#'   zi_logit_1   <- (fit$mu_f0_post + fit$tau_f0_post) -
#'                   (fit$mu_f1_post + fit$tau_f1_post)
#'   mu0 <- (1 - sigmoid(zi_logit_0)) * exp(log_lambda_0)
#'   mu1 <- (1 - sigmoid(zi_logit_1)) * exp(log_lambda_1)
#'   cate <- mu1 - mu0
#'   ate  <- rowMeans(cate)   # one ATE draw per MCMC iteration
#' }
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' n  <- 500
#' x  <- matrix(rnorm(n * 3), n, 3)
#' pi <- plogis(0.5 * x[, 1] - 0.4 * x[, 2])
#' z  <- rbinom(n, 1, pi)
#' log_lambda <- 1 + 0.5 * x[, 1] - 0.3 * x[, 2] + z * (0.3 + 0.2 * x[, 1])
#' p_zi       <- plogis(-1 + 0.5 * x[, 2] + z * (-0.3 + 0.1 * x[, 3]))
#' y <- ifelse(rbinom(n, 1, p_zi) == 1, 0L, rpois(n, exp(log_lambda)))
#'
#' fit <- countbcf(
#'   y = y, z = z, x_control = x, pihat = pi,
#'   nburn = 500, nsim = 500,
#'   count_model = "zipoisson"
#' )
#' }
#'
#' @seealso \code{\link{count_bart}} for the non-causal log-linear BART
#'   count regression that countbcf extends.
#'
#' @export
#'
countbcf <- function(
  y,
  z,
  x_control,
  x_moderate = x_control,
  x_zero = x_control,
  x_pos = x_control,
  pihat = rep(0.5, length(y)),
  offset = NULL,
  nburn,
  nsim,
  nthin = 1,
  update_interval = 100,
  ntree_control = 250,
  ntree_moderate = 50,
  nztree_control = 100,
  nztree_moderate = 50,
  a0 = NA,
  z_conc = 3.5 / sqrt(2),
  base_control = 0.95,
  power_control = 2,
  base_moderate = 0.25,
  power_moderate = 3,
  kappa_a = 5,
  kappa_b = 3,
  kappa_prop_sd = 0.21,
  count_model = "poisson",
  include_pihat = "control",
  randeff_design = matrix(1),
  randeff_variance_component_design = matrix(1),
  randeff_scales = 1,
  randeff_df = 3,
  return_trees = FALSE,
  debug = FALSE,
  x_control_predict = NULL,
  x_moderate_predict = x_control_predict,
  pihat_predict = NULL
) {

  if (!(count_model %in% c("poisson", "nb", "zipoisson", "zinb"))) {
    stop("Undefined 'count_model'; must be 'poisson', 'nb', 'zipoisson', or 'zinb'")
  }
  model_type <- which(count_model == c("poisson", "nb", "zipoisson", "zinb"))
  is_zi <- (model_type == 3) || (model_type == 4)

  if (any(is.na(y))) stop("Missing values in y")
  if (any(is.na(z))) stop("Missing values in z")
  if (any(is.na(x_control))) stop("Missing values in x_control")
  if (any(is.na(x_moderate))) stop("Missing values in x_moderate")
  if (any(is.na(pihat))) stop("Missing values in pihat")

  if (any(!is.finite(y))) stop("Non-numeric values in y")
  if (any(!is.finite(z))) stop("Non-numeric values in z")
  if (any(!is.finite(x_control))) stop("Non-numeric values in x_control")
  if (any(!is.finite(x_moderate))) stop("Non-numeric values in x_moderate")
  if (any(!is.finite(pihat))) stop("Non-numeric values in pihat")

  if (any(y < 0)) stop("Negative values in y")
  if (!isTRUE(all(y == floor(y)))) stop("y must only contain integer values")

  if (!all(sort(unique(z)) %in% c(0, 1))) stop("z must contain only 0 and 1")

  if (!.ident(length(y), nrow(x_control), nrow(x_moderate), length(z), length(pihat))) {
    stop("Data size mismatch among y, z, x_control, x_moderate, pihat")
  }

  if (is_zi) {
    if (any(is.na(x_zero))) stop("Missing values in x_zero")
    if (any(is.na(x_pos))) stop("Missing values in x_pos")
    if (any(!is.finite(x_zero))) stop("Non-numeric values in x_zero")
    if (any(!is.finite(x_pos))) stop("Non-numeric values in x_pos")
    if (!.ident(length(y), nrow(x_zero), nrow(x_pos))) {
      stop("Data size mismatch with x_zero / x_pos")
    }
  }

  if (nburn < 0) stop("nburn must be positive")
  if (nsim < 0) stop("nsim must be positive")
  if (nthin < 0) stop("nthin must be positive")
  if (nthin > nsim + 1) stop("nthin must be < nsim")
  if (nburn < 100) warning("A low (<100) value for nburn was supplied")

  pihat <- as.matrix(pihat)
  if (!(include_pihat %in% c("control", "moderate", "both", "none"))) {
    stop("include_pihat must be 'control', 'moderate', 'both', or 'none'")
  }
  if (include_pihat != "none" && length(unique(pihat)) == 1) {
    warning("All values of pihat are equal. pihat will not be included among covariates")
    include_pihat <- "none"
  }

  predict_count <- !is.null(x_control_predict) || !is.null(x_moderate_predict)
  if (predict_count) {
    if (is.null(x_control_predict) || is.null(x_moderate_predict))
      stop("Both count-component prediction matrices must be supplied")
    X_c_predict <- as.matrix(x_control_predict)
    X_m_predict <- as.matrix(x_moderate_predict)
    if (!is.numeric(X_c_predict) || !is.numeric(X_m_predict) ||
        any(!is.finite(X_c_predict)) || any(!is.finite(X_m_predict)))
      stop("Prediction covariates must be finite numeric values")
    if (nrow(X_c_predict) < 1L || nrow(X_c_predict) != nrow(X_m_predict) ||
        ncol(X_c_predict) != ncol(x_control) ||
        ncol(X_m_predict) != ncol(x_moderate))
      stop("Prediction matrices must match training columns and have equal positive row counts")
    if (include_pihat != "none") {
      if (!is.numeric(pihat_predict) || length(pihat_predict) != nrow(X_c_predict) ||
          any(!is.finite(pihat_predict)))
        stop("Finite pihat_predict values are required for each prediction row")
      if (include_pihat %in% c("control", "both"))
        X_c_predict <- cbind(X_c_predict, pihat_predict)
      if (include_pihat %in% c("moderate", "both"))
        X_m_predict <- cbind(X_m_predict, pihat_predict)
    }
  }

  ### offset
  if (is.null(offset)) {
    mu0 <- rep(0, length(y))
  } else if (length(offset) == 1) {
    mu0 <- rep(log(offset), length(y))
  } else if (length(offset) == length(y)) {
    mu0 <- log(offset)
  } else {
    stop("offset incorrectly specified")
  }

  ### sort by y so the zeros come first (required by countbcf C++ code)
  order_vec <- order(y)
  y_sort <- y[order_vec]
  z_sort <- z[order_vec]
  mu0_sort <- mu0[order_vec]
  pihat_sort <- pihat[order_vec, , drop = FALSE]

  ### the random-effects design is indexed by observation, so its rows must
  ### follow the same permutation as y/z/X/offset
  randeff_variance_component_design <- as.matrix(randeff_variance_component_design)
  randeff_design_sort <- .sort_randeff_design(
    randeff_design, randeff_variance_component_design, order_vec
  )

  ### the C++ sampler draws the random effects from a Gaussian working model on
  ### y - allfit, but never adds allfit_random into log_mean (see the log_mean
  ### updates in src/countbcf.cpp), so they do not enter the count likelihood
  if (length(randeff_variance_component_design) > 1L) {
    warning(paste0(
      "Random effects are sampled but do not enter the count likelihood: the ",
      "fitted log-mean is offset + mu_f + z * tau_f only. 'random_effects' is ",
      "therefore a Gaussian working-model artefact, and group structure is ",
      "absorbed by the mu forest instead. To model groups, pass group ",
      "indicators in x_control."
    ), call. = FALSE)
  }

  X_c <- matrix(x_control, ncol = ncol(x_control))[order_vec, , drop = FALSE]
  X_m <- matrix(x_moderate, ncol = ncol(x_moderate))[order_vec, , drop = FALSE]
  if (is_zi) {
    X_z_c <- matrix(x_zero, ncol = ncol(x_zero))[order_vec, , drop = FALSE]
    X_z_m <- matrix(x_zero, ncol = ncol(x_zero))[order_vec, , drop = FALSE]
    X_p_c <- matrix(x_pos, ncol = ncol(x_pos))[order_vec, , drop = FALSE]
    X_p_m <- matrix(x_pos, ncol = ncol(x_pos))[order_vec, , drop = FALSE]
  }

  if (include_pihat %in% c("control", "both")) {
    X_c <- cbind(X_c, pihat_sort)
    if (is_zi) {
      X_z_c <- cbind(X_z_c, pihat_sort)
      X_p_c <- cbind(X_p_c, pihat_sort)
    }
  }
  if (include_pihat %in% c("moderate", "both")) {
    X_m <- cbind(X_m, pihat_sort)
    if (is_zi) {
      X_z_m <- cbind(X_z_m, pihat_sort)
      X_p_m <- cbind(X_p_m, pihat_sort)
    }
  }

  n <- length(y_sort)

  ### leaf hyperparameters

  # count-component mu prior (a0)
  if (is.na(a0)) {
    y_star <- quantile(y_sort, probs = 0.95)
    a0 <- 0.5 * (log(max(y_star, 1)) - mean(mu0_sort))
    if (!is.finite(a0) || a0 <= 0) a0 <- 1
  }
  # The mu and tau forests within a component share this leaf concentration.
  # Stronger shrinkage of the tau forests is delivered by their smaller tree
  # count and stub-biased tree-shape prior, not by a separate concentration:
  # with the leaf scale calibrated for ntree_control trees but applied to the
  # ntree_moderate-tree tau forest, the implied tau-function prior SD is
  # a0 * sqrt(ntree_moderate / ntree_control).
  log_leaf_var_mu <- a0^2 / ntree_control
  c_start_mu <- ntree_control / a0^2 + 0.5
  leaf_c_mu <- optim(c_start_mu, fn = function(x) (trigamma(x) - log_leaf_var_mu)^2,
                     method = "Brent", lower = 0, upper = 1e9)$par
  leaf_d_mu <- exp(digamma(leaf_c_mu))

  if (is_zi) {
    z_leaf_var_mu <- z_conc^2 / nztree_control
    z_c_start_mu <- nztree_control / z_conc^2 + 0.5
    z_c_mu <- optim(z_c_start_mu, fn = function(x) (trigamma(x) - z_leaf_var_mu)^2,
                    method = "Brent", lower = 0, upper = 1e9)$par
    z_d_mu <- exp(digamma(z_c_mu))
  } else {
    z_c_mu <- 0; z_d_mu <- 0
  }

  ### leaf-scale (Sigma0) is unused in the loglinear leaf draw but required by the spec.
  Sigma0_dummy <- matrix(1, nrow = 1, ncol = 1)

  ### designs and specs
  X_list <- list(X_c, X_m)
  basis_list <- list(matrix(rep(1, n), ncol = 1), matrix(z_sort, ncol = 1))

  if (is_zi) {
    X_list <- c(X_list, list(X_z_c, X_z_m, X_p_c, X_p_m))
    basis_list <- c(basis_list,
                    list(matrix(rep(1, n), ncol = 1), matrix(z_sort, ncol = 1),
                         matrix(rep(1, n), ncol = 1), matrix(z_sort, ncol = 1)))
  }

  designs <- countbcf:::make_bart_designs(X_list, basis_list)
  if (predict_count) {
    # These rows are evaluated using the training cutpoints and current trees;
    # they never enter the likelihood, tree proposals or propensity fitting.
    designs[[1]]$X_predict <- t(X_c_predict)
    designs[[2]]$X_predict <- t(X_m_predict)
  }

  make_spec <- function(design_idx, ntree, alpha, beta, vanilla, function_group) {
    design <- designs[[design_idx]]
    spec <- countbcf:::make_bart_spec(
      design = design,
      ntree = ntree,
      Sigma0 = Sigma0_dummy,
      scale_df = -1,
      vanilla = vanilla,
      alpha = alpha,
      beta = beta,
      update_leaf_scale = FALSE
    )
    spec$function_group <- function_group
    spec
  }

  specs <- list(
    make_spec(1, ntree_control,  base_control,  power_control,  TRUE,  0),  # mu_f
    make_spec(2, ntree_moderate, base_moderate, power_moderate, FALSE, 0)   # tau_f
  )
  if (is_zi) {
    specs <- c(specs, list(
      make_spec(3, nztree_control,  base_control,  power_control,  TRUE,  1),  # mu_f0
      make_spec(4, nztree_moderate, base_moderate, power_moderate, FALSE, 1),  # tau_f0
      make_spec(5, nztree_control,  base_control,  power_control,  TRUE,  2),  # mu_f1
      make_spec(6, nztree_moderate, base_moderate, power_moderate, FALSE, 2)   # tau_f1
    ))
  }

  ### The C++ sampler picks the leaf-prior (c, d) purely from function_group:
  ### (leaf_c, leaf_d) for group 0 (count) and (z_c, z_d) for groups 1/2 (ZI).
  ### The mu and tau forests within a group therefore share the same (c, d);
  ### the mu/tau scale separation is delivered by the tree-count difference
  ### and the stub-biased tree-shape prior, not by a per-forest concentration.

  ### ghost lambda/nu (unused for count models)
  lambda <- 1
  nu <- 3

  ### The C++ countbcf selects (c, d) from these function args using
  ### function_group; mu and tau within a group share the same (c, d), and the
  ### leaf-scale separation between mu and tau is handled by the ntree
  ### difference and the tree-shape prior.
  const_args <- list(
    y_     = y_sort,
    offset_ = mu0_sort,
    bart_specs = specs,
    bart_designs = designs,
    random_des = randeff_design_sort,
    random_var_ix = randeff_variance_component_design,
    random_var = matrix(rep(1e-8, ncol(randeff_variance_component_design)), ncol = 1),
    random_var_df = randeff_df,
    randeff_scales = randeff_scales,
    burn = nburn, nd = nsim, thin = nthin,
    count_model = model_type,
    lambda = lambda, nu = nu,
    kappa_a = kappa_a, kappa_b = kappa_b,
    leaf_c = leaf_c_mu, leaf_d = leaf_d_mu,
    z_c = z_c_mu, z_d = z_d_mu,
    kappa_prop_sd = kappa_prop_sd,
    status_interval = update_interval,
    text_trace = TRUE
  )

  ### NOTE: the Rcpp-generated R binding for the C++ countbcf is masked in
  ### this package's namespace by the R wrapper (this function) of the same
  ### name. Calling `countbcf:::countbcf` resolves to the wrapper itself,
  ### producing an "argument 20 matches multiple formal arguments" error
  ### because `z_c` partial-matches `z_conc` in the wrapper's
  ### formal-argument list. Bypass the namespace lookup by invoking the
  ### compiled function directly through .Call.
  fit_raw <- with(const_args, .Call(
    "_countbcf_countbcf",
    PACKAGE = "countbcf",
    y_, offset_, bart_specs, bart_designs,
    random_des, random_var, random_var_ix, random_var_df, randeff_scales,
    burn, nd, thin,
    count_model, lambda, nu,
    kappa_a, kappa_b, leaf_c, leaf_d, z_c, z_d,
    kappa_prop_sd,
    as.logical(return_trees),  # return_trees
    FALSE,                # save_trees
    FALSE,                # est_mod_fits
    FALSE,                # est_con_fits
    FALSE,                # prior_sample
    as.integer(status_interval),
    as.numeric(c(0.0)),   # lower_bd
    as.numeric(c(0.0)),   # upper_bd
    FALSE,                # probit
    text_trace,
    FALSE                 # R_trace
  ))

  ### organize per-forest fits into mu/tau components
  forest_fits <- fit_raw$forest_fits
  fg <- fit_raw$function_group

  out <- list(
    yhat_log = fit_raw$yhat_post,
    yhat = exp(fit_raw$yhat_post),
    order_vec = order_vec,
    sigma = fit_raw$sigma,
    random_effects = fit_raw$gamma,
    random_effects_sd = fit_raw$random_sd_post
  )

  ### per-forest raw coefficients at the training X.
  ### fit_raw$coefs[[s]] is a (basis_dim=1, n_sorted, nsim) arma::cube. How
  ### RcppArmadillo materializes it on the R side is version-dependent (3D
  ### array, 2D matrix with the singleton dropped, or a flat numeric vector),
  ### but the underlying memory layout is column-major (slice j contains
  ### n_sorted consecutive doubles), so we reshape from the flat values.
  inv <- order(order_vec)
  n_sorted <- length(order_vec)
  .extract_coef <- function(cc) {
    nsim <- length(cc) / n_sorted
    arr <- matrix(as.numeric(cc), nrow = n_sorted, ncol = nsim)
    t(arr[inv, , drop = FALSE])         # nsim x n (original order)
  }

  ### group 0 (count): forests 1 (mu) and 2 (tau)
  out$mu_f_log  <- forest_fits[[1]]
  out$tau_f_log <- forest_fits[[2]]   # already = z_trt_i * tau_f(x_i) per obs
  out$mu_f_post  <- .extract_coef(fit_raw$coefs[[1]])
  out$tau_f_post <- .extract_coef(fit_raw$coefs[[2]])
  if (predict_count) {
    np <- nrow(X_c_predict)
    prediction_matrix <- function(cc) t(matrix(as.numeric(cc), nrow = np))
    out$mu_f_predict <- prediction_matrix(fit_raw$prediction_coefs[[1]])
    out$tau_f_predict <- prediction_matrix(fit_raw$prediction_coefs[[2]])
  }

  if (is_zi) {
    out$mu_f0_log  <- forest_fits[[3]]
    out$tau_f0_log <- forest_fits[[4]]
    out$mu_f1_log  <- forest_fits[[5]]
    out$tau_f1_log <- forest_fits[[6]]
    out$mu_f0_post  <- .extract_coef(fit_raw$coefs[[3]])
    out$tau_f0_post <- .extract_coef(fit_raw$coefs[[4]])
    out$mu_f1_post  <- .extract_coef(fit_raw$coefs[[5]])
    out$tau_f1_post <- .extract_coef(fit_raw$coefs[[6]])
  }

  if (return_trees) {
    out$control_fit  <- list(tree_samples = fit_raw$tree_trace[[1]])
    out$moderate_fit <- list(tree_samples = fit_raw$tree_trace[[2]])
    if (is_zi) {
      out$f0_control_fit  <- list(tree_samples = fit_raw$tree_trace[[3]])
      out$f0_moderate_fit <- list(tree_samples = fit_raw$tree_trace[[4]])
      out$f1_control_fit  <- list(tree_samples = fit_raw$tree_trace[[5]])
      out$f1_moderate_fit <- list(tree_samples = fit_raw$tree_trace[[6]])
    }
  }

  if (model_type %in% c(2, 4)) {
    out$kappa <- fit_raw$kappa
    out$kappa_acpt <- fit_raw$kappa_acceptance
    out$kappa_a <- kappa_a
    out$kappa_b <- kappa_b
  }

  if (debug) {
    out$const_args <- const_args
    out$raw_fit <- fit_raw
  }

  class(out) <- "countbcf_fit"
  return(out)
}
