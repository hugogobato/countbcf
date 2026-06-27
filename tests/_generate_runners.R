################################################################################
##  Generates the full set of self-contained Colab runners.
##
##  Each emitted file contains everything needed to reproduce ONE cell of the
##  simulation grid for ONE (study, method, DGP) combination, for N_SIM=100
##  Monte Carlo replicates.  The only external dependency is `library(countbart)`
##  (and base R).  No source() of any repo file.
##
##  Design grid (CATE-focused DGPs only):
##    Reference cell ......... N=250, P=5,  ATE=1.25
##    Ablation: N ............ N=100, P=5,  ATE=1.25   (N=500 added separately)
##                              N=500, P=5,  ATE=1.25
##    Ablation: P ............ N=250, P=50,  ATE=1.25
##                              N=250, P=250, ATE=1.25
##    Ablation: ATE .......... N=250, P=5,  ATE=0.5
##                              N=250, P=5,  ATE=2.5
##    (One-at-a-time; no combinatorial cells.)
##
##  Output:
##    runners_count/run_count_<method>_<dgp>_N<N>_P<P>_ATE<ATEx100>.R
##    runners_zi/run_zi_<method>_<dgp>_N<N>_P<P>_ATE<ATEx100>.R
##
##  Usage:
##    Rscript tests/_generate_runners.R
################################################################################

OUT_DIR_COUNT <- "runners_count"
OUT_DIR_ZI    <- "runners_zi"
dir.create(OUT_DIR_COUNT, showWarnings = FALSE, recursive = TRUE)
dir.create(OUT_DIR_ZI,    showWarnings = FALSE, recursive = TRUE)

## ----------------------------------------------------------------------------
##  Cells: 7 per (study, method, DGP).
## ----------------------------------------------------------------------------

CELLS <- data.frame(
  factor_varied = c("REF", "N",   "N",   "P",   "P",   "ATE", "ATE"),
  N             = c(250L,  100L,  500L,  250L,  250L,  250L,  250L),
  P             = c(5L,    5L,    5L,    50L,   250L,  5L,    5L),
  ATE           = c(1.25,  1.25,  1.25,  1.25,  1.25,  0.50,  2.50),
  stringsAsFactors = FALSE
)

## ----------------------------------------------------------------------------
##  DGP-specific definitions (CATE-focused only).
## ----------------------------------------------------------------------------

F_MU_LINEAR        <- "1.0 + 0.5 * X[, 1] - 0.3 * X[, 2] + 0.2 * X[, 3]"
F_MU_NONLINEAR     <- "0.8 + 0.6 * sin(X[, 1]) + 0.4 * X[, 2] * X[, 3] + 0.5 * (X[, 4] > 0)"
F_TAU_LIN_CATE     <- "0.60 * X[, 1] - 0.45 * X[, 2] + 0.30 * X[, 3]"
F_TAU_NONLIN_CATE  <- "0.80 * sin(X[, 1]) + 0.50 * (X[, 2] > 0) * X[, 3] - 0.40 * X[, 5]"

F_MU_ZI_LIN        <- "-1.0 + 0.5 * X[, 2]"
F_MU_ZI_NONLIN     <- "-1.0 + 0.7 * abs(X[, 2]) - 0.5 * exp(-X[, 3]^2)"
F_TAU_ZI_LIN_CATE  <- "-0.50 * X[, 1] + 0.30 * X[, 3] - 0.40 * X[, 4]"
F_TAU_ZI_NL_CATE   <- "-0.60 * X[, 1] * X[, 4] + 0.40 * X[, 2] * (X[, 5] > 0)"

DGPS_COUNT <- list(
  linear_poisson_cate    = list(f_mu = F_MU_LINEAR,    f_tau = F_TAU_LIN_CATE,    kappa = "NA_real_", count_model = "poisson"),
  nonlinear_poisson_cate = list(f_mu = F_MU_NONLINEAR, f_tau = F_TAU_NONLIN_CATE, kappa = "NA_real_", count_model = "poisson"),
  linear_nb_cate         = list(f_mu = F_MU_LINEAR,    f_tau = F_TAU_LIN_CATE,    kappa = "2",        count_model = "nb"),
  nonlinear_nb_cate      = list(f_mu = F_MU_NONLINEAR, f_tau = F_TAU_NONLIN_CATE, kappa = "2",        count_model = "nb")
)

DGPS_ZI <- list(
  linear_zip_cate     = list(f_mu = F_MU_LINEAR,    f_tau = F_TAU_LIN_CATE,    f_mu_zi = F_MU_ZI_LIN,    f_tau_zi = F_TAU_ZI_LIN_CATE, kappa = "NA_real_", count_model = "zipoisson"),
  nonlinear_zip_cate  = list(f_mu = F_MU_NONLINEAR, f_tau = F_TAU_NONLIN_CATE, f_mu_zi = F_MU_ZI_NONLIN, f_tau_zi = F_TAU_ZI_NL_CATE,  kappa = "NA_real_", count_model = "zipoisson"),
  linear_zinb_cate    = list(f_mu = F_MU_LINEAR,    f_tau = F_TAU_LIN_CATE,    f_mu_zi = F_MU_ZI_LIN,    f_tau_zi = F_TAU_ZI_LIN_CATE, kappa = "2",        count_model = "zinb"),
  nonlinear_zinb_cate = list(f_mu = F_MU_NONLINEAR, f_tau = F_TAU_NONLIN_CATE, f_mu_zi = F_MU_ZI_NONLIN, f_tau_zi = F_TAU_ZI_NL_CATE,  kappa = "2",        count_model = "zinb")
)

## ----------------------------------------------------------------------------
##  Method-specific fit bodies.
## ----------------------------------------------------------------------------

FIT_COUNTBCF_COUNT <- '
fit_method <- function(d) {
  fit <- countbcf(
    y = d$y, z = d$z,
    x_control = d$x, x_moderate = d$x,
    pihat = d$pihat,
    nburn = NBURN, nsim = NSIM, nthin = NTHIN,
    count_model = COUNT_MODEL,
    include_pihat = "control",
    update_interval = max(1L, (NBURN + NSIM) %/% 5),
    return_trees = FALSE)
  mu0_post  <- exp(fit$mu_f_post)
  mu1_post  <- exp(fit$mu_f_post + fit$tau_f_post)
  cate_post <- mu1_post - mu0_post
  z_mat <- matrix(d$z, nrow = nrow(mu0_post),
                  ncol = length(d$z), byrow = TRUE)
  yhat  <- (1 - z_mat) * mu0_post + z_mat * mu1_post
  list(ate = rowMeans(cate_post),
       cate_post = cate_post,
       yhat_mean = colMeans(yhat))
}
'

FIT_BCF_GAUSS_COUNT <- '
fit_method <- function(d) {
  fit <- bcf_binary(
    y = d$y, z = d$z,
    x_control = d$x, x_moderate = d$x,
    pihat = d$pihat,
    nburn = NBURN, nsim = NSIM, nthin = NTHIN,
    include_pi = "control",
    sighat = sd(d$y),
    update_interval = max(1L, (NBURN + NSIM) %/% 5))
  tau_post <- get_forest_fit(fit$moderate_fit, d$x)
  list(ate = rowMeans(tau_post),
       cate_post = tau_post,
       yhat_mean = colMeans(fit$yhat))
}
'

FIT_COUNTBCF_ZI <- '
fit_method <- function(d) {
  fit <- countbcf(
    y = d$y, z = d$z,
    x_control = d$x, x_moderate = d$x,
    x_zero = d$x, x_pos = d$x,
    pihat = d$pihat,
    nburn = NBURN, nsim = NSIM, nthin = NTHIN,
    count_model = COUNT_MODEL,
    include_pihat = "control",
    update_interval = max(1L, (NBURN + NSIM) %/% 5),
    return_trees = FALSE)
  log_lambda_0 <- fit$mu_f_post
  log_lambda_1 <- fit$mu_f_post + fit$tau_f_post
  zi_logit_0   <- fit$mu_f0_post              - fit$mu_f1_post
  zi_logit_1   <- (fit$mu_f0_post + fit$tau_f0_post) -
                   (fit$mu_f1_post + fit$tau_f1_post)
  mu0_post <- (1 - sigmoid(zi_logit_0)) * exp(log_lambda_0)
  mu1_post <- (1 - sigmoid(zi_logit_1)) * exp(log_lambda_1)
  cate_post <- mu1_post - mu0_post
  z_mat <- matrix(d$z, nrow = nrow(mu0_post),
                  ncol = length(d$z), byrow = TRUE)
  yhat <- (1 - z_mat) * mu0_post + z_mat * mu1_post
  list(ate = rowMeans(cate_post),
       cate_post = cate_post,
       yhat_mean = colMeans(yhat))
}
'

FIT_BCF_GAUSS_ZI <- FIT_BCF_GAUSS_COUNT  # identical signature for zi study

## ----------------------------------------------------------------------------
##  Parametric / classical benchmark fit bodies.
##
##  These produce the same (ate, cate_post, yhat_mean) contract expected by
##  summarize_fit().  "Posterior" draws for the parametric GLMs are obtained by
##  simulating the regression coefficients from their asymptotic sampling
##  distribution N(beta_hat, vcov) -- the King-Tomz-Wittenberg / arm::sim
##  approach -- and propagating them to the per-unit CATE
##  CATE(x) = E[Y | z = 1, x] - E[Y | z = 0, x].
## ----------------------------------------------------------------------------

##  Parametric Poisson GLM:  Y ~ Z * X with log link.
FIT_GLM_POISSON <- '
fit_method <- function(d) {
  p  <- ncol(d$x)
  xn <- paste0("X", seq_len(p))
  df <- data.frame(y = d$y, z = d$z, d$x)
  colnames(df) <- c("y", "z", xn)
  rhs  <- paste0("z * (", paste(xn, collapse = " + "), ")")
  form <- stats::as.formula(paste("y ~", rhs))
  fit  <- stats::glm(form, family = stats::poisson(link = "log"), data = df)
  beta <- stats::coef(fit)
  V    <- stats::vcov(fit)
  if (any(!is.finite(beta)) || any(!is.finite(V)))
    stop("glm_poisson: non-finite (aliased) coefficients; p >= n?")
  df0 <- df; df0$z <- 0
  df1 <- df; df1$z <- 1
  X0  <- stats::model.matrix(form, df0)
  X1  <- stats::model.matrix(form, df1)
  draws <- MASS::mvrnorm(n = NSIM, mu = beta, Sigma = V)
  if (is.null(dim(draws))) draws <- matrix(draws, nrow = NSIM)
  mu0_post  <- exp(draws %*% t(X0))
  mu1_post  <- exp(draws %*% t(X1))
  cate_post <- mu1_post - mu0_post
  z_mat <- matrix(d$z, nrow = NSIM, ncol = length(d$z), byrow = TRUE)
  yhat  <- (1 - z_mat) * mu0_post + z_mat * mu1_post
  list(ate = rowMeans(cate_post),
       cate_post = cate_post,
       yhat_mean = colMeans(yhat))
}
'

##  Parametric Negative-Binomial GLM:  Y ~ Z * X (MASS::glm.nb), log link.
FIT_GLM_NB <- '
fit_method <- function(d) {
  p  <- ncol(d$x)
  xn <- paste0("X", seq_len(p))
  df <- data.frame(y = d$y, z = d$z, d$x)
  colnames(df) <- c("y", "z", xn)
  rhs  <- paste0("z * (", paste(xn, collapse = " + "), ")")
  form <- stats::as.formula(paste("y ~", rhs))
  fit  <- suppressWarnings(MASS::glm.nb(form, data = df))
  beta <- stats::coef(fit)
  V    <- stats::vcov(fit)
  if (any(!is.finite(beta)) || any(!is.finite(V)))
    stop("glm_nb: non-finite (aliased) coefficients; p >= n?")
  df0 <- df; df0$z <- 0
  df1 <- df; df1$z <- 1
  X0  <- stats::model.matrix(form, df0)
  X1  <- stats::model.matrix(form, df1)
  draws <- MASS::mvrnorm(n = NSIM, mu = beta, Sigma = V)
  if (is.null(dim(draws))) draws <- matrix(draws, nrow = NSIM)
  mu0_post  <- exp(draws %*% t(X0))
  mu1_post  <- exp(draws %*% t(X1))
  cate_post <- mu1_post - mu0_post
  z_mat <- matrix(d$z, nrow = NSIM, ncol = length(d$z), byrow = TRUE)
  yhat  <- (1 - z_mat) * mu0_post + z_mat * mu1_post
  list(ate = rowMeans(cate_post),
       cate_post = cate_post,
       yhat_mean = colMeans(yhat))
}
'

##  Zero-Inflated Poisson (pscl::zeroinfl, dist = "poisson").
##  E[Y | x] = (1 - pi(x)) * lambda(x); both parts share the Z * X design.
FIT_GLM_ZIP <- '
fit_method <- function(d) {
  p  <- ncol(d$x)
  xn <- paste0("X", seq_len(p))
  df <- data.frame(y = d$y, z = d$z, d$x)
  colnames(df) <- c("y", "z", xn)
  rhs  <- paste0("z * (", paste(xn, collapse = " + "), ")")
  form <- stats::as.formula(paste("y ~", rhs, "|", rhs))
  fit  <- pscl::zeroinfl(form, data = df, dist = "poisson")
  bc <- fit$coefficients$count
  bz <- fit$coefficients$zero
  pc <- length(bc); pz <- length(bz)
  b  <- c(bc, bz)
  V  <- stats::vcov(fit)[seq_len(pc + pz), seq_len(pc + pz), drop = FALSE]
  if (any(!is.finite(b)) || any(!is.finite(V)))
    stop("glm_zip: non-finite coefficients; p >= n?")
  df0 <- df; df0$z <- 0
  df1 <- df; df1$z <- 1
  rhs_form <- stats::as.formula(paste("~", rhs))
  X0 <- stats::model.matrix(rhs_form, df0)
  X1 <- stats::model.matrix(rhs_form, df1)
  if (ncol(X0) != pc) stop("glm_zip: design/coef mismatch")
  draws <- MASS::mvrnorm(n = NSIM, mu = b, Sigma = V)
  if (is.null(dim(draws))) draws <- matrix(draws, nrow = NSIM)
  dc <- draws[, seq_len(pc), drop = FALSE]
  dz <- draws[, pc + seq_len(pz), drop = FALSE]
  lam0 <- exp(dc %*% t(X0)); lam1 <- exp(dc %*% t(X1))
  pi0  <- stats::plogis(dz %*% t(X0)); pi1 <- stats::plogis(dz %*% t(X1))
  mu0_post  <- (1 - pi0) * lam0
  mu1_post  <- (1 - pi1) * lam1
  cate_post <- mu1_post - mu0_post
  z_mat <- matrix(d$z, nrow = NSIM, ncol = length(d$z), byrow = TRUE)
  yhat  <- (1 - z_mat) * mu0_post + z_mat * mu1_post
  list(ate = rowMeans(cate_post),
       cate_post = cate_post,
       yhat_mean = colMeans(yhat))
}
'

##  Zero-Inflated Negative-Binomial (pscl::zeroinfl, dist = "negbin").
##  Identical mean structure to ZIP; vcov has a trailing log(theta) entry that
##  we drop (theta does not enter the mean (1 - pi) * lambda).
FIT_GLM_ZINB <- '
fit_method <- function(d) {
  p  <- ncol(d$x)
  xn <- paste0("X", seq_len(p))
  df <- data.frame(y = d$y, z = d$z, d$x)
  colnames(df) <- c("y", "z", xn)
  rhs  <- paste0("z * (", paste(xn, collapse = " + "), ")")
  form <- stats::as.formula(paste("y ~", rhs, "|", rhs))
  fit  <- pscl::zeroinfl(form, data = df, dist = "negbin")
  bc <- fit$coefficients$count
  bz <- fit$coefficients$zero
  pc <- length(bc); pz <- length(bz)
  b  <- c(bc, bz)
  V  <- stats::vcov(fit)[seq_len(pc + pz), seq_len(pc + pz), drop = FALSE]
  if (any(!is.finite(b)) || any(!is.finite(V)))
    stop("glm_zinb: non-finite coefficients; p >= n?")
  df0 <- df; df0$z <- 0
  df1 <- df; df1$z <- 1
  rhs_form <- stats::as.formula(paste("~", rhs))
  X0 <- stats::model.matrix(rhs_form, df0)
  X1 <- stats::model.matrix(rhs_form, df1)
  if (ncol(X0) != pc) stop("glm_zinb: design/coef mismatch")
  draws <- MASS::mvrnorm(n = NSIM, mu = b, Sigma = V)
  if (is.null(dim(draws))) draws <- matrix(draws, nrow = NSIM)
  dc <- draws[, seq_len(pc), drop = FALSE]
  dz <- draws[, pc + seq_len(pz), drop = FALSE]
  lam0 <- exp(dc %*% t(X0)); lam1 <- exp(dc %*% t(X1))
  pi0  <- stats::plogis(dz %*% t(X0)); pi1 <- stats::plogis(dz %*% t(X1))
  mu0_post  <- (1 - pi0) * lam0
  mu1_post  <- (1 - pi1) * lam1
  cate_post <- mu1_post - mu0_post
  z_mat <- matrix(d$z, nrow = NSIM, ncol = length(d$z), byrow = TRUE)
  yhat  <- (1 - z_mat) * mu0_post + z_mat * mu1_post
  list(ate = rowMeans(cate_post),
       cate_post = cate_post,
       yhat_mean = colMeans(yhat))
}
'

##  Causal Forest (grf::causal_forest) -- heterogeneous treatment effects.
##  Per-unit CATE draws are simulated from N(tau_hat(x), var_hat(x)) using the
##  forest variance estimates; the ATE (and its CI) come from the doubly-robust
##  AIPW estimate average_treatment_effect().  yhat uses the R-learner identity
##  E[Y | x, z] = m(x) + (z - e(x)) * tau(x), with m = Y.hat, e = W.hat.
FIT_GRF <- '
fit_method <- function(d) {
  forest <- grf::causal_forest(X = d$x, Y = d$y, W = d$z,
                               num.trees = 2000L, seed = 1L)
  pred <- predict(forest, estimate.variance = TRUE)
  tau  <- pred$predictions
  v    <- pmax(pred$variance.estimates, 0)
  n    <- length(tau)
  cate_post <- matrix(rnorm(NSIM * n,
                            mean = rep(tau,     each = NSIM),
                            sd   = rep(sqrt(v), each = NSIM)),
                      nrow = NSIM, ncol = n)
  ate_est  <- grf::average_treatment_effect(forest, target.sample = "all")
  ate_post <- rnorm(NSIM, mean = ate_est[["estimate"]],
                          sd   = ate_est[["std.err"]])
  yhat_mean <- as.numeric(forest$Y.hat + (d$z - forest$W.hat) * tau)
  list(ate = ate_post,
       cate_post = cate_post,
       yhat_mean = yhat_mean)
}
'

## ----------------------------------------------------------------------------
##  Per-method package-setup snippets (emitted as the notebook install cell).
##  Only the dependencies actually used by the fit body are installed.
## ----------------------------------------------------------------------------

PKG_COUNTBCF <- 'install.packages("remotes")
if (!require("devtools")) {
  install.packages("devtools")
}
devtools::install_github("hugogobato/countbcf")
library(countbcf)'

PKG_MASS <- 'if (!requireNamespace("MASS", quietly = TRUE)) install.packages("MASS")'

PKG_PSCL <- 'if (!requireNamespace("pscl", quietly = TRUE)) install.packages("pscl")
if (!requireNamespace("MASS", quietly = TRUE)) install.packages("MASS")'

PKG_GRF  <- 'if (!requireNamespace("grf", quietly = TRUE)) install.packages("grf")'

## Each method carries its install snippet (pkg) and its fit body (fit).
METHODS_COUNT <- list(
  countbcf      = list(pkg = PKG_COUNTBCF, fit = FIT_COUNTBCF_COUNT),
  bcf_gauss     = list(pkg = PKG_COUNTBCF, fit = FIT_BCF_GAUSS_COUNT),
  glm_poisson   = list(pkg = PKG_MASS,     fit = FIT_GLM_POISSON),
  glm_nb        = list(pkg = PKG_MASS,     fit = FIT_GLM_NB),
  causal_forest = list(pkg = PKG_GRF,      fit = FIT_GRF)
)
METHODS_ZI <- list(
  countbcf      = list(pkg = PKG_COUNTBCF, fit = FIT_COUNTBCF_ZI),
  bcf_gauss     = list(pkg = PKG_COUNTBCF, fit = FIT_BCF_GAUSS_ZI),
  glm_zip       = list(pkg = PKG_PSCL,     fit = FIT_GLM_ZIP),
  glm_zinb      = list(pkg = PKG_PSCL,     fit = FIT_GLM_ZINB),
  causal_forest = list(pkg = PKG_GRF,      fit = FIT_GRF)
)

## ----------------------------------------------------------------------------
##  Common boilerplate (helpers, summarize_fit, empty_metrics).
## ----------------------------------------------------------------------------

HEADER <- '################################################################################
##  Self-contained Colab runner.
##  Study   : %s
##  Method  : %s
##  DGP     : %s
##  Cell    : N=%d  P=%d  ATE=%.2f  (factor varied: %s)
##  N_SIM   : 100 (seeds 1..100 are independent draws)
##
##  Only external dependency: library(countbart).
##  Outputs %s and %s to the working directory.
################################################################################

%s

N_SIM   <- 100L
NBURN   <- 1000L
NSIM    <- 1000L
NTHIN   <- 1L

METHOD       <- "%s"
DGP_NAME     <- "%s"
FACTOR_VAR   <- "%s"
N            <- %dL
P            <- %dL
ATE_TARGET   <- %.2f
KAPPA        <- %s
COUNT_MODEL  <- "%s"

OUT_CSV  <- "%s"
LOG_FILE <- "%s"
if (file.exists(OUT_CSV))  file.remove(OUT_CSV)
if (file.exists(LOG_FILE)) file.remove(LOG_FILE)

log_msg <- function(fmt, ...) {
  msg <- sprintf(paste0("[%%s] ", fmt),
                 format(Sys.time(), "%%Y-%%m-%%d %%H:%%M:%%S"), ...)
  cat(msg, "\\n")
  cat(msg, "\\n", file = LOG_FILE, append = TRUE)
}

sigmoid <- function(x) 1 / (1 + exp(-x))
.true_propensity <- function(X) sigmoid(0.5 * X[, 1] - 0.4 * X[, 2])

.f_mu_count <- function(X) %s
.f_tau_het  <- function(X) %s
'

HEADER_ZI_EXTRA <- '.f_mu_zi    <- function(X) %s
.f_tau_zi   <- function(X) %s
'

SIMULATE_COUNT <- '
.simulate_count <- function(n, p, seed, tau_const) {
  set.seed(seed)
  X <- matrix(rnorm(n * p), nrow = n, ncol = p)
  X[, p]     <- rbinom(n, 1, 0.4)
  X[, p - 1] <- rbinom(n, 1, 0.6)
  pi_x <- .true_propensity(X); Z <- rbinom(n, 1, pi_x)
  log_lambda_0 <- .f_mu_count(X)
  log_lambda_1 <- log_lambda_0 + tau_const + .f_tau_het(X)
  log_lambda   <- ifelse(Z == 1, log_lambda_1, log_lambda_0)
  lambda <- exp(log_lambda)
  Y <- if (is.na(KAPPA)) rpois(n, lambda)
       else              rnbinom(n, size = KAPPA, mu = lambda)
  mu0 <- exp(log_lambda_0); mu1 <- exp(log_lambda_1); cate <- mu1 - mu0
  list(y = Y, z = Z, x = X, pihat = pi_x,
       mu0 = mu0, mu1 = mu1, cate = cate, ate = mean(cate),
       pct_zero = mean(Y == 0))
}

.calibrate_tau_const <- function(p, ate_target,
                                 n_cal = 5000L, seed_cal = 9999991L) {
  old_seed <- if (exists(".Random.seed", envir = globalenv()))
    get(".Random.seed", envir = globalenv()) else NULL
  on.exit({
    if (!is.null(old_seed))
      assign(".Random.seed", old_seed, envir = globalenv())
  })
  set.seed(seed_cal)
  Xc <- matrix(rnorm(n_cal * p), nrow = n_cal, ncol = p)
  Xc[, p]     <- rbinom(n_cal, 1, 0.4)
  Xc[, p - 1] <- rbinom(n_cal, 1, 0.6)
  mu0_c <- exp(.f_mu_count(Xc)); het_c <- .f_tau_het(Xc)
  diff_fn <- function(tc) mean(mu0_c * (exp(tc + het_c) - 1)) - ate_target
  uniroot(diff_fn, lower = -5, upper = 5,
          extendInt = "yes", tol = 1e-4)$root
}

draw_dgp <- function(n, p, seed, ate_target) {
  tc <- .calibrate_tau_const(p, ate_target)
  .simulate_count(n, p, seed, tc)
}
'

SIMULATE_ZI <- '
.simulate_zi <- function(n, p, seed, tau_const) {
  set.seed(seed)
  X <- matrix(rnorm(n * p), nrow = n, ncol = p)
  X[, p]     <- rbinom(n, 1, 0.4)
  X[, p - 1] <- rbinom(n, 1, 0.6)
  pi_x <- .true_propensity(X); Z <- rbinom(n, 1, pi_x)
  log_lambda_0 <- .f_mu_count(X)
  log_lambda_1 <- log_lambda_0 + tau_const + .f_tau_het(X)
  zi_logit_0   <- .f_mu_zi(X)
  zi_logit_1   <- zi_logit_0 + .f_tau_zi(X)
  p_zi_0 <- sigmoid(zi_logit_0); p_zi_1 <- sigmoid(zi_logit_1)
  log_lambda <- ifelse(Z == 1, log_lambda_1, log_lambda_0)
  p_zi       <- ifelse(Z == 1, p_zi_1,       p_zi_0)
  lambda     <- exp(log_lambda)
  is_struct_zero <- rbinom(n, 1, p_zi)
  Y_count <- if (is.na(KAPPA)) rpois(n, lambda)
             else              rnbinom(n, size = KAPPA, mu = lambda)
  Y <- ifelse(is_struct_zero == 1, 0L, Y_count)
  mu0  <- (1 - p_zi_0) * exp(log_lambda_0)
  mu1  <- (1 - p_zi_1) * exp(log_lambda_1)
  cate <- mu1 - mu0
  list(y = Y, z = Z, x = X, pihat = pi_x,
       mu0 = mu0, mu1 = mu1, cate = cate, ate = mean(cate),
       pct_zero = mean(Y == 0),
       pct_struct_zero = mean(is_struct_zero))
}

.calibrate_tau_const_zi <- function(p, ate_target,
                                    n_cal = 5000L, seed_cal = 9999992L) {
  old_seed <- if (exists(".Random.seed", envir = globalenv()))
    get(".Random.seed", envir = globalenv()) else NULL
  on.exit({
    if (!is.null(old_seed))
      assign(".Random.seed", old_seed, envir = globalenv())
  })
  set.seed(seed_cal)
  Xc <- matrix(rnorm(n_cal * p), nrow = n_cal, ncol = p)
  Xc[, p]     <- rbinom(n_cal, 1, 0.4)
  Xc[, p - 1] <- rbinom(n_cal, 1, 0.6)
  mu_c <- .f_mu_count(Xc); het <- .f_tau_het(Xc)
  pz0  <- sigmoid(.f_mu_zi(Xc))
  pz1  <- sigmoid(.f_mu_zi(Xc) + .f_tau_zi(Xc))
  diff_fn <- function(tc)
    mean((1 - pz1) * exp(mu_c + tc + het) - (1 - pz0) * exp(mu_c)) - ate_target
  uniroot(diff_fn, lower = -5, upper = 5,
          extendInt = "yes", tol = 1e-4)$root
}

draw_dgp <- function(n, p, seed, ate_target) {
  tc <- .calibrate_tau_const_zi(p, ate_target)
  .simulate_zi(n, p, seed, tc)
}
'

METRICS_BLOCK <- '
summarize_fit <- function(res, d) {
  cate_post <- res$cate_post
  cate_mean <- colMeans(cate_post)
  cate_q025 <- apply(cate_post, 2, quantile, 0.025, names = FALSE)
  cate_q975 <- apply(cate_post, 2, quantile, 0.975, names = FALSE)
  cate_q250 <- apply(cate_post, 2, quantile, 0.250, names = FALSE)
  cate_q750 <- apply(cate_post, 2, quantile, 0.750, names = FALSE)
  true_mean_obs <- ifelse(d$z == 1, d$mu1, d$mu0)

  ate_mean <- mean(res$ate)
  ate_q025 <- quantile(res$ate, 0.025, names = FALSE)
  ate_q975 <- quantile(res$ate, 0.975, names = FALSE)
  ate_cov  <- as.integer(ate_q025 <= d$ate && d$ate <= ate_q975)

  rmse_yhat <- sqrt(mean((res$yhat_mean - true_mean_obs)^2))
  bias_yhat <- mean(res$yhat_mean - true_mean_obs)

  cate_bias  <- mean(cate_mean - d$cate)
  pehe       <- sqrt(mean((cate_mean - d$cate)^2))
  mae_cate   <- mean(abs(cate_mean - d$cate))
  cate_cor   <- suppressWarnings(stats::cor(cate_mean, d$cate))
  cate_cov95 <- mean(d$cate >= cate_q025 & d$cate <= cate_q975)
  cate_cov50 <- mean(d$cate >= cate_q250 & d$cate <= cate_q750)
  cate_ci_width95 <- mean(cate_q975 - cate_q025)

  data.frame(
    rmse_yhat = rmse_yhat, bias_yhat = bias_yhat,
    ate_mean  = ate_mean,  ate_q025  = ate_q025,
    ate_q975  = ate_q975,  ate_coverage = ate_cov,
    pehe      = pehe,      mae_cate  = mae_cate,
    cate_bias = cate_bias, cate_cor  = cate_cor,
    cate_cov95 = cate_cov95, cate_cov50 = cate_cov50,
    cate_ci_width95 = cate_ci_width95,
    stringsAsFactors = FALSE)
}

empty_metrics <- function() {
  data.frame(
    rmse_yhat = NA_real_, bias_yhat = NA_real_,
    ate_mean  = NA_real_, ate_q025  = NA_real_,
    ate_q975  = NA_real_, ate_coverage = NA_integer_,
    pehe      = NA_real_, mae_cate  = NA_real_,
    cate_bias = NA_real_, cate_cor  = NA_real_,
    cate_cov95 = NA_real_, cate_cov50 = NA_real_,
    cate_ci_width95 = NA_real_,
    stringsAsFactors = FALSE)
}
'

DRIVER_COUNT <- '
run_one_fit <- function(sim_id) {
  d <- draw_dgp(N, P, seed = sim_id, ate_target = ATE_TARGET)
  t0  <- Sys.time()
  res <- tryCatch(fit_method(d),
                  error = function(e) list(error = conditionMessage(e)))
  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  base <- data.frame(
    factor_varied = FACTOR_VAR,
    method = METHOD, dgp = DGP_NAME,
    N = N, P = P, ate_target = ATE_TARGET, sim = sim_id,
    true_ate = d$ate, pct_zero = d$pct_zero,
    cate_sd_true = sd(d$cate), elapsed_sec = elapsed,
    stringsAsFactors = FALSE)
  if (!is.null(res$error))
    return(cbind(base, empty_metrics(),
                 error = res$error, stringsAsFactors = FALSE))
  cbind(base, summarize_fit(res, d),
        error = NA_character_, stringsAsFactors = FALSE)
}

log_msg("BEGIN study=count method=%s dgp=%s N=%d P=%d ATE=%.2f n_sim=%d",
        METHOD, DGP_NAME, N, P, ATE_TARGET, N_SIM)
header_written <- FALSE
for (sim_id in seq_len(N_SIM)) {
  log_msg("  sim=%d/%d", sim_id, N_SIM)
  df <- run_one_fit(sim_id)
  write.table(df, OUT_CSV, sep = ",", row.names = FALSE,
              col.names = !header_written, append = header_written)
  header_written <- TRUE
  gc(verbose = FALSE)
}
log_msg("END study=count method=%s dgp=%s -> %s", METHOD, DGP_NAME, OUT_CSV)
'

DRIVER_ZI <- '
run_one_fit <- function(sim_id) {
  d <- draw_dgp(N, P, seed = sim_id, ate_target = ATE_TARGET)
  t0  <- Sys.time()
  res <- tryCatch(fit_method(d),
                  error = function(e) list(error = conditionMessage(e)))
  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  base <- data.frame(
    factor_varied = FACTOR_VAR,
    method = METHOD, dgp = DGP_NAME,
    N = N, P = P, ate_target = ATE_TARGET, sim = sim_id,
    true_ate = d$ate, pct_zero = d$pct_zero,
    pct_struct_zero = d$pct_struct_zero,
    cate_sd_true = sd(d$cate), elapsed_sec = elapsed,
    stringsAsFactors = FALSE)
  if (!is.null(res$error))
    return(cbind(base, empty_metrics(),
                 error = res$error, stringsAsFactors = FALSE))
  cbind(base, summarize_fit(res, d),
        error = NA_character_, stringsAsFactors = FALSE)
}

log_msg("BEGIN study=zi method=%s dgp=%s N=%d P=%d ATE=%.2f n_sim=%d",
        METHOD, DGP_NAME, N, P, ATE_TARGET, N_SIM)
header_written <- FALSE
for (sim_id in seq_len(N_SIM)) {
  log_msg("  sim=%d/%d", sim_id, N_SIM)
  df <- run_one_fit(sim_id)
  write.table(df, OUT_CSV, sep = ",", row.names = FALSE,
              col.names = !header_written, append = header_written)
  header_written <- TRUE
  gc(verbose = FALSE)
}
log_msg("END study=zi method=%s dgp=%s -> %s", METHOD, DGP_NAME, OUT_CSV)
'

## ----------------------------------------------------------------------------
##  Builder.
## ----------------------------------------------------------------------------

make_filename <- function(study, method, dgp, N, P, ATE) {
  sprintf("run_%s_%s_%s_N%d_P%d_ATE%d.R",
          study, method, dgp, N, P, as.integer(round(ATE * 100)))
}

make_csv_name <- function(study, method, dgp, N, P, ATE) {
  sprintf("sim_%s__%s__%s__N%d_P%d_ATE%d.csv",
          study, method, dgp, N, P, as.integer(round(ATE * 100)))
}

make_log_name <- function(study, method, dgp, N, P, ATE) {
  sprintf("sim_%s__%s__%s__N%d_P%d_ATE%d.log",
          study, method, dgp, N, P, as.integer(round(ATE * 100)))
}

write_one <- function(study, method, dgp_name, dgp_def, method_def, cell) {
  N <- cell$N; P <- cell$P; ATE <- cell$ATE; FV <- cell$factor_varied
  out_csv <- make_csv_name(study, method, dgp_name, N, P, ATE)
  out_log <- make_log_name(study, method, dgp_name, N, P, ATE)

  header <- sprintf(HEADER,
                    study, method, dgp_name,
                    N, P, ATE, FV,
                    out_csv, out_log,
                    method_def$pkg,
                    method, dgp_name, FV,
                    N, P, ATE,
                    dgp_def$kappa, dgp_def$count_model,
                    out_csv, out_log,
                    dgp_def$f_mu, dgp_def$f_tau)
  if (study == "zi") {
    header <- paste0(header,
                     sprintf(HEADER_ZI_EXTRA, dgp_def$f_mu_zi, dgp_def$f_tau_zi))
  }

  sim_block    <- if (study == "count") SIMULATE_COUNT else SIMULATE_ZI
  driver_block <- if (study == "count") DRIVER_COUNT   else DRIVER_ZI

  body <- paste0(header, sim_block, method_def$fit, METRICS_BLOCK, driver_block)

  out_dir <- if (study == "count") OUT_DIR_COUNT else OUT_DIR_ZI
  fname   <- file.path(out_dir, make_filename(study, method, dgp_name, N, P, ATE))
  writeLines(body, fname)
  fname
}

## ----------------------------------------------------------------------------
##  Generate all files.
## ----------------------------------------------------------------------------

written <- character(0)

for (i in seq_len(nrow(CELLS))) {
  cell <- as.list(CELLS[i, ])
  for (method in names(METHODS_COUNT)) {
    method_def <- METHODS_COUNT[[method]]
    for (dgp_name in names(DGPS_COUNT)) {
      f <- write_one("count", method, dgp_name,
                     DGPS_COUNT[[dgp_name]], method_def, cell)
      written <- c(written, f)
    }
  }
  for (method in names(METHODS_ZI)) {
    method_def <- METHODS_ZI[[method]]
    for (dgp_name in names(DGPS_ZI)) {
      f <- write_one("zi", method, dgp_name,
                     DGPS_ZI[[dgp_name]], method_def, cell)
      written <- c(written, f)
    }
  }
}

cat(sprintf("Wrote %d runner files.\n", length(written)))
cat(sprintf("  count: %d\n", sum(grepl("^runners_count/", written))))
cat(sprintf("  zi   : %d\n", sum(grepl("^runners_zi/",    written))))
