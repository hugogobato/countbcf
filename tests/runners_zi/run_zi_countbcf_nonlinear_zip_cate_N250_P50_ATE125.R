################################################################################
##  Self-contained Colab runner.
##  Study   : zi
##  Method  : countbcf
##  DGP     : nonlinear_zip_cate
##  Cell    : N=250  P=50  ATE=1.25  (factor varied: P)
##  N_SIM   : 100 (seeds 1..100 are independent draws)
##
##  Only external dependency: library(countbart).
##  Outputs sim_zi__countbcf__nonlinear_zip_cate__N250_P50_ATE125.csv and sim_zi__countbcf__nonlinear_zip_cate__N250_P50_ATE125.log to the working directory.
################################################################################

install.packages("remotes")
if (!require("devtools")) {
  install.packages("devtools")
}
devtools::install_github("hugogobato/countbcf")
library(countbcf)

N_SIM   <- 100L
NBURN   <- 1000L
NSIM    <- 1000L
NTHIN   <- 1L

METHOD       <- "countbcf"
DGP_NAME     <- "nonlinear_zip_cate"
FACTOR_VAR   <- "P"
N            <- 250L
P            <- 50L
ATE_TARGET   <- 1.25
KAPPA        <- NA_real_
COUNT_MODEL  <- "zipoisson"

OUT_CSV  <- "sim_zi__countbcf__nonlinear_zip_cate__N250_P50_ATE125.csv"
LOG_FILE <- "sim_zi__countbcf__nonlinear_zip_cate__N250_P50_ATE125.log"
if (file.exists(OUT_CSV))  file.remove(OUT_CSV)
if (file.exists(LOG_FILE)) file.remove(LOG_FILE)

log_msg <- function(fmt, ...) {
  msg <- sprintf(paste0("[%s] ", fmt),
                 format(Sys.time(), "%Y-%m-%d %H:%M:%S"), ...)
  cat(msg, "\n")
  cat(msg, "\n", file = LOG_FILE, append = TRUE)
}

sigmoid <- function(x) 1 / (1 + exp(-x))
.true_propensity <- function(X) sigmoid(0.5 * X[, 1] - 0.4 * X[, 2])

.f_mu_count <- function(X) 0.8 + 0.6 * sin(X[, 1]) + 0.4 * X[, 2] * X[, 3] + 0.5 * (X[, 4] > 0)
.f_tau_het  <- function(X) 0.80 * sin(X[, 1]) + 0.50 * (X[, 2] > 0) * X[, 3] - 0.40 * X[, 5]
.f_mu_zi    <- function(X) -1.0 + 0.7 * abs(X[, 2]) - 0.5 * exp(-X[, 3]^2)
.f_tau_zi   <- function(X) -0.60 * X[, 1] * X[, 4] + 0.40 * X[, 2] * (X[, 5] > 0)

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

