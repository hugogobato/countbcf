################################################################################
##  One run per DGP: fit countbcf() with count_model="zipoisson" and save the
##  posterior of the response-scale ATE plus per-unit summaries to CSV.
##
##  Response-scale CATE is obtained by counterfactual prediction at z=0 and z=1
##  using the saved tree_samples (six forests: mu_f, tau_f, mu_f0, tau_f0,
##  mu_f1, tau_f1). See get_forest_fit().
##
##  Outputs (under tests/results/):
##    countbcf_ate_posterior.csv   long format: dgp, iter, ate
##    countbcf_unit_estimates.csv  long format: dgp, unit, z, true_cate,
##                                              cate_post_mean, yhat_post_mean
##    countbcf_summary.csv         per-DGP: rmse, bias, ate_mean, ate_q025,
##                                          ate_q975, true_ate, elapsed_sec
##
##  Usage (from package root):
##    Rscript tests/run_countbcf_only.R
################################################################################

suppressMessages({
  library(countbart)
})

.this_dir <- tryCatch(
  dirname(sys.frame(1)$ofile),
  error = function(e) "tests"
)
source(file.path(.this_dir, "dgps.R"))

RESULTS_DIR <- file.path(.this_dir, "results")
if (!dir.exists(RESULTS_DIR)) dir.create(RESULTS_DIR, recursive = TRUE)

## ---- MCMC settings --------------------------------------------------------
NBURN <- 500
NSIM  <- 500
NTHIN <- 1
N     <- 1000
P     <- 5

sigmoid <- function(x) 1 / (1 + exp(-x))

## ---- per-DGP fit ----------------------------------------------------------
fit_one <- function(dgp_fn, dgp_seed) {

  d <- dgp_fn(n = N, p = P, seed = dgp_seed)
  cat("\n========================================================\n")
  summarize_dgp(d)

  cat("\n--- Fitting countbcf (zipoisson) ---\n")
  t0 <- Sys.time()
  fit <- countbcf(
    y           = d$y,
    z           = d$z,
    x_control   = d$x,
    x_moderate  = d$x,
    x_zero      = d$x,
    x_pos       = d$x,
    pihat       = d$pihat,
    nburn       = NBURN,
    nsim        = NSIM,
    nthin       = NTHIN,
    count_model = "zipoisson",
    include_pihat = "control",
    update_interval = max(1L, (NBURN + NSIM) %/% 5),
    return_trees = FALSE
  )
  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  cat(sprintf("  done in %.1f sec\n", elapsed))

  ## ---- counterfactual prediction for every unit on the response scale -----
  ##
  ## For a ZIP model:
  ##   log_lambda(x, z) = mu_f(x, pihat) + z * tau_f(x)
  ##   zi_logit(x, z)   = (mu_f0(x,pihat) + z*tau_f0(x)) - (mu_f1(x,pihat) + z*tau_f1(x))
  ##   P(SZ | x, z)     = sigmoid(zi_logit)
  ##   E[Y | x, z]      = (1 - P(SZ)) * exp(log_lambda)
  ##
  ## With return_trees=FALSE we pull the per-iteration raw coefficients
  ## directly from the fit (already in-sample predictions, original order,
  ## dim nsim x n).
  mu_f_post   <- fit$mu_f_post
  tau_f_post  <- fit$tau_f_post
  mu_f0_post  <- fit$mu_f0_post
  tau_f0_post <- fit$tau_f0_post
  mu_f1_post  <- fit$mu_f1_post
  tau_f1_post <- fit$tau_f1_post

  ## potential outcomes (nsim x n)
  log_lambda_0 <- mu_f_post
  log_lambda_1 <- mu_f_post + tau_f_post
  zi_logit_0   <- mu_f0_post              - mu_f1_post
  zi_logit_1   <- (mu_f0_post + tau_f0_post) - (mu_f1_post + tau_f1_post)
  p_zi_0       <- sigmoid(zi_logit_0)
  p_zi_1       <- sigmoid(zi_logit_1)
  mu0_post     <- (1 - p_zi_0) * exp(log_lambda_0)
  mu1_post     <- (1 - p_zi_1) * exp(log_lambda_1)
  cate_post    <- mu1_post - mu0_post                          # nsim x n

  ## per-iter ATE = average over units
  ate_per_iter <- rowMeans(cate_post)

  ## per-unit posterior means
  cate_mean <- colMeans(cate_post)
  z_mat     <- matrix(d$z, nrow = NSIM, ncol = length(d$z), byrow = TRUE)
  yhat_post <- (1 - z_mat) * mu0_post + z_mat * mu1_post
  yhat_mean <- colMeans(yhat_post)

  true_mean_obs <- ifelse(d$z == 1, d$mu1, d$mu0)
  rmse_yhat <- sqrt(mean((yhat_mean - true_mean_obs)^2))
  bias_yhat <- mean(yhat_mean - true_mean_obs)
  ate_mean  <- mean(ate_per_iter)
  ate_q025  <- quantile(ate_per_iter, 0.025, names = FALSE)
  ate_q975  <- quantile(ate_per_iter, 0.975, names = FALSE)

  cat(sprintf("  RMSE yhat=%.3f  bias=%+.3f  ATE_post=%.3f (95%% [%.3f, %.3f])  true ATE=%.3f\n",
              rmse_yhat, bias_yhat, ate_mean, ate_q025, ate_q975, d$ate))

  list(
    dgp_name      = d$dgp_name,
    elapsed_sec   = elapsed,
    true_ate      = d$ate,
    ate_per_iter  = ate_per_iter,
    cate_mean     = cate_mean,
    yhat_mean     = yhat_mean,
    true_cate     = d$cate,
    z             = d$z,
    rmse_yhat     = rmse_yhat,
    bias_yhat     = bias_yhat
  )
}

## ---- run all four DGPs ----------------------------------------------------
DGP_SEEDS <- list(
  linear_zi            = 11,
  linear_extreme_zi    = 12,
  nonlinear_zi         = 13,
  nonlinear_extreme_zi = 14
)
DGP_FNS <- list(
  linear_zi            = dgp_linear_zi,
  linear_extreme_zi    = dgp_linear_extreme_zi,
  nonlinear_zi         = dgp_nonlinear_zi,
  nonlinear_extreme_zi = dgp_nonlinear_extreme_zi
)

results <- list()
for (nm in names(DGP_FNS)) {
  results[[nm]] <- fit_one(DGP_FNS[[nm]], DGP_SEEDS[[nm]])
  gc(verbose = FALSE)
}

## ---- write CSVs -----------------------------------------------------------
ate_long <- do.call(rbind, lapply(results, function(r) data.frame(
  dgp  = r$dgp_name,
  iter = seq_along(r$ate_per_iter),
  ate  = r$ate_per_iter,
  stringsAsFactors = FALSE
)))
write.csv(ate_long,
          file.path(RESULTS_DIR, "countbcf_ate_posterior.csv"),
          row.names = FALSE)

unit_long <- do.call(rbind, lapply(results, function(r) data.frame(
  dgp            = r$dgp_name,
  unit           = seq_along(r$z),
  z              = r$z,
  true_cate      = r$true_cate,
  cate_post_mean = r$cate_mean,
  yhat_post_mean = r$yhat_mean,
  stringsAsFactors = FALSE
)))
write.csv(unit_long,
          file.path(RESULTS_DIR, "countbcf_unit_estimates.csv"),
          row.names = FALSE)

summary_df <- do.call(rbind, lapply(results, function(r) data.frame(
  dgp         = r$dgp_name,
  rmse_yhat   = r$rmse_yhat,
  bias_yhat   = r$bias_yhat,
  ate_mean    = mean(r$ate_per_iter),
  ate_q025    = quantile(r$ate_per_iter, 0.025, names = FALSE),
  ate_q975    = quantile(r$ate_per_iter, 0.975, names = FALSE),
  true_ate    = r$true_ate,
  elapsed_sec = r$elapsed_sec,
  stringsAsFactors = FALSE
)))
write.csv(summary_df,
          file.path(RESULTS_DIR, "countbcf_summary.csv"),
          row.names = FALSE)

cat("\n=================== CountBCF SUMMARY ===================\n")
print(summary_df, row.names = FALSE, digits = 4)
cat("\nWrote:\n",
    "  ", file.path(RESULTS_DIR, "countbcf_ate_posterior.csv"),  "\n",
    "  ", file.path(RESULTS_DIR, "countbcf_unit_estimates.csv"), "\n",
    "  ", file.path(RESULTS_DIR, "countbcf_summary.csv"),        "\n",
    sep = "")

invisible(results)
