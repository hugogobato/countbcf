################################################################################
##  One run per DGP: fit bcf_binary (Gaussian BCF via multibart on counts) and
##  save the ATE posterior + per-unit summaries as CSV files for later plotting.
##
##  Outputs (under tests/results/):
##    bcf_ate_posterior.csv      long format: dgp, iter, ate
##    bcf_unit_estimates.csv     long format: dgp, unit, z, true_cate,
##                                            tau_post_mean, yhat_post_mean
##    bcf_summary.csv            per-DGP: rmse, bias, ate_mean, ate_q025, ate_q975,
##                                        true_ate, elapsed_sec
##
##  Usage (from package root):
##    Rscript tests/run_bcf_only.R
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
NBURN <- 1000
NSIM  <- 1000
NTHIN <- 1
N     <- 1000
P     <- 5

## ---- helpers --------------------------------------------------------------

## Recover the posterior of tau(x) from a bcf_binary fit.
## moderate_fit$tree_samples has basis_dim==1 (Z = z indicator), so get_forest_fit
## with Z=NULL returns the per-unit per-iteration tau coefficient. The wrapper
## multiplies by moderate_fit$scale = sdy so the output is on the response scale.
.tau_posterior <- function(fit, x_moderate, include_pi = "control", pihat = NULL) {
  X <- if (include_pi %in% c("moderate", "both")) {
    cbind(x_moderate, pihat)
  } else x_moderate
  get_forest_fit(fit$moderate_fit, X)            # nsim x n
}

## Recover mu(x) on the response scale (mu = mean function under z=0).
.mu_posterior <- function(fit, x_control, include_pi = "control", pihat = NULL) {
  X <- if (include_pi %in% c("control", "both")) {
    cbind(x_control, pihat)
  } else x_control
  get_forest_fit(fit$control_fit, X)             # nsim x n
}

## ---- per-DGP fit ----------------------------------------------------------
fit_one <- function(dgp_fn, dgp_seed) {

  d <- dgp_fn(n = N, p = P, seed = dgp_seed)
  cat("\n========================================================\n")
  summarize_dgp(d)

  cat("\n--- Fitting bcf_binary (Gaussian BCF on counts) ---\n")
  t0 <- Sys.time()
  fit <- bcf_binary(
    y          = d$y,
    z          = d$z,
    x_control  = d$x,
    x_moderate = d$x,
    pihat      = d$pihat,
    nburn      = NBURN,
    nsim       = NSIM,
    nthin      = NTHIN,
    include_pi = "control",
    update_interval = max(1L, (NBURN + NSIM) %/% 5)
  )
  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  cat(sprintf("  done in %.1f sec\n", elapsed))

  ## tau posterior (nsim x n, response scale)
  tau_post  <- .tau_posterior(fit, d$x, "control", d$pihat)
  mu_post   <- .mu_posterior(fit, d$x, "control", d$pihat)

  ## per-iter ATE = average tau across units
  ate_per_iter <- rowMeans(tau_post)

  ## per-unit posterior means
  tau_mean  <- colMeans(tau_post)
  yhat_mean <- colMeans(fit$yhat)
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
    tau_mean      = tau_mean,
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
}

## ---- write CSVs -----------------------------------------------------------
ate_long <- do.call(rbind, lapply(results, function(r) data.frame(
  dgp  = r$dgp_name,
  iter = seq_along(r$ate_per_iter),
  ate  = r$ate_per_iter,
  stringsAsFactors = FALSE
)))
write.csv(ate_long,
          file.path(RESULTS_DIR, "bcf_ate_posterior.csv"),
          row.names = FALSE)

unit_long <- do.call(rbind, lapply(results, function(r) data.frame(
  dgp            = r$dgp_name,
  unit           = seq_along(r$z),
  z              = r$z,
  true_cate      = r$true_cate,
  tau_post_mean  = r$tau_mean,
  yhat_post_mean = r$yhat_mean,
  stringsAsFactors = FALSE
)))
write.csv(unit_long,
          file.path(RESULTS_DIR, "bcf_unit_estimates.csv"),
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
          file.path(RESULTS_DIR, "bcf_summary.csv"),
          row.names = FALSE)

cat("\n=================== BCF SUMMARY ===================\n")
print(summary_df, row.names = FALSE, digits = 4)
cat("\nWrote:\n",
    "  ", file.path(RESULTS_DIR, "bcf_ate_posterior.csv"),  "\n",
    "  ", file.path(RESULTS_DIR, "bcf_unit_estimates.csv"), "\n",
    "  ", file.path(RESULTS_DIR, "bcf_summary.csv"),        "\n",
    sep = "")

invisible(results)
