################################################################################
##  Simulation Study 1 -- Count-only DGPs (Poisson + Negative Binomial).
##
##  This file is both a *library* (DGPs, method wrappers, run_pair()) and an
##  *entry point* (main() runs every (method, DGP) pair sequentially).  The
##  per-pair runner files in runners_count/ source this file with the sentinel
##  SIM_NO_AUTO_RUN <- TRUE so that main() does *not* fire automatically; they
##  call run_pair() with their own (method, dgp).
##
##  Usage:
##    Rscript tests/sim_count.R                          # run everything serial
##    Rscript tests/runners_count/run_<method>_<dgp>.R   # run a single pair
##
##  Each (method, dgp) pair writes its own CSV under tests/results_count/, so
##  pairs can run safely in parallel (one R process per CPU).  See sim_count.md
##  for the full specification.
################################################################################

suppressMessages({
  library(countbart)
})

## ============================================================================
##  USER CONFIG
## ============================================================================

## Number of Monte Carlo replicates per (cell x dgp).  Seed of replicate s is
## exactly s (set.seed(sim_id)) so individual draws are reproducible from
## their (dgp, N, p, ATE, sim_id) tuple.
N_SIM   <- 100L

NBURN   <- 1000L
NSIM    <- 1000L
NTHIN   <- 1L

N_LEVELS   <- c(100L, 250L, 500L, 1000L, 2000L)
P_LEVELS   <- c(5L, 50L, 250L, 500L)
ATE_LEVELS <- c(0.5, 0.8, 1.25, 1.75, 2.5)

N_REF      <- 500L
P_REF      <- 5L
ATE_REF    <- 1.25

NB_KAPPA   <- 2

## 8 DGPs: 4 "ATE-focused" (mild heterogeneity, current spec) and 4
## "CATE-focused" (strong, sign-changing heterogeneity).  Both flavours share
## the same control-arm mean structure so the only difference is in tau_h(X).
DGP_NAMES_ATE  <- c("linear_poisson",      "nonlinear_poisson",
                    "linear_nb",           "nonlinear_nb")
DGP_NAMES_CATE <- c("linear_poisson_cate", "nonlinear_poisson_cate",
                    "linear_nb_cate",      "nonlinear_nb_cate")
DGP_NAMES      <- c(DGP_NAMES_ATE, DGP_NAMES_CATE)

METHOD_NAMES   <- c("countbcf", "bcf_gauss")

## ============================================================================
##  PATHS / LOGGING
## ============================================================================

.this_dir <- tryCatch(dirname(sys.frame(1)$ofile), error = function(e) "tests")
RESULTS_DIR <- file.path(.this_dir, "results_count")
if (!dir.exists(RESULTS_DIR)) dir.create(RESULTS_DIR, recursive = TRUE)

make_log_fun <- function(log_file) {
  function(fmt, ...) {
    msg <- sprintf(paste0("[%s] ", fmt),
                   format(Sys.time(), "%Y-%m-%d %H:%M:%S"), ...)
    cat(msg, "\n")
    cat(msg, "\n", file = log_file, append = TRUE)
  }
}

sigmoid <- function(x) 1 / (1 + exp(-x))

## ============================================================================
##  DATA-GENERATING PROCESSES
## ============================================================================

.true_propensity <- function(X) sigmoid(0.5 * X[, 1] - 0.4 * X[, 2])

.simulate_count <- function(n, p, seed,
                            f_mu_count, f_tau_het,
                            tau_const, kappa, dgp_name) {
  set.seed(seed)
  X <- matrix(rnorm(n * p), nrow = n, ncol = p)
  X[, p]     <- rbinom(n, 1, 0.4)
  X[, p - 1] <- rbinom(n, 1, 0.6)

  pi_x <- .true_propensity(X)
  Z    <- rbinom(n, 1, pi_x)

  log_lambda_0 <- f_mu_count(X)
  log_lambda_1 <- log_lambda_0 + tau_const + f_tau_het(X)
  log_lambda   <- ifelse(Z == 1, log_lambda_1, log_lambda_0)
  lambda       <- exp(log_lambda)

  Y <- if (is.na(kappa)) rpois(n, lambda)
       else              rnbinom(n, size = kappa, mu = lambda)

  mu0  <- exp(log_lambda_0)
  mu1  <- exp(log_lambda_1)
  cate <- mu1 - mu0

  list(
    y = Y, z = Z, x = X, pihat = pi_x,
    log_lambda_0 = log_lambda_0, log_lambda_1 = log_lambda_1,
    mu0 = mu0, mu1 = mu1, cate = cate, ate = mean(cate),
    tau_count = tau_const + f_tau_het(X), tau_const = tau_const,
    kappa = kappa, pct_zero = mean(Y == 0),
    pi_true = pi_x, dgp_name = dgp_name
  )
}

## Solve tau_const so mean(E[Y(1)] - E[Y(0)]) ~ ate_target on a held-out X.
.calibrate_tau_const <- function(p, f_mu_count, f_tau_het,
                                 ate_target, n_cal = 5000L,
                                 seed_cal = 9999991L) {
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
  mu0_c <- exp(f_mu_count(Xc))
  het_c <- f_tau_het(Xc)
  diff_fn <- function(tc) mean(mu0_c * (exp(tc + het_c) - 1)) - ate_target
  uniroot(diff_fn, lower = -5, upper = 5,
          extendInt = "yes", tol = 1e-4)$root
}

## ---- Mean and treatment-effect functions ---------------------------------
##
## ATE-focused (mild heterogeneity): per-unit tau varies a little around the
## mean.  CATE-focused (strong heterogeneity): per-unit tau is sign-changing
## and depends on multiple covariates -- the same calibrated marginal ATE is
## reachable, but unit-level estimation is much harder.

.f_mu_linear    <- function(X) 1.0 + 0.5 * X[, 1] - 0.3 * X[, 2] + 0.2 * X[, 3]
.f_mu_nonlinear <- function(X)
  0.8 + 0.6 * sin(X[, 1]) + 0.4 * X[, 2] * X[, 3] + 0.5 * (X[, 4] > 0)

## ATE-focused tau_h: small, single-covariate heterogeneity
.f_tau_linear      <- function(X) 0.20 * X[, 1]
.f_tau_nonlinear   <- function(X) 0.40 * (X[, 1] > 0) - 0.20 * X[, 5]

## CATE-focused tau_h: strong, multi-covariate, sign-changing heterogeneity
.f_tau_linear_cate    <- function(X) 0.60 * X[, 1] - 0.45 * X[, 2] + 0.30 * X[, 3]
.f_tau_nonlinear_cate <- function(X)
  0.80 * sin(X[, 1]) + 0.50 * (X[, 2] > 0) * X[, 3] - 0.40 * X[, 5]

## ---- DGP factory ---------------------------------------------------------
.dgp_meta <- list(
  linear_poisson         = list(f_mu = .f_mu_linear,    f_tau = .f_tau_linear,         kappa = NA),
  nonlinear_poisson      = list(f_mu = .f_mu_nonlinear, f_tau = .f_tau_nonlinear,      kappa = NA),
  linear_nb              = list(f_mu = .f_mu_linear,    f_tau = .f_tau_linear,         kappa = NB_KAPPA),
  nonlinear_nb           = list(f_mu = .f_mu_nonlinear, f_tau = .f_tau_nonlinear,      kappa = NB_KAPPA),
  linear_poisson_cate    = list(f_mu = .f_mu_linear,    f_tau = .f_tau_linear_cate,    kappa = NA),
  nonlinear_poisson_cate = list(f_mu = .f_mu_nonlinear, f_tau = .f_tau_nonlinear_cate, kappa = NA),
  linear_nb_cate         = list(f_mu = .f_mu_linear,    f_tau = .f_tau_linear_cate,    kappa = NB_KAPPA),
  nonlinear_nb_cate      = list(f_mu = .f_mu_nonlinear, f_tau = .f_tau_nonlinear_cate, kappa = NB_KAPPA)
)

draw_dgp <- function(dgp_name, n, p, seed, ate_target) {
  meta <- .dgp_meta[[dgp_name]]
  if (is.null(meta)) stop("Unknown DGP: ", dgp_name)
  tc <- .calibrate_tau_const(p, meta$f_mu, meta$f_tau, ate_target)
  .simulate_count(n, p, seed, meta$f_mu, meta$f_tau,
                  tc, kappa = meta$kappa, dgp_name = dgp_name)
}

dgp_to_count_model <- function(dgp_name) {
  if (grepl("poisson", dgp_name)) "poisson" else "nb"
}

## ============================================================================
##  METHOD WRAPPERS
##
##  Each returns a list with: ate (nsim), cate_post (nsim x n), yhat_mean (n).
##  The runner converts these to summary metrics.
## ============================================================================

fit_countbcf <- function(d, count_model) {
  fit <- countbcf(
    y = d$y, z = d$z,
    x_control = d$x, x_moderate = d$x,
    pihat = d$pihat,
    nburn = NBURN, nsim = NSIM, nthin = NTHIN,
    count_model = count_model,
    include_pihat = "control",
    update_interval = max(1L, (NBURN + NSIM) %/% 5),
    return_trees = FALSE
  )
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

fit_bcf_gauss <- function(d, count_model = NULL) {
  fit <- bcf_binary(
    y = d$y, z = d$z,
    x_control = d$x, x_moderate = d$x,
    pihat = d$pihat,
    nburn = NBURN, nsim = NSIM, nthin = NTHIN,
    include_pi = "control",
    sighat = sd(d$y),
    update_interval = max(1L, (NBURN + NSIM) %/% 5)
  )
  tau_post <- get_forest_fit(fit$moderate_fit, d$x)
  list(ate = rowMeans(tau_post),
       cate_post = tau_post,
       yhat_mean = colMeans(fit$yhat))
}

METHOD_FUNS <- list(
  countbcf  = fit_countbcf,
  bcf_gauss = fit_bcf_gauss
)

## ============================================================================
##  METRIC HELPERS
## ============================================================================

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

  cate_bias <- mean(cate_mean - d$cate)
  pehe      <- sqrt(mean((cate_mean - d$cate)^2))
  mae_cate  <- mean(abs(cate_mean - d$cate))
  cate_cor  <- suppressWarnings(stats::cor(cate_mean, d$cate))
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
    stringsAsFactors = FALSE
  )
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
    stringsAsFactors = FALSE
  )
}

## ============================================================================
##  ONE-FIT RUNNER
## ============================================================================

run_one_fit <- function(method, dgp_name, N, P, ate_target, sim_id) {
  d <- draw_dgp(dgp_name, n = N, p = P, seed = sim_id, ate_target = ate_target)
  count_model <- dgp_to_count_model(dgp_name)
  fit_fn <- METHOD_FUNS[[method]]
  if (is.null(fit_fn)) stop("Unknown method: ", method)

  t0  <- Sys.time()
  res <- tryCatch(fit_fn(d, count_model),
                  error = function(e) list(error = conditionMessage(e)))
  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  base <- data.frame(
    method = method, dgp = dgp_name,
    N = N, P = P, ate_target = ate_target,
    sim = sim_id, true_ate = d$ate, pct_zero = d$pct_zero,
    cate_sd_true = sd(d$cate),
    elapsed_sec = elapsed,
    stringsAsFactors = FALSE
  )
  if (!is.null(res$error)) {
    return(cbind(base, empty_metrics(),
                 error = res$error, stringsAsFactors = FALSE))
  }
  cbind(base, summarize_fit(res, d),
        error = NA_character_, stringsAsFactors = FALSE)
}

## ============================================================================
##  GRID + PER-PAIR DRIVER
## ============================================================================

build_cells <- function() {
  cells <- list()
  for (n in N_LEVELS)
    cells[[length(cells) + 1]] <- data.frame(
      factor_varied = "N", N = n, P = P_REF, ate_target = ATE_REF,
      stringsAsFactors = FALSE)
  for (p in P_LEVELS)
    cells[[length(cells) + 1]] <- data.frame(
      factor_varied = "P", N = N_REF, P = p, ate_target = ATE_REF,
      stringsAsFactors = FALSE)
  for (a in ATE_LEVELS)
    cells[[length(cells) + 1]] <- data.frame(
      factor_varied = "ATE", N = N_REF, P = P_REF, ate_target = a,
      stringsAsFactors = FALSE)
  cells_df <- do.call(rbind, cells)
  key <- paste(cells_df$N, cells_df$P, cells_df$ate_target, sep = "|")
  cells_df[!duplicated(key), , drop = FALSE]
}

## Run all cells x n_sim for ONE (method, dgp) pair.  Writes results to a
## pair-specific CSV under results_dir, plus a per-pair progress log.
run_pair <- function(method, dgp,
                     n_sim = N_SIM,
                     results_dir = RESULTS_DIR,
                     tag = "count") {
  if (!(method %in% METHOD_NAMES))
    stop("Unknown method: ", method)
  if (!(dgp %in% DGP_NAMES))
    stop("Unknown dgp: ", dgp)
  cells_df <- build_cells()

  base_name <- sprintf("sim_%s__%s__%s", tag, method, dgp)
  out_path  <- file.path(results_dir, paste0(base_name, ".csv"))
  log_path  <- file.path(results_dir, paste0(base_name, ".log"))
  if (file.exists(out_path)) file.remove(out_path)
  if (file.exists(log_path)) file.remove(log_path)
  log_msg   <- make_log_fun(log_path)
  header_written <- FALSE

  log_msg("BEGIN pair method=%s dgp=%s n_sim=%d cells=%d total_fits=%d",
          method, dgp, n_sim, nrow(cells_df), nrow(cells_df) * n_sim)

  for (i in seq_len(nrow(cells_df))) {
    row <- cells_df[i, ]
    for (sim_id in seq_len(n_sim)) {
      log_msg("[cell %d/%d %s] N=%d p=%d ate=%.2f sim=%d/%d",
              i, nrow(cells_df), row$factor_varied,
              row$N, row$P, row$ate_target, sim_id, n_sim)
      df <- tryCatch(
        run_one_fit(method, dgp, row$N, row$P, row$ate_target, sim_id),
        error = function(e) {
          log_msg("  TOP-LEVEL ERROR: %s", conditionMessage(e))
          NULL
        })
      if (is.null(df)) next
      df$factor_varied <- row$factor_varied
      df <- df[, c("factor_varied", setdiff(names(df), "factor_varied"))]
      write.table(df, out_path, sep = ",", row.names = FALSE,
                  col.names = !header_written, append = header_written)
      header_written <- TRUE
      gc(verbose = FALSE)
    }
  }
  log_msg("END pair method=%s dgp=%s -> %s", method, dgp, out_path)
  invisible(out_path)
}

## ============================================================================
##  MAIN (sequential -- one process, all pairs)
## ============================================================================

main <- function() {
  log_msg <- make_log_fun(file.path(RESULTS_DIR, "sim_count_main.log"))
  log_msg("MAIN start. Methods=%d DGPs=%d N_sim=%d",
          length(METHOD_NAMES), length(DGP_NAMES), N_SIM)
  for (m in METHOD_NAMES) for (d in DGP_NAMES) {
    log_msg("running pair: method=%s dgp=%s", m, d)
    run_pair(m, d, n_sim = N_SIM, tag = "count")
  }
  log_msg("MAIN done.")
}

## Auto-run unless a per-pair runner has set the sentinel.
if (!isTRUE(SIM_NO_AUTO_RUN <- ifelse(exists("SIM_NO_AUTO_RUN"),
                                       SIM_NO_AUTO_RUN, FALSE)))
  main()
