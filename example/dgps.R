################################################################################
##  Data-generating processes for testing Zero-Inflated Bayesian Causal Forests
##  (countbcf) and competing ZI / count benchmark models.
##
##  Four scenarios:
##    1. dgp_linear_zi          -- linear effects, moderate zero-inflation (~30%)
##    2. dgp_linear_extreme_zi  -- linear effects, heavy zero-inflation (~70%)
##    3. dgp_nonlinear_zi       -- nonlinear effects, moderate zero-inflation
##    4. dgp_nonlinear_extreme_zi -- nonlinear effects, heavy zero-inflation
##
##  Each function returns a named list with:
##    y, z, x, pihat              -- observed data + propensity score estimate
##    log_lambda_0, log_lambda_1  -- true log mean of the count component
##                                   under control / treatment
##    p_zi_0, p_zi_1              -- true P(structural zero) under control /
##                                   treatment
##    mu0, mu1                    -- true E[Y(0)|X], E[Y(1)|X]
##                                       = (1 - p_zi) * exp(log_lambda)
##    cate                        -- true CATE on the response scale (mu1 - mu0)
##    ate                         -- true ATE = mean(cate)
##    tau_count, tau_zi           -- true treatment effects on the log-rate
##                                   and ZI-logit scales (for diagnostics)
##    pct_struct_zero             -- realized fraction of structural zeros
##    pct_zero                    -- realized fraction of total zeros
##    dgp_name                    -- string identifier
##
##  All four DGPs share the same propensity model so results are comparable.
################################################################################

## ---- helpers ---------------------------------------------------------------

sigmoid <- function(z) 1 / (1 + exp(-z))

## propensity: confounded selection on x1 and x2
.true_propensity <- function(X) sigmoid(0.5 * X[, 1] - 0.4 * X[, 2])

## simulate one data set given the four functions:
##   f_mu_count(X)   : log lambda for control units
##   f_tau_count(X)  : treatment shift on log lambda
##   f_mu_zi(X)      : zero-inflation logit for control units
##   f_tau_zi(X)     : treatment shift on the ZI logit
.simulate <- function(n, p, seed,
                      f_mu_count, f_tau_count,
                      f_mu_zi,    f_tau_zi,
                      dgp_name) {
  set.seed(seed)

  ## covariates: mix of continuous + binary so models that prefer one or the
  ## other are not unfairly disadvantaged
  X <- matrix(rnorm(n * p), nrow = n, ncol = p)
  X[, p]     <- rbinom(n, 1, 0.4)             # last column binary
  X[, p - 1] <- rbinom(n, 1, 0.6)             # second-last column binary

  pi_x <- .true_propensity(X)
  Z    <- rbinom(n, 1, pi_x)

  log_lambda_0 <- f_mu_count(X)
  log_lambda_1 <- log_lambda_0 + f_tau_count(X)

  zi_logit_0 <- f_mu_zi(X)
  zi_logit_1 <- zi_logit_0 + f_tau_zi(X)
  p_zi_0 <- sigmoid(zi_logit_0)
  p_zi_1 <- sigmoid(zi_logit_1)

  ## realized
  log_lambda <- ifelse(Z == 1, log_lambda_1, log_lambda_0)
  p_zi       <- ifelse(Z == 1, p_zi_1,       p_zi_0)

  is_struct_zero <- rbinom(n, 1, p_zi)
  Y_count        <- rpois(n, exp(log_lambda))
  Y              <- ifelse(is_struct_zero == 1, 0L, Y_count)

  ## true potential-outcome means and CATE on the response scale
  mu0  <- (1 - p_zi_0) * exp(log_lambda_0)
  mu1  <- (1 - p_zi_1) * exp(log_lambda_1)
  cate <- mu1 - mu0

  ## "estimated" propensity: we use the true one here so DGP comparisons
  ## focus on the outcome model.  Replace by a fitted logit/RF in real
  ## benchmarks if you want to evaluate sensitivity to pi_hat.
  pihat <- pi_x

  list(
    y       = Y,
    z       = Z,
    x       = X,
    pihat   = pihat,
    log_lambda_0 = log_lambda_0,
    log_lambda_1 = log_lambda_1,
    p_zi_0  = p_zi_0,
    p_zi_1  = p_zi_1,
    mu0     = mu0,
    mu1     = mu1,
    cate    = cate,
    ate     = mean(cate),
    tau_count = f_tau_count(X),
    tau_zi    = f_tau_zi(X),
    pct_struct_zero = mean(is_struct_zero),
    pct_zero        = mean(Y == 0),
    pi_true = pi_x,
    dgp_name = dgp_name
  )
}

## ---- DGP 1: linear, moderate ZI -------------------------------------------
##
##  log lambda          = 1.0 + 0.5 X1 - 0.3 X2 + 0.2 X3
##  tau on log lambda   = 0.30 + 0.20 X1                  (linear CATE)
##  zi logit (control)  = -1.0 + 0.5 X2                    (~ 27% baseline ZI)
##  tau on zi logit     = -0.30 + 0.10 X3                  (treatment lowers ZI)
##
dgp_linear_zi <- function(n = 1000, p = 5, seed = 1) {
  .simulate(
    n = n, p = p, seed = seed,
    f_mu_count  = function(X) 1.0 + 0.5 * X[, 1] - 0.3 * X[, 2] + 0.2 * X[, 3],
    f_tau_count = function(X) 0.30 + 0.20 * X[, 1],
    f_mu_zi     = function(X) -1.0 + 0.5 * X[, 2],
    f_tau_zi    = function(X) -0.30 + 0.10 * X[, 3],
    dgp_name    = "linear_zi"
  )
}

## ---- DGP 2: linear, extreme ZI --------------------------------------------
##
##  Same count component, but the ZI logit is shifted up so that
##  P(structural zero) ~ 70% on average.  The treatment effect on the
##  ZI logit is also stronger.
##
dgp_linear_extreme_zi <- function(n = 1000, p = 5, seed = 2) {
  .simulate(
    n = n, p = p, seed = seed,
    f_mu_count  = function(X) 1.2 + 0.5 * X[, 1] - 0.3 * X[, 2] + 0.2 * X[, 3],
    f_tau_count = function(X) 0.30 + 0.20 * X[, 1],
    f_mu_zi     = function(X)  1.2 + 0.5 * X[, 2],          # ~77% ZI baseline
    f_tau_zi    = function(X) -0.60 + 0.20 * X[, 3],
    dgp_name    = "linear_extreme_zi"
  )
}

## ---- DGP 3: nonlinear, moderate ZI ----------------------------------------
##
##  log lambda            includes sin(X1), an interaction X2*X3 and a
##                        threshold indicator on X4
##  tau on log lambda     piecewise: stronger effect for X1 > 0
##  zi logit (control)    nonlinear via |X2| and exp-bumped X3
##  tau on zi logit       nonlinear interaction X1*X4
##
dgp_nonlinear_zi <- function(n = 1000, p = 5, seed = 3) {
  .simulate(
    n = n, p = p, seed = seed,
    f_mu_count  = function(X) {
      0.8 + 0.6 * sin(X[, 1]) + 0.4 * X[, 2] * X[, 3] + 0.5 * (X[, 4] > 0)
    },
    f_tau_count = function(X) {
      0.20 + 0.40 * (X[, 1] > 0) - 0.20 * X[, 5]
    },
    f_mu_zi = function(X) {
      -1.0 + 0.7 * abs(X[, 2]) - 0.5 * exp(-X[, 3]^2)
    },
    f_tau_zi = function(X) {
      -0.30 * X[, 1] * X[, 4]
    },
    dgp_name = "nonlinear_zi"
  )
}

## ---- DGP 4: nonlinear, extreme ZI -----------------------------------------
##
##  Same nonlinear forms as DGP 3, but ZI logit shifted up + stronger
##  treatment effect on the ZI part to push the structural-zero rate ~ 70-80%.
##
dgp_nonlinear_extreme_zi <- function(n = 1000, p = 5, seed = 4) {
  .simulate(
    n = n, p = p, seed = seed,
    f_mu_count  = function(X) {
      1.0 + 0.6 * sin(X[, 1]) + 0.4 * X[, 2] * X[, 3] + 0.5 * (X[, 4] > 0)
    },
    f_tau_count = function(X) {
      0.20 + 0.40 * (X[, 1] > 0) - 0.20 * X[, 5]
    },
    f_mu_zi = function(X) {
      1.4 + 0.7 * abs(X[, 2]) - 0.5 * exp(-X[, 3]^2)
    },
    f_tau_zi = function(X) {
      -0.60 - 0.30 * X[, 1] * X[, 4]
    },
    dgp_name = "nonlinear_extreme_zi"
  )
}

## ---- convenience: list all four -------------------------------------------

all_dgps <- function() {
  list(
    linear_zi             = dgp_linear_zi,
    linear_extreme_zi     = dgp_linear_extreme_zi,
    nonlinear_zi          = dgp_nonlinear_zi,
    nonlinear_extreme_zi  = dgp_nonlinear_extreme_zi
  )
}

## quick sanity printer ------------------------------------------------------

summarize_dgp <- function(d) {
  cat("DGP:                ", d$dgp_name,                    "\n", sep = "")
  cat("  n                : ", length(d$y),                   "\n", sep = "")
  cat("  p (covariates)   : ", ncol(d$x),                     "\n", sep = "")
  cat("  treatment frac   : ", round(mean(d$z), 3),           "\n", sep = "")
  cat("  total zeros      : ", round(d$pct_zero, 3),          "\n", sep = "")
  cat("  structural zeros : ", round(d$pct_struct_zero, 3),   "\n", sep = "")
  cat("  E[Y]             : ", round(mean(d$y), 3),           "\n", sep = "")
  cat("  max(Y)           : ", max(d$y),                       "\n", sep = "")
  cat("  true ATE         : ", round(d$ate, 3),               "\n", sep = "")
  cat("  CATE range       : ", round(min(d$cate), 3), " to ",
                                round(max(d$cate), 3),         "\n", sep = "")
  invisible(d)
}
