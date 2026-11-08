# Regression check: evaluating extra rows must not change the fitted chain.
# Run after installing the current source. Uses small, fixed-seed fits.
library(countbcf)
set.seed(902)
x <- cbind(a = rnorm(50), b = rep(0:1, 25))
ps <- plogis(x[, 1] / 2)
z <- rbinom(50, 1, ps)
y <- rpois(50, exp(.3 + x[, 1] / 3 + z * .2))
ix <- c(1L, 7L, 26L, 50L)
for (model in c("poisson", "nb")) {
  for (inc in c("none", "control", "moderate", "both")) {
    args <- list(y = y, z = z, x_control = x, x_moderate = x[, 1, drop = FALSE],
      pihat = ps, include_pihat = inc, count_model = model,
      nburn = 100, nsim = 12, ntree_control = 4, ntree_moderate = 2,
      update_interval = 1000, return_trees = FALSE)
    set.seed(722)
    original <- do.call(countbcf, args)
    set.seed(722)
    predicted <- do.call(countbcf, c(args, list(x_control_predict = x[ix, , drop = FALSE],
      x_moderate_predict = x[ix, 1, drop = FALSE], pihat_predict = ps[ix])))
    stopifnot(identical(original$mu_f_post, predicted$mu_f_post),
      identical(original$tau_f_post, predicted$tau_f_post),
      identical(original$kappa, predicted$kappa),
      isTRUE(all.equal(unname(predicted$mu_f_predict),
        unname(original$mu_f_post[, ix, drop = FALSE]), tolerance = 1e-12)),
      isTRUE(all.equal(unname(predicted$tau_f_predict),
        unname(original$tau_f_post[, ix, drop = FALSE]), tolerance = 1e-12)))
  }
}
cat("Prediction-row regression checks passed for both count models and all propensity placements.\n")
