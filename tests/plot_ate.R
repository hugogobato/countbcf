################################################################################
##  Read the ATE posteriors saved by run_bcf_only.R and run_countbcf_only.R
##  and produce density plots of the posterior ATE for each DGP, with a
##  vertical line at the true ATE.
##
##  Outputs (under tests/results/):
##    ate_density.png             4-panel facet, one DGP per panel
##    ate_density_table.csv       merged per-DGP summary across both models
##
##  Usage (from package root):
##    Rscript tests/plot_ate.R
################################################################################

suppressMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

.this_dir <- tryCatch(
  dirname(sys.frame(1)$ofile),
  error = function(e) "tests"
)
RESULTS_DIR <- file.path(.this_dir, "results")

bcf_post      <- read.csv(file.path(RESULTS_DIR, "bcf_ate_posterior.csv"),
                          stringsAsFactors = FALSE)
countbcf_post <- read.csv(file.path(RESULTS_DIR, "countbcf_ate_posterior.csv"),
                          stringsAsFactors = FALSE)
bcf_sum       <- read.csv(file.path(RESULTS_DIR, "bcf_summary.csv"),
                          stringsAsFactors = FALSE)
countbcf_sum  <- read.csv(file.path(RESULTS_DIR, "countbcf_summary.csv"),
                          stringsAsFactors = FALSE)

bcf_post$model      <- "BCF (Gaussian)"
countbcf_post$model <- "CountBCF (ZIP)"
ate_long <- rbind(bcf_post, countbcf_post)

## DGP order for plotting
dgp_levels <- c("linear_zi", "linear_extreme_zi",
                "nonlinear_zi", "nonlinear_extreme_zi")
ate_long$dgp <- factor(ate_long$dgp, levels = dgp_levels)

true_ate_df <- bcf_sum[, c("dgp", "true_ate")]
true_ate_df$dgp <- factor(true_ate_df$dgp, levels = dgp_levels)

## ---- plot -----------------------------------------------------------------
p <- ggplot(ate_long, aes(x = ate, fill = model, colour = model)) +
  geom_density(alpha = 0.35, linewidth = 0.5) +
  geom_vline(data = true_ate_df,
             aes(xintercept = true_ate),
             linetype = "dashed", colour = "black", linewidth = 0.6) +
  geom_text(data = true_ate_df,
            aes(x = true_ate, y = 0, label = sprintf("true ATE = %.2f", true_ate)),
            inherit.aes = FALSE,
            angle = 90, hjust = -0.1, vjust = -0.4, size = 3.0) +
  facet_wrap(~ dgp, scales = "free", ncol = 2) +
  labs(
    title    = "Posterior of the ATE: BCF (Gaussian) vs CountBCF (ZIP)",
    subtitle = "1 simulated run per DGP; n = 1000, nburn = nsim = 1000",
    x        = "ATE (response scale)",
    y        = "Posterior density",
    fill     = NULL,
    colour   = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "top")

out_png <- file.path(RESULTS_DIR, "ate_density.png")
ggsave(out_png, p, width = 9, height = 6.5, dpi = 150)

## ---- merged summary table ------------------------------------------------
merged <- merge(
  bcf_sum[, c("dgp", "ate_mean", "ate_q025", "ate_q975", "rmse_yhat", "bias_yhat", "true_ate")],
  countbcf_sum[, c("dgp", "ate_mean", "ate_q025", "ate_q975", "rmse_yhat", "bias_yhat")],
  by = "dgp", suffixes = c("_bcf", "_countbcf")
)
merged <- merged[match(dgp_levels, merged$dgp), ]
write.csv(merged, file.path(RESULTS_DIR, "ate_density_table.csv"),
          row.names = FALSE)

cat("\n=================== MERGED SUMMARY ===================\n")
print(merged, row.names = FALSE, digits = 4)
cat("\nSaved plot to ", out_png, "\n", sep = "")
cat("Saved table to ", file.path(RESULTS_DIR, "ate_density_table.csv"), "\n",
    sep = "")
