#!/usr/bin/env Rscript
# Reconstruct final city/profile tables from compact archived draws; no fitting.
# CBCF_REPLICATION_INPUT defaults to reproducibility/applied (models + tables).
# CBCF_OUTDIR chooses an output analysis directory, as in the other scripts.
source(file.path(Sys.getenv("CBCF_ROOT", "."), "Applied_Study/code/common.R"))
stopifnot(requireNamespace("posterior", quietly = TRUE), requireNamespace("coda", quietly = TRUE))
INPUT <- Sys.getenv("CBCF_REPLICATION_INPUT", file.path(ROOT, "reproducibility", "applied"))
artifact <- function(name, table = FALSE) {
  candidates <- c(file.path(INPUT, if (table) "tables" else "models", name),
                  file.path(INPUT, name), file.path(if (table) TAB else OUTDIR, name))
  found <- candidates[file.exists(candidates)]
  if (!length(found)) stop("Missing archived artifact: ", name)
  found[[1]]
}

# Profile draws preserve the exact grid, city weights and fitting identities.
profile_archive <- readRDS(artifact("profile_age30_70_draws.rds"))
chains <- profile_archive$chains
grid <- profile_archive$grid
stopifnot(length(chains) >= 2L, all(grid$agex == 30), all(grid$age == 70))
pf <- readRDS(artifact("preston_fit.rds"))$preston_par
pt <- read_csv(artifact("T08_preston_err_ear_refit.csv", table = TRUE), show_col_types = FALSE)
profile_draws <- lapply(chains, function(ch) {
  out <- list()
  for (s in 1:2) {
    ix <- which(grid$sex == s)
    nm <- if (s == 1) "Male" else "Female"
    w <- grid$city_weight[ix]
    r0 <- exp(ch$logmu[, ix, drop = FALSE])
    r1 <- exp(ch$logmu[, ix, drop = FALSE] + ch$logrr[, ix, drop = FALSE])
    out[[paste(nm, "City mixture", sep = ":")]] <- as.numeric(r1 %*% w / (r0 %*% w)) - 1
    for (i in ix) out[[paste(nm, as.character(grid$city_f[i]), sep = ":")]] <- exp(ch$logrr[, i]) - 1
  }
  out
})
res <- bind_rows(lapply(names(profile_draws[[1]]), function(nm) {
  ss <- strsplit(nm, ":", fixed = TRUE)[[1]]
  sx <- if (ss[1] == "Male") 1L else 2L
  dr <- lapply(profile_draws, `[[`, nm)
  x <- unlist(dr); dc <- unique(grid$dose_contrast_Gy[grid$sex == sx])
  pr <- pt %>% filter(model == "ERR", grepl(if (sx == 1) "^beta_male" else "^beta_female", parameter))
  tibble(sex = ss[1], city = ss[2], exposure_age = 30, attained_age = 70,
         dose_contrast_Gy = dc,
         band_err = mean(x), band_err_lo = unname(quantile(x, .025)),
         band_err_hi = unname(quantile(x, .975)),
         err_per_Gy_rescaled = mean(x) / dc,
         err_per_Gy_rescaled_lo = unname(quantile(x, .025)) / dc,
         err_per_Gy_rescaled_hi = unname(quantile(x, .975)) / dc,
         preston_err_per_Gy = pf$ERR[pf$p0 + sx],
         preston_se = pr$se, preston_lo90 = pr$lo, preston_hi90 = pr$hi,
         interval_note = "CountBCF pooled draw 2.5%/97.5% quantiles; inspect Rhat/ESS before interpreting as CrI; Preston 90% CI; CountBCF/Gy rescales a binary group contrast")
}))
write_csv(res, file.path(TAB, "T42_profile_age30_70.csv"))
diag <- bind_rows(lapply(names(profile_draws[[1]]), function(nm) {
  X <- do.call(cbind, lapply(profile_draws, `[[`, nm))
  tibble(quantity = nm, rhat = posterior::rhat(X),
         ess_bulk = posterior::ess_bulk(X), ess_tail = posterior::ess_tail(X),
         mcse_mean = posterior::mcse_mean(X),
         mcse_q025 = posterior::mcse_quantile(X, probs = .025, names = FALSE),
         mcse_q975 = posterior::mcse_quantile(X, probs = .975, names = FALSE))
}))
write_csv(diag, file.path(TAB, "T42b_profile_mcmc.csv"))
if (any(!is.finite(diag$rhat) | diag$rhat > 1.01 | diag$ess_bulk < 400))
  warning("Fixed-profile MCMC diagnostics need attention; retained-draw quantiles are not validated credible intervals. See T42b_profile_mcmc.csv.")
write_csv(bind_rows(lapply(chains, function(ch) mutate(ch$main_summary, chain = ch$chain, seed = ch$seed))),
          file.path(TAB, "T42c_profile_main_fit_check.csv"))
write_csv(grid, file.path(TAB, "T42d_profile_grid.csv"))

# City diagnostics retain both within-city summaries and their draw-wise contrasts.
city_archive <- readRDS(artifact("city_statistic_draws.rds"))
chains <- city_archive$chains
sub_low <- contrast_ds(build_ds(read_lss(), "solid", "cola02w10"), "low")
QOI <- names(chains[[1]])
mcmc_tab <- bind_rows(lapply(QOI, function(q) {
  ch <- lapply(chains, `[[`, q)
  re <- rhat_ess(ch)
  x  <- unlist(ch)
  # Quantile uncertainty is different from uncertainty in the posterior mean.
  # posterior expects iterations x chains and estimates the ESS for each tail.
  qmcse <- posterior::mcse_quantile(do.call(cbind, ch), probs = c(.025, .975))
  tibble(quantity = q, mean = mean(x), lo = quantile(x, .025),
         med = quantile(x, .5), hi = quantile(x, .975), sd = sd(x),
         rhat = re[["rhat"]], ess = re[["ess"]],
         mcse = sd(x) / sqrt(re[["ess"]]),
         mcse_mean = sd(x) / sqrt(re[["ess"]]),
         mcse_q025 = unname(qmcse[1]), mcse_q975 = unname(qmcse[2]),
         # Descriptive numerical margin for the lower endpoint, not a test.
         lo_over_mcse = quantile(x, .025) / unname(qmcse[1]))
}))
write_csv(mcmc_tab, file.path(TAB, "T34_city_footing_mcmc.csv"))
print(as.data.frame(mcmc_tab), digits = 4)

pooled <- lapply(QOI, function(q) unlist(lapply(chains, `[[`, q)))
names(pooled) <- QOI

# =====================================================================
# 2. E-values and descriptive per-Gy rescaling for within-city contrasts
# =====================================================================
dose_contrast <- function(s, idx_city = NULL) {
  k <- if (is.null(idx_city)) rep(TRUE, nrow(s)) else idx_city
  weighted.mean(s$dose[k & s$z == 1], s$pyr[k & s$z == 1]) -
    weighted.mean(s$dose[k & s$z == 0], s$pyr[k & s$z == 0])
}
dc_all <- dose_contrast(sub_low)
dc_hi  <- dose_contrast(sub_low, sub_low$city_f == "Hiroshima")
dc_na  <- dose_contrast(sub_low, sub_low$city_f == "Nagasaki")

ev_row <- function(nm, draws, dc, difference = FALSE) {
  e   <- as.numeric(c(mean(draws), quantile(draws, c(.025, .975))))
  # A difference of ERRs is not a rate ratio, so it has no E-value obtained
  # by adding one. Nor is dividing the city difference by pooled dose a slope.
  R   <- if (difference) rep(NA_real_, 3) else 1 + e
  evp <- if (difference) NA_real_ else evalue_pt(R[1])
  evc <- if (difference) NA_real_ else evalue_ci(R[2], R[3])
  scaled <- if (difference) rep(NA_real_, 3) else e / dc
  tibble(quantity = nm, dose_contrast_Gy = if (difference) NA_real_ else dc,
         err = e[1], err_lo = e[2], err_hi = e[3],
         rr = R[1], rr_lo = R[2], rr_hi = R[3],
         err_per_Gy = scaled[1], err_per_Gy_lo = scaled[2],
         err_per_Gy_hi = scaled[3],
         evalue_point = evp, evalue_ci = evc,
         interpretation = if (difference) "Difference of ERRs; RR, E-value and dose slope not defined" else
           "City-specific band ERR; ERR/Gy is a descriptive dose-contrast rescaling")
}
ev_city <- bind_rows(
  ev_row("Hiroshima, ratio of rates", pooled$serr_hi, dc_hi),
  ev_row("Nagasaki, ratio of rates",  pooled$serr_na, dc_na),
  ev_row("Nagasaki - Hiroshima, ratio of rates", pooled$d_naga_hiro_serr, dc_all, difference = TRUE),
  ev_row("Hiroshima, per-cell mean", pooled$err_hi, dc_hi),
  ev_row("Nagasaki, per-cell mean",  pooled$err_na, dc_na),
  ev_row("Nagasaki - Hiroshima, per-cell mean", pooled$d_naga_hiro_err, dc_all, difference = TRUE))
write_csv(ev_city, file.path(TAB, "T35_city_evalues.csv"))
print(as.data.frame(ev_city), digits = 4)


dir.create(file.path(OUTDIR, "logs"), recursive = TRUE, showWarnings = FALSE)
writeLines(capture.output(sessionInfo()), file.path(OUTDIR, "logs", "sessionInfo_11.txt"))
stamp("11 archived posterior tables rebuilt -> ", TAB)
