#!/usr/bin/env Rscript
# =====================================================================
# 04  Propensity score, overlap/positivity, and the repeated-person-year
#     clustering question.
#
#  (a) what model produced the propensity score, and how well do the arms
#      overlap?  Five specifications are compared (constant, simple
#      logistic, the main spline logistic, a saturated spline logistic and
#      a probability random forest), with the counts of cells near 0 and 1
#      and the weighted balance they achieve.
#  (b) does the answer move if the score model changes, or if low-overlap
#      cells are trimmed away?
#  (c) are the featured subgroups sitting in the low-overlap tail?
#  (d) the LSS person-year table repeats each survivor across attained-age
#      and calendar-period cells.  The independent-Poisson likelihood
#      therefore understates variance.  Two checks: a cluster bootstrap of
#      the parametric model (design effect) and a CountBCF fit with random
#      intercepts at the survivor-cohort level.
# =====================================================================

source(file.path(Sys.getenv("CBCF_ROOT", "."), "Applied_Study_v2/code/common.R"))
suppressMessages({library(countbcf); library(ranger)})

d <- read_lss()
ds_solid <- build_ds(d, "solid", "cola02w10")

PS_LIST <- c("const", "simple", "main", "rich", "forest")

# =====================================================================
# (a) propensity model comparison and overlap diagnostics
# =====================================================================
stamp("(a) propensity / overlap diagnostics")

auc <- function(p, z) {
  r <- rank(p); n1 <- sum(z == 1); n0 <- sum(z == 0)
  (sum(r[z == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}
wsmd <- function(x, z, w) {              # IPW-weighted standardised difference
  w1 <- w[z == 1]; w0 <- w[z == 0]; x1 <- x[z == 1]; x0 <- x[z == 0]
  m1 <- sum(w1 * x1) / sum(w1); m0 <- sum(w0 * x0) / sum(w0)
  v1 <- sum(w1 * (x1 - m1)^2) / sum(w1); v0 <- sum(w0 * (x0 - m0)^2) / sum(w0)
  (m1 - m0) / sqrt((v1 + v0) / 2)
}

ps_diag <- list(); ps_dens <- list()
for (lv in c("low", "medium", "high")) {
  sub <- contrast_ds(ds_solid, lv)
  for (spec in PS_LIST) {
    p <- fit_ps(sub, spec)
    p <- pmin(pmax(p, 1e-6), 1 - 1e-6)
    w <- ifelse(sub$z == 1, 1 / p, 1 / (1 - p))          # IPW (ATE) weights
    ps_diag[[paste(lv, spec)]] <- tibble(
      contrast = lv, ps_model = spec, n = nrow(sub), treated_share = mean(sub$z),
      auc = auc(p, sub$z),
      min = min(p), q01 = quantile(p, .01), q50 = median(p),
      q99 = quantile(p, .99), max = max(p),
      n_lt_001 = sum(p < 0.01), n_lt_005 = sum(p < 0.05),
      n_gt_095 = sum(p > 0.95), n_gt_099 = sum(p > 0.99),
      pct_outside_005_095 = mean(p < 0.05 | p > 0.95),
      smd_agex_raw = wsmd(sub$agex, sub$z, rep(1, nrow(sub))),
      smd_agex_ipw = wsmd(sub$agex, sub$z, w),
      smd_age_raw = wsmd(sub$age, sub$z, rep(1, nrow(sub))),
      smd_age_ipw = wsmd(sub$age, sub$z, w),
      smd_female_raw = wsmd(as.numeric(sub$sex == 2), sub$z, rep(1, nrow(sub))),
      smd_female_ipw = wsmd(as.numeric(sub$sex == 2), sub$z, w),
      smd_hiro_raw = wsmd(as.numeric(sub$city == 1), sub$z, rep(1, nrow(sub))),
      smd_hiro_ipw = wsmd(as.numeric(sub$city == 1), sub$z, w))
    if (spec == "main")
      ps_dens[[lv]] <- tibble(contrast = lv, z = sub$z, pihat = p, pyr = sub$pyr,
                              sex_f = sub$sex_f, city_f = sub$city_f,
                              agex = sub$agex, age = sub$age)
  }
}
ps_tab <- bind_rows(ps_diag)
write_csv(ps_tab, file.path(TAB, "T23_propensity_overlap.csv"))
print(as.data.frame(ps_tab %>% dplyr::select(contrast, ps_model, auc, min, q01, q50,
  q99, max, n_lt_005, n_gt_095, smd_hiro_raw, smd_hiro_ipw, smd_agex_raw, smd_agex_ipw)),
  digits = 3)
overlap_cells <- bind_rows(ps_dens)
write_csv(overlap_cells, file.path(TAB, "T24_overlap_cells.csv.gz"))

# (c) is any featured subgroup concentrated in the low-overlap tail?
sub_hi <- contrast_ds(ds_solid, "high")
p_hi <- fit_ps(sub_hi, "main")
sg_ov <- add_bands(sub_hi) %>% mutate(pihat = p_hi) %>% filter(z == 1) %>%
  group_by(sex_f, agex_b) %>%
  summarise(n = n(), cases = sum(outcome), pyr = sum(pyr),
            pihat_min = min(pihat), pihat_q05 = quantile(pihat, .05),
            pihat_med = median(pihat), pihat_q95 = quantile(pihat, .95),
            pct_lt_005 = mean(pihat < 0.05), pct_gt_095 = mean(pihat > 0.95),
            .groups = "drop")
write_csv(sg_ov, file.path(TAB, "T25_subgroup_overlap.csv"))
print(as.data.frame(sg_ov), digits = 3)

# =====================================================================
# (b) does the answer move with the score model / with trimming?
# =====================================================================
stamp("(b) CountBCF refits across propensity specifications and trimming")

score_fit <- function(r, spec) {
  cl <- add_bands(r$cells); w <- cl$pyr; tr <- which(cl$z == 1)
  err <- exp(r$logrr) - 1
  ear <- ear_mat(r$logmu, r$logrr)
  gi <- function(f) intersect(which(f), tr)
  m <- draw_wmean(err, w, gi(cl$sex_f == "Male"))
  f <- draw_wmean(err, w, gi(cl$sex_f == "Female"))
  y <- draw_wmean(err, w, gi(cl$agex < 20)); o <- draw_wmean(err, w, gi(cl$agex >= 40))
  as_tibble(c(list(spec = spec, n_cells = nrow(cl), cases = sum(cl$outcome)),
              qsum(draw_wmean(err, w, tr), "att_err"),
              qsum(draw_wmean(ear, w, tr), "att_ear"),
              qsum(draw_serr(r$logmu, r$logrr, w, tr), "att_serr"),
              qsum(f - m, "d_female_male"), qsum(y - o, "d_young_old")))
}

ps_grid <- bind_rows(
  tibble(spec = paste0("PS = ", PS_LIST), ps = PS_LIST, trim = NA_real_,
         re = FALSE, re_var = "cluster"),
  tibble(spec = c("trim pihat outside [0.05,0.95]", "trim pihat outside [0.10,0.90]"),
         ps = "main", trim = c(0.05, 0.10), re = FALSE, re_var = "cluster"),
  tibble(spec = c("random intercepts: city x sex x age-at-exposure",
                  "random intercepts: finest survivor cohort"),
         ps = "main", trim = NA_real_, re = TRUE,
         re_var = c("cohort", "cluster")))

res <- mclapply(seq_len(nrow(ps_grid)), function(k) {
  g <- ps_grid[k, ]
  r <- tryCatch(fit_cbcf(sub_hi, label = g$spec, ps_spec = g$ps,
                         trim = if (is.na(g$trim)) NULL else g$trim,
                         cluster_re = g$re, re_var = g$re_var,
                         nburn = NBURN_SENS, nsim = NSIM_SENS),
                error = function(e) { stamp("FAILED ", g$spec, ": ", conditionMessage(e)); NULL })
  if (is.null(r)) return(NULL)
  out <- score_fit(r, g$spec) %>% mutate(re_sd = mean(r$re_sd), .after = spec)
  rm(r); gc(verbose = FALSE); out
}, mc.cores = n_workers(5), mc.preschedule = FALSE)
ps_sens <- bind_rows(res)
write_csv(ps_sens, file.path(TAB, "T26_propensity_cluster_sensitivity.csv"))
print(as.data.frame(ps_sens %>% dplyr::select(spec, n_cells, re_sd, att_err, att_err_lo,
  att_err_hi, d_female_male, d_female_male_lo, d_female_male_hi,
  d_young_old, d_young_old_lo, d_young_old_hi)), digits = 3)

# interval-width comparison: independent cells vs random intercepts
widths <- ps_sens %>% filter(grepl("^PS = main|^random intercepts", spec)) %>%
  transmute(spec, att_err, width_err = att_err_hi - att_err_lo,
            att_ear, width_ear = att_ear_hi - att_ear_lo)
write_csv(widths, file.path(TAB, "T27_cluster_interval_widths.csv"))
print(as.data.frame(widths), digits = 4)

# =====================================================================
# (d) cluster bootstrap of the parametric model -> design effect
# =====================================================================
stamp("(d) cluster bootstrap design effect (parametric)")
ic <- d %>% filter(dcatn > 0, band %in% BANDS) %>%
  mutate(outcome = solid, sex_f = factor(sex, 1:2, c("Male", "Female")),
         city_f = factor(city, 1:2, c("Hiroshima", "Nagasaki")),
         band = droplevels(factor(as.character(band), levels = BANDS)))
FBOOT <- outcome ~ band + sex_f * (ns(log(age), 4) + ns(agex, 3)) + city_f +
                   offset(log(pyr))
m0 <- glm(FBOOT, poisson("log"), ic)
b0 <- coef(m0)[grep("^band", names(coef(m0)))]
se0 <- sqrt(diag(vcov(m0)))[names(b0)]

cl_list <- split(seq_len(nrow(ic)), ic$cluster)
NBOOT <- 300L
boot <- mclapply(seq_len(NBOOT), function(bi) {
  set.seed(SEED + bi)
  idx <- unlist(cl_list[sample(length(cl_list), replace = TRUE)], use.names = FALSE)
  bb <- tryCatch(coef(glm(FBOOT, poisson("log"), ic[idx, ])),
                 error = function(e) NULL)
  if (is.null(bb)) return(NULL)
  bb[names(b0)]
}, mc.cores = n_workers(8))
Bm <- do.call(rbind, boot[!vapply(boot, is.null, logical(1))])
design <- tibble(term = names(b0), estimate = b0, se_model = se0,
                 se_cluster_boot = apply(Bm, 2, sd, na.rm = TRUE),
                 design_effect = (se_cluster_boot / se_model)^2,
                 se_ratio = se_cluster_boot / se_model, n_boot = nrow(Bm))
write_csv(design, file.path(TAB, "T28_cluster_bootstrap_design_effect.csv"))
print(as.data.frame(design), digits = 4)

stamp("04 done -> ", TAB)
