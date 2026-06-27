#!/usr/bin/env Rscript
# =====================================================================
# CountBCF applied study on the RERF Life Span Study (LSS) cancer
# incidence person-year tables (lssinc07.csv).
#
# Design (addresses the issues raised in CountBCF_Applied_Study_Review.md):
#   * Three-level treatment from the continuous organ dose:
#       control (low) = [0, 0.1) Gy
#       moderate      = [0.1, 0.5) Gy
#       high          = [0.5, 2.0] Gy          (>2 Gy dropped: dose reliability)
#     Analysis restricted to un4gy == 1 (standard analysis cohort).
#   * Parametric Poisson / NB GLMs are used for the ATE ONLY (rate ratio +
#     marginal excess rate), no propensity score -- exactly what an applied
#     radiation epidemiologist would run. This is the convergence reference.
#   * CountBCF (NB) carries the CATE: per-cell excess-rate heterogeneity by
#     age-at-exposure, sex and attained age, with full posterior intervals.
#     CountBCF gets a propensity score (it does the RIC adjustment for us).
#
# Estimands (all on the response/rate scale, comparable across methods):
#   * CATE_i  = excess cases per 100,000 person-years for cell i
#             = ( exp(mu_f + tau_f) - exp(mu_f) ) * 1e5     [EAR-like]
#   * log RR  = tau_f                                       [ERR-like]
#   * ATE     = person-year-weighted mean of CATE_i  = total excess cases /
#               total person-years * 1e5
#
# Runs end-to-end locally; also pasteable into Colab (set INSTALL_CBCF=TRUE).
# =====================================================================

suppressMessages({
  library(MASS)
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(tibble)
})

INSTALL_CBCF <- FALSE
if (INSTALL_CBCF && !requireNamespace("countbcf", quietly = TRUE)) {
  install.packages("remotes"); remotes::install_github("hugogobato/countbcf")
}
suppressMessages(library(countbcf))

set.seed(20260622)

# ---- paths ----------------------------------------------------------
find_base <- function() {
  cands <- c("Applied_Study/RERF_data",
             "RERF_data",
             file.path(dirname(sys.frame(1)$ofile %||% "."), "RERF_data"))
  for (p in cands) if (file.exists(file.path(p, "lssinc07.csv"))) return(normalizePath(p))
  stop("Could not locate lssinc07.csv; run from the repo root or Applied_Study/.")
}
`%||%` <- function(a, b) if (is.null(a)) b else a
BASE   <- tryCatch(find_base(), error = function(e) "Applied_Study/RERF_data")
OUTDIR <- file.path(dirname(BASE), "countbcf_applied_analysis")
dir.create(OUTDIR,                 showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(OUTDIR, "plots"), showWarnings = FALSE)
dir.create(file.path(OUTDIR, "tables"), showWarnings = FALSE)

solid <- read_csv(file.path(BASE, "lssinc07.csv"), show_col_types = FALSE)
names(solid) <- trimws(names(solid))

NBURN <- 1000L; NSIM <- 1000L          # MCMC budget (paper-consistent)
NSIM_PARAM <- 4000L                    # coef draws for parametric ATE intervals
DOSE_CAP <- 2.0
LOW_HI <- 0.1; MOD_HI <- 0.5           # dose-band cut points (Gy)

theme_set(theme_bw(base_size = 13))
sex_lab <- c(`1` = "Male", `2` = "Female")

# =====================================================================
# 1. Dataset builder + EDA
# =====================================================================
build_dataset <- function(data, endpoint, dose_var) {
  stopifnot(endpoint %in% names(data), dose_var %in% names(data))
  d <- data %>%
    mutate(dose = ifelse(.data[[dose_var]] < 0, NA_real_, .data[[dose_var]])) %>%
    filter(un4gy == 1, !is.na(dose), dose >= 0, dose <= DOSE_CAP,
           !is.na(.data[[endpoint]]), !is.na(pyr), pyr > 0) %>%
    mutate(
      outcome    = .data[[endpoint]],
      dose_level = factor(case_when(
                     dose <  LOW_HI ~ "control",
                     dose <  MOD_HI ~ "moderate",
                     TRUE           ~ "high"),
                   levels = c("control", "moderate", "high")),
      sex_f = factor(sex), city_f = factor(city),
      agx_f = factor(agxcat), age_f = factor(agecat)) %>%
    select(outcome, dose, pyr, dose_level, sex, city, sex_f, city_f,
           agx_f, age_f, agxcat, agecat, time, gdist, agex, age, year)
  d
}

eda_summary <- function(d, endpoint, dose_var) {
  d %>% group_by(dose_level) %>%
    summarise(endpoint = endpoint, dose_var = dose_var,
      n_cells = n(), cases = sum(outcome), pyr = sum(pyr),
      rate_per_100k = cases / pyr * 1e5,
      mean_dose = mean(dose), median_dose = median(dose),
      pct_zero_cells = mean(outcome == 0),
      var_mean_ratio = var(outcome) / mean(outcome), .groups = "drop")
}

# =====================================================================
# 2. Parametric ATE (Poisson + NB), no propensity -- the convergence ref
# =====================================================================
# Adjusted mean structure: dose level + sex + city + age-at-exposure +
# attained-age (factors). 'time' omitted to keep the design full rank
# (attained age = age-at-exposure + time since exposure).
F_ADJ <- outcome ~ dose_level + sex_f + city_f + agx_f + age_f
F_MIN <- outcome ~ dose_level

# Marginal excess rate (per 1e5 PY) and rate ratio for a treated level vs
# control, by g-computation (standardisation over the covariate distribution).
# Because the mean structure is multiplicative, the marginal excess rate is
# base_rate * (RR - 1); we carry uncertainty through the (stable) treatment
# coefficient. Simulating the full coefficient vector instead is numerically
# unstable here -- sparse age x age-at-exposure factor levels produce extreme
# draws -- which is itself the parametric fragility the paper highlights.
marginal_ate <- function(model, d, treated, nsim = NSIM_PARAM) {
  b <- coef(model); V <- vcov(model)
  tt <- delete.response(terms(model))
  d0 <- d; d0$dose_level <- factor("control", levels = levels(d$dose_level))
  X0 <- model.matrix(tt, d0)
  keep <- intersect(colnames(X0), names(b))
  rate0 <- exp(as.numeric(X0[, keep, drop = FALSE] %*% b[keep]))   # per-PY control rate
  base_rate_100k <- weighted.mean(rate0, d$pyr) * 1e5             # standardised baseline
  lv <- paste0("dose_level", treated)
  bl <- b[lv]; se <- sqrt(V[lv, lv])
  rr_draws <- exp(rnorm(nsim, bl, se))
  excess_draws <- base_rate_100k * (rr_draws - 1)
  tibble(treated = treated, base_rate_100k = base_rate_100k,
         excess_rate_100k = base_rate_100k * (exp(bl) - 1),
         excess_lo = quantile(excess_draws, .025),
         excess_hi = quantile(excess_draws, .975),
         rr = exp(bl), rr_lo = exp(bl - 1.96 * se), rr_hi = exp(bl + 1.96 * se))
}

fit_parametric <- function(d, endpoint) {
  models <- list(
    `Poisson (minimal)`  = glm(F_MIN, poisson("log"), d, offset = log(pyr)),
    `Poisson (adjusted)` = glm(F_ADJ, poisson("log"), d, offset = log(pyr)),
    `NB (minimal)`       = glm.nb(update(F_MIN, . ~ . + offset(log(pyr))), d),
    `NB (adjusted)`      = glm.nb(update(F_ADJ, . ~ . + offset(log(pyr))), d))
  overdisp <- function(m) if (inherits(m, "negbin")) NA_real_ else
    sum(residuals(m, "pearson")^2) / m$df.residual
  out <- lapply(names(models), function(nm) {
    m <- models[[nm]]
    bind_rows(lapply(c("moderate", "high"), function(lv) marginal_ate(m, d, lv))) %>%
      mutate(endpoint = endpoint, model = nm, aic = AIC(m),
             poisson_overdisp = overdisp(m), .before = 1)
  })
  bind_rows(out)
}

# =====================================================================
# 3. CountBCF (NB) -- the CATE engine
# =====================================================================
fit_cbcf <- function(d, endpoint, treated, label) {
  sub <- d %>% filter(dose_level %in% c("control", treated)) %>% droplevels()
  z   <- as.integer(sub$dose_level == treated)

  # Propensity from the demographic confounders only. Distance-from-hypocentre
  # is deliberately EXCLUDED: it is a near-deterministic proxy for dose, so
  # conditioning on it collapses positivity/overlap and flips the estimate.
  # A-bomb dose is quasi-random given demographics, so this score is well
  # overlapped; CountBCF folds it into the prognostic forest (the BCF RIC fix).
  ps    <- glm(z ~ city_f + sex_f + agex + age, binomial, data = sub)
  pihat <- as.numeric(predict(ps, type = "response"))

  xc <- model.matrix(~ agex + age + sex + city + year + time, sub)[, -1]
  xm <- model.matrix(~ agex + age + sex, sub)[, -1]      # CATE drivers (shrunk hard)

  fit <- countbcf(y = sub$outcome, z = z, x_control = xc, x_moderate = xm,
                  pihat = pihat, offset = sub$pyr, count_model = "nb",
                  include_pihat = "control", nburn = NBURN, nsim = NSIM)

  rate0 <- exp(fit$mu_f_post)                    # nsim x n, per-PY rate at control
  rate1 <- exp(fit$mu_f_post + fit$tau_f_post)   # per-PY rate at treated
  cate  <- (rate1 - rate0) * 1e5                 # excess per 100k PY
  logRR <- fit$tau_f_post

  w <- sub$pyr / sum(sub$pyr)
  ate_draws <- as.numeric(cate %*% w)            # PY-weighted ATE per draw
  rr_draws  <- as.numeric((rate1 %*% sub$pyr) / (rate0 %*% sub$pyr))

  sub$cate      <- colMeans(cate)
  sub$cate_lo   <- apply(cate, 2, quantile, .025)
  sub$cate_hi   <- apply(cate, 2, quantile, .975)
  sub$logRR     <- colMeans(logRR)
  sub$rr        <- colMeans(exp(logRR))
  sub$pihat     <- pihat
  sub$contrast  <- label
  sub$endpoint  <- endpoint

  list(
    cells = sub,
    ate = tibble(endpoint = endpoint, contrast = label, method = "CountBCF (NB)",
                 excess_rate_100k = mean(ate_draws),
                 excess_lo = quantile(ate_draws, .025),
                 excess_hi = quantile(ate_draws, .975),
                 rr = mean(rr_draws), rr_lo = quantile(rr_draws, .025),
                 rr_hi = quantile(rr_draws, .975),
                 kappa = mean(fit$kappa)))
}

# =====================================================================
# 4. Run
# =====================================================================
ENDPOINTS <- tibble::tribble(
  ~endpoint, ~dose_var,    ~nice,
  "solid",   "cola02w10",  "All solid cancers (colon dose)",
  "stomach", "stoa02w10",  "Stomach cancer (stomach dose)")

eda_all <- list(); param_all <- list(); cbcf_ate <- list(); cells_all <- list()

for (i in seq_len(nrow(ENDPOINTS))) {
  ep <- ENDPOINTS$endpoint[i]; dv <- ENDPOINTS$dose_var[i]
  message(">>> ", ENDPOINTS$nice[i])
  d <- build_dataset(solid, ep, dv)

  eda_all[[ep]]   <- eda_summary(d, ep, dv)
  param_all[[ep]] <- fit_parametric(d, ep)

  for (lv in c("moderate", "high")) {
    lab <- paste0("control_vs_", lv)
    message("    CountBCF: ", lab)
    r <- fit_cbcf(d, ep, lv, lab)
    cbcf_ate[[paste(ep, lv)]]   <- r$ate
    cells_all[[paste(ep, lv)]]  <- r$cells
  }
}

eda_tab   <- bind_rows(eda_all)
param_tab <- bind_rows(param_all)
cbcf_tab  <- bind_rows(cbcf_ate)
cells     <- bind_rows(cells_all)

write_csv(eda_tab,   file.path(OUTDIR, "tables", "eda_group_summary.csv"))
write_csv(param_tab, file.path(OUTDIR, "tables", "parametric_ate.csv"))
write_csv(cbcf_tab,  file.path(OUTDIR, "tables", "countbcf_ate.csv"))
write_csv(cells %>% select(endpoint, contrast, sex, city, agxcat, agecat,
                           agex, age, year, pyr, dose, dose_level, pihat,
                           cate, cate_lo, cate_hi, logRR, rr),
          file.path(OUTDIR, "tables", "countbcf_cate_cells.csv"))

# ATE convergence table: parametric (adjusted) vs CountBCF, per endpoint x level
ate_compare <- bind_rows(
  param_tab %>% filter(grepl("adjusted", model)) %>%
    transmute(endpoint, contrast = paste0("control_vs_", treated),
              method = model, excess_rate_100k, excess_lo, excess_hi, rr, rr_lo, rr_hi),
  cbcf_tab %>% select(endpoint, contrast, method, excess_rate_100k,
                      excess_lo, excess_hi, rr, rr_lo, rr_hi)) %>%
  arrange(endpoint, contrast, method)
write_csv(ate_compare, file.path(OUTDIR, "tables", "ate_convergence.csv"))

# =====================================================================
# 5. CATE graphs (CountBCF)
# =====================================================================
gg_save <- function(p, file, w = 7.5, h = 5)
  ggsave(file.path(OUTDIR, "plots", file), p, width = w, height = h, dpi = 300)

plot_endpoint <- function(df, tag, nice) {
  df$Sex <- factor(sex_lab[as.character(df$sex)], levels = c("Male", "Female"))

  p1 <- ggplot(df, aes(agex, cate, colour = Sex)) +
    geom_hline(yintercept = 0, linetype = 2, colour = "grey50") +
    geom_point(alpha = .25, size = .8) +
    geom_smooth(se = TRUE, method = "loess", span = .9) +
    labs(title = paste0("CountBCF CATE by age at exposure — ", nice),
         subtitle = "Excess cases per 100,000 person-years (control vs high dose)",
         x = "Age at exposure (years)", y = "Excess rate / 100k PY")
  gg_save(p1, paste0(tag, "_cate_by_agex_sex.png"))

  p2 <- ggplot(df, aes(age, cate, colour = Sex)) +
    geom_hline(yintercept = 0, linetype = 2, colour = "grey50") +
    geom_point(alpha = .25, size = .8) +
    geom_smooth(se = TRUE, method = "loess", span = .9) +
    labs(title = paste0("CountBCF CATE by attained age — ", nice),
         subtitle = "Excess cases per 100,000 person-years (control vs high dose)",
         x = "Attained age (years)", y = "Excess rate / 100k PY")
  gg_save(p2, paste0(tag, "_cate_by_age_sex.png"))

  p3 <- ggplot(df, aes(agex, exp(logRR) - 1, colour = Sex)) +
    geom_hline(yintercept = 0, linetype = 2, colour = "grey50") +
    geom_smooth(se = TRUE, method = "loess", span = .9) +
    labs(title = paste0("CountBCF excess relative risk by age at exposure — ", nice),
         subtitle = "ERR = RR - 1 on the rate scale (control vs high dose)",
         x = "Age at exposure (years)", y = "Excess relative risk (RR - 1)")
  gg_save(p3, paste0(tag, "_err_by_agex_sex.png"))

  # subgroup posterior CATE with 95% CrI: sex x age-at-exposure category,
  # PY-weighted within subgroup
  sg <- df %>%
    mutate(agx_band = cut(agex, c(0, 20, 40, Inf),
                          labels = c("0-20", "20-40", "40+"), right = FALSE)) %>%
    group_by(Sex, agx_band) %>%
    summarise(cate = weighted.mean(cate, pyr),
              lo = weighted.mean(cate_lo, pyr),
              hi = weighted.mean(cate_hi, pyr), .groups = "drop")
  p4 <- ggplot(sg, aes(agx_band, cate, colour = Sex)) +
    geom_hline(yintercept = 0, linetype = 2, colour = "grey50") +
    geom_pointrange(aes(ymin = lo, ymax = hi),
                    position = position_dodge(.4), size = .7) +
    labs(title = paste0("CountBCF subgroup CATE — ", nice),
         subtitle = "PY-weighted posterior mean and ~95% interval",
         x = "Age-at-exposure band", y = "Excess rate / 100k PY")
  gg_save(p4, paste0(tag, "_cate_subgroups.png"), w = 7, h = 5)

  # heatmap over age-at-exposure x attained-age
  hm <- df %>%
    mutate(agx_b = cut(agex, seq(0, 80, 10)),
           age_b = cut(age,  seq(0, 110, 10))) %>%
    group_by(agx_b, age_b) %>%
    summarise(cate = weighted.mean(cate, pyr), .groups = "drop") %>%
    filter(!is.na(agx_b), !is.na(age_b))
  p5 <- ggplot(hm, aes(agx_b, age_b, fill = cate)) +
    geom_tile() +
    scale_fill_gradient2(low = "steelblue", mid = "white", high = "firebrick",
                         midpoint = 0, name = "Excess /100k PY") +
    labs(title = paste0("CountBCF CATE surface — ", nice),
         x = "Age at exposure", y = "Attained age") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  gg_save(p5, paste0(tag, "_cate_heatmap.png"), w = 7.5, h = 5.5)

  invisible(NULL)
}

# per-endpoint cells frames (re-derive from list so we keep the endpoint tag)
plot_endpoint(cells_all[["solid high"]],   "solid",   "All solid cancers")
plot_endpoint(cells_all[["stomach high"]], "stomach", "Stomach cancer")

# ATE convergence figure (NB GLM vs CountBCF)
ate_fig <- ate_compare %>% filter(grepl("NB", method)) %>%
  mutate(lvl = sub("control_vs_", "", contrast),
         m = ifelse(grepl("CountBCF", method), "CountBCF (NB)", "Parametric (NB GLM)"))
p_ate <- ggplot(ate_fig, aes(interaction(endpoint, lvl), excess_rate_100k, colour = m)) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey50") +
  geom_pointrange(aes(ymin = excess_lo, ymax = excess_hi),
                  position = position_dodge(.5)) +
  coord_flip() +
  labs(title = "ATE convergence: parametric NB vs CountBCF",
       subtitle = "Marginal excess rate (cases / 100k PY), 95% intervals",
       x = "endpoint . dose level", y = "Excess rate / 100k PY", colour = NULL)
gg_save(p_ate, "ate_convergence.png", w = 8, h = 5)

# =====================================================================
print(eda_tab,   n = 50)
print(ate_compare, n = 50)
message("\nOutputs written to: ", OUTDIR)
