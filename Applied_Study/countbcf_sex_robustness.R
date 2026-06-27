#!/usr/bin/env Rscript
# =====================================================================
# Phase 2 (B-1..B-5) + Phase 3 (B-6, B-7) for the Applied Study revision
# (response to Igor & Eric, 2026-06).
#
#   B-1  All-solid construction check (breast/thyroid in? nmskin out?)
#   B-2  Breast & thyroid fit with their OWN organ doses (recover sex ERR)
#   B-3  All-solid WITH vs WITHOUT breast+thyroid, still colon-dosed
#   B-4  Design-factor sensitivity of the solid sex ERR
#        (i) 2 Gy vs 4 Gy cap; (ii) dose-band cutpoints; (iii) PY-weighted vs not
#   B-5  Packaging: this script + the CSVs it writes ARE the deliverable.
#   B-6  E-values for the four ATE contrasts (no new fit).
#   B-7  Site-stratified age-at-exposure shape (Little-2009 aggregation check).
#
# Modelling spec is IDENTICAL to countbcf_applied_analysis.R:
#   prognostic covars  agex, age, sex, city, year, time
#   moderating  covars agex, age, sex
#   count_model = "nb", offset = person-years, include_pihat = "control"
#   propensity from demographics only (distance EXCLUDED).
#
# Writes CSV tables to countbcf_applied_analysis/tables/ and prints a
# console summary used to finalise A-1a/A-1b and A-5e in the draft.
# =====================================================================

suppressMessages({
  library(MASS); library(dplyr); library(readr); library(tibble); library(countbcf)
})
set.seed(20260622)

BASE <- "Applied_Study/RERF_data"
OUT  <- "Applied_Study/countbcf_applied_analysis"
TAB  <- file.path(OUT, "tables")
dir.create(TAB, showWarnings = FALSE, recursive = TRUE)
stamp <- function(...) message(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ...)

solid <- read_csv(file.path(BASE, "lssinc07.csv"), show_col_types = FALSE)
names(solid) <- trimws(names(solid))
solid$solid_excl_bt <- pmax(solid$solid - solid$breast - solid$thyroid, 0)

NBURN <- 1000L; NSIM <- 1000L
DEF_CAP <- 2.0; DEF_LOW <- 0.1; DEF_MOD <- 0.5

# ---------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------
build_ds <- function(data, outcome_col, dose_col,
                     cap = DEF_CAP, low = DEF_LOW, mod = DEF_MOD, sexes = c(1, 2)) {
  data %>%
    mutate(dose = ifelse(.data[[dose_col]] < 0, NA_real_, .data[[dose_col]]),
           outcome = .data[[outcome_col]]) %>%
    filter(un4gy == 1, sex %in% sexes, !is.na(dose), dose >= 0, dose <= cap,
           !is.na(outcome), !is.na(pyr), pyr > 0) %>%
    mutate(dose_level = factor(case_when(dose <  low ~ "control",
                                         dose <  mod ~ "moderate",
                                         TRUE        ~ "high"),
                               levels = c("control", "moderate", "high")),
           sex_f = factor(sex), city_f = factor(city)) %>%
    select(outcome, dose, pyr, dose_level, sex, city, sex_f, city_f,
           agex, age, year, time)
}

# CountBCF fit, exactly as in the main analysis script.
fit_cbcf <- function(d, treated = "high", label = "") {
  sub <- d %>% filter(dose_level %in% c("control", treated)) %>% droplevels()
  z   <- as.integer(sub$dose_level == treated)
  # drop single-level factors (e.g. sex for a female-only endpoint like breast)
  nuniq <- function(v) length(unique(sub[[v]]))
  pterms <- c("agex", "age")
  if (nlevels(factor(sub$sex))  > 1) pterms <- c(pterms, "sex_f")
  if (nlevels(factor(sub$city)) > 1) pterms <- c(pterms, "city_f")
  ps    <- glm(reformulate(pterms, response = "z"), binomial, data = sub)
  pihat <- as.numeric(predict(ps, type = "response"))
  cterms <- c("agex", "age", "sex", "city", "year", "time")
  cterms <- cterms[vapply(cterms, nuniq, integer(1)) > 1]
  mterms <- c("agex", "age", "sex")
  mterms <- mterms[vapply(mterms, nuniq, integer(1)) > 1]
  xc <- model.matrix(reformulate(cterms), sub)[, -1, drop = FALSE]
  xm <- model.matrix(reformulate(mterms), sub)[, -1, drop = FALSE]
  fit <- countbcf(y = sub$outcome, z = z, x_control = xc, x_moderate = xm,
                  pihat = pihat, offset = sub$pyr, count_model = "nb",
                  include_pihat = "control", nburn = NBURN, nsim = NSIM)
  rate0 <- exp(fit$mu_f_post); rate1 <- exp(fit$mu_f_post + fit$tau_f_post)
  cate  <- (rate1 - rate0) * 1e5; logRR <- fit$tau_f_post
  sub$cate  <- colMeans(cate)
  sub$rr    <- colMeans(exp(logRR))
  sub$logRR <- colMeans(logRR)
  w <- sub$pyr / sum(sub$pyr)
  ate_draws <- as.numeric(cate %*% w)
  rr_draws  <- as.numeric((rate1 %*% sub$pyr) / (rate0 %*% sub$pyr))
  list(cells = sub,
       ate = tibble(label = label, n_cells = nrow(sub),
                    cases_treated = sum(sub$outcome[sub$dose_level == treated]),
                    excess = mean(ate_draws),
                    excess_lo = quantile(ate_draws, .025),
                    excess_hi = quantile(ate_draws, .975),
                    rr = mean(rr_draws),
                    rr_lo = quantile(rr_draws, .025),
                    rr_hi = quantile(rr_draws, .975)))
}

# PY-weighted (and unweighted) ERR/EAR by sex from CountBCF cells.
sex_err <- function(cells, label) {
  cells %>% group_by(sex) %>%
    summarise(err_pyw = weighted.mean(rr - 1, pyr),
              err_unw = mean(rr - 1),
              ear_pyw = weighted.mean(cate, pyr),
              cases   = sum(outcome), pyr = sum(pyr), .groups = "drop") %>%
    mutate(label = label, sexlab = ifelse(sex == 1, "male", "female"))
}

# Sex-stratified NB-GLM ERR (high vs control) with own organ dose -- the
# direct, robust comparator for B-2.
glm_err <- function(d, label) {
  res <- lapply(sort(unique(d$sex)), function(sx) {
    ds <- d %>% filter(sex == sx, dose_level %in% c("control", "high")) %>% droplevels()
    nh <- sum(ds$outcome[ds$dose_level == "high"])
    if (nh < 3 || nlevels(droplevels(ds$dose_level)) < 2) return(NULL)
    f <- outcome ~ dose_level + city_f + agex + age + offset(log(pyr))
    if (nlevels(factor(ds$city)) < 2) f <- update(f, . ~ . - city_f)
    m <- tryCatch(glm.nb(f, ds), error = function(e) NULL)
    if (is.null(m) || !"dose_levelhigh" %in% names(coef(m))) return(NULL)
    b <- coef(m)["dose_levelhigh"]; se <- sqrt(vcov(m)["dose_levelhigh", "dose_levelhigh"])
    tibble(label = label, sex = sx, sexlab = ifelse(sx == 1, "male", "female"),
           err = exp(b) - 1, err_lo = exp(b - 1.96 * se) - 1,
           err_hi = exp(b + 1.96 * se) - 1, cases_high = nh)
  })
  bind_rows(res)
}

sexratio <- function(se_tab, label) {
  m <- se_tab$err_pyw[se_tab$sexlab == "male"]
  f <- se_tab$err_pyw[se_tab$sexlab == "female"]
  mu <- se_tab$err_unw[se_tab$sexlab == "male"]
  fu <- se_tab$err_unw[se_tab$sexlab == "female"]
  tibble(label = label,
         err_male_pyw = m, err_female_pyw = f,
         ratio_fm_pyw = ifelse(length(m) && length(f), f / m, NA_real_),
         err_male_unw = mu, err_female_unw = fu,
         ratio_fm_unw = ifelse(length(mu) && length(fu), fu / mu, NA_real_))
}

# =====================================================================
# B-1  All-solid construction table
# =====================================================================
stamp("B-1 all-solid construction")
S <- function(x) sum(solid[[x]], na.rm = TRUE)
disj <- c("oralca", "digestca", "respca", "femgenca", "malgenca", "urinca")  # parent=children EXACT
base <- rowSums(solid[, disj], na.rm = TRUE)
named_incl_bt <- base + solid$breast + solid$thyroid + solid$cnsca +
                 solid$bone + solid$connect + solid$thymus + solid$othsol
b1 <- tibble(
  group = c(disj, "breast", "thyroid", "cnsca",
            "bone", "connect", "thymus", "othsol",
            "(named solid incl. breast+thyroid)", "solid (RERF all-solid total)",
            "nmskin (non-melanoma skin)"),
  cases = c(sapply(disj, S), S("breast"), S("thyroid"), S("cnsca"),
            S("bone"), S("connect"), S("thymus"), S("othsol"),
            sum(named_incl_bt), S("solid"), S("nmskin")),
  in_all_solid = c(rep("yes (exact group)", length(disj)),
                   "yes", "yes", "yes", "yes", "yes", "yes", "yes",
                   "—", "—", "NO (RERF convention)"))
b1_tests <- tibble(
  test = c("cells where solid < (named solid incl. breast+thyroid)",
           "share of solid accounted by named solid sites incl. breast+thyroid",
           "residual (grouped 'other & ill-defined' slack)"),
  value = c(as.character(sum(solid$solid < named_incl_bt, na.rm = TRUE)),
            sprintf("%.1f%%", 100 * sum(named_incl_bt) / S("solid")),
            as.character(S("solid") - sum(named_incl_bt))))
write_csv(b1, file.path(TAB, "B1_allsolid_construction.csv"))
write_csv(b1_tests, file.path(TAB, "B1_allsolid_tests.csv"))
print(b1); print(b1_tests)

# =====================================================================
# B-2  Breast & thyroid with their OWN organ doses
# =====================================================================
b2_eda <- list(); b2_cbcf <- list(); b2_glm <- list()
b2_specs <- tibble::tribble(
  ~name,     ~outcome,  ~dose,        ~sexes,
  "breast",  "breast",  "brea02w10",  "2",       # female breast
  "thyroid", "thyroid", "thya02w10",  "1,2")
for (i in seq_len(nrow(b2_specs))) {
  nm <- b2_specs$name[i]
  sx <- as.integer(strsplit(b2_specs$sexes[i], ",")[[1]])
  d  <- build_ds(solid, b2_specs$outcome[i], b2_specs$dose[i], sexes = sx)
  ed <- d %>% group_by(dose_level) %>%
    summarise(cells = n(), cases = sum(outcome), pyr = sum(pyr), .groups = "drop") %>%
    mutate(endpoint = nm)
  b2_eda[[nm]] <- ed
  stamp("B-2 ", nm, "  high-dose cases = ",
        ed$cases[ed$dose_level == "high"], " | GLM ERR")
  b2_glm[[nm]] <- glm_err(d, nm)
  stamp("B-2 ", nm, "  CountBCF")
  r <- tryCatch(fit_cbcf(d, "high", nm), error = function(e) { stamp("  cbcf failed: ", conditionMessage(e)); NULL })
  if (!is.null(r)) b2_cbcf[[nm]] <- sex_err(r$cells, nm)
}
b2_eda_tab  <- bind_rows(b2_eda)
b2_glm_tab  <- bind_rows(b2_glm)
b2_cbcf_tab <- bind_rows(b2_cbcf)
write_csv(b2_eda_tab,  file.path(TAB, "B2_breast_thyroid_eda.csv"))
write_csv(b2_glm_tab,  file.path(TAB, "B2_breast_thyroid_glm_err.csv"))
write_csv(b2_cbcf_tab, file.path(TAB, "B2_breast_thyroid_cbcf_sexerr.csv"))
print(b2_eda_tab); print(b2_glm_tab); print(b2_cbcf_tab)

# =====================================================================
# B-3  All-solid WITH vs WITHOUT breast+thyroid, still colon-dosed
# =====================================================================
b3_cells <- list(); b3_ate <- list(); b3_sex <- list()
b3_specs <- tibble::tribble(
  ~name,                  ~outcome,
  "solid_full",           "solid",
  "solid_excl_breast_thy","solid_excl_bt")
for (i in seq_len(nrow(b3_specs))) {
  nm <- b3_specs$name[i]
  stamp("B-3 ", nm, " (colon-dosed)")
  d <- build_ds(solid, b3_specs$outcome[i], "cola02w10")
  r <- tryCatch(fit_cbcf(d, "high", nm),
                error = function(e) { stamp("  B-3 ", nm, " failed: ", conditionMessage(e)); NULL })
  if (is.null(r)) next
  b3_cells[[nm]] <- r$cells; b3_ate[[nm]] <- r$ate
  b3_sex[[nm]]   <- sex_err(r$cells, nm)
}
b3_sex_tab   <- bind_rows(b3_sex)
b3_ratio_tab <- bind_rows(lapply(names(b3_sex), function(nm) sexratio(b3_sex[[nm]], nm)))
write_csv(b3_sex_tab,   file.path(TAB, "B3_solid_excl_bt_sexerr.csv"))
write_csv(b3_ratio_tab, file.path(TAB, "B3_solid_excl_bt_sexratio.csv"))
write_csv(bind_rows(b3_ate), file.path(TAB, "B3_solid_excl_bt_ate.csv"))
print(b3_ratio_tab)

# =====================================================================
# B-4  Design-factor sensitivity of the solid sex ERR (colon-dosed)
# =====================================================================
# Re-use the B-3 full-solid fit as the baseline spec (cap 2, bands .1/.5).
b4_specs <- tibble::tribble(
  ~name,            ~cap, ~low,  ~mod,
  "cap4Gy",          4.0, 0.10,  0.50,
  "bands_hi1Gy",     2.0, 0.10,  1.00,
  "bands_low.05",    2.0, 0.05,  0.50)
b4_sex <- list(); b4_ratio <- list()
# baseline (from B-3)
b4_sex[["baseline_cap2_b.1.5"]]   <- b3_sex[["solid_full"]] %>% mutate(label = "baseline_cap2_b.1.5")
b4_ratio[["baseline_cap2_b.1.5"]] <- sexratio(b3_sex[["solid_full"]], "baseline_cap2_b.1.5")
for (i in seq_len(nrow(b4_specs))) {
  nm <- b4_specs$name[i]
  stamp("B-4 ", nm)
  d <- build_ds(solid, "solid", "cola02w10",
                cap = b4_specs$cap[i], low = b4_specs$low[i], mod = b4_specs$mod[i])
  r <- tryCatch(fit_cbcf(d, "high", nm),
                error = function(e) { stamp("  B-4 ", nm, " failed: ", conditionMessage(e)); NULL })
  if (is.null(r)) next
  se <- sex_err(r$cells, nm)
  b4_sex[[nm]] <- se; b4_ratio[[nm]] <- sexratio(se, nm)
}
b4_ratio_tab <- bind_rows(b4_ratio)
write_csv(bind_rows(b4_sex), file.path(TAB, "B4_design_sensitivity_sexerr.csv"))
write_csv(b4_ratio_tab,      file.path(TAB, "B4_design_sensitivity_sexratio.csv"))
print(b4_ratio_tab)

# =====================================================================
# B-6  E-values for the four ATE contrasts (uses existing CountBCF RR/CI)
# =====================================================================
stamp("B-6 E-values")
evalue_pt <- function(rr) { if (rr < 1) rr <- 1 / rr; rr + sqrt(rr * (rr - 1)) }
evalue_ci <- function(lo, hi) {                 # E-value for the CI limit nearest the null
  if (lo <= 1 && hi >= 1) return(1)             # CI crosses null -> E-value = 1
  b <- if (hi < 1) hi else lo
  if (b < 1) b <- 1 / b
  b + sqrt(b * (b - 1))
}
ate <- read_csv(file.path(TAB, "countbcf_ate.csv"), show_col_types = FALSE)
b6 <- ate %>% transmute(endpoint, contrast, rr, rr_lo, rr_hi,
                        evalue_point = vapply(rr, evalue_pt, numeric(1)),
                        evalue_ci    = mapply(evalue_ci, rr_lo, rr_hi))
write_csv(b6, file.path(TAB, "B6_evalues.csv"))
print(b6)

# =====================================================================
# B-7  Site-stratified age-at-exposure shape (Little-2009 aggregation check)
# =====================================================================
agex_band <- function(cells, label) {
  cells %>%
    mutate(band = cut(agex, c(0, 10, 20, 30, 40, 50, 60, Inf),
                      labels = c("<10","10-20","20-30","30-40","40-50","50-60","60+"),
                      right = FALSE)) %>%
    group_by(band) %>%
    summarise(err_pyw = weighted.mean(rr - 1, pyr), pyr = sum(pyr),
              cases = sum(outcome), .groups = "drop") %>%
    mutate(endpoint = label)
}
b7_specs <- tibble::tribble(
  ~name,    ~outcome,  ~dose,
  "solid",  "solid",   "cola02w10",
  "stomach","stomach", "stoa02w10",
  "lung",   "lung",    "luna02w10",
  "colon",  "colon",   "cola02w10")
b7 <- list()
for (i in seq_len(nrow(b7_specs))) {
  nm <- b7_specs$name[i]
  stamp("B-7 ", nm)
  d <- build_ds(solid, b7_specs$outcome[i], b7_specs$dose[i])
  r <- tryCatch(fit_cbcf(d, "high", nm), error = function(e) { stamp("  failed: ", conditionMessage(e)); NULL })
  if (!is.null(r)) b7[[nm]] <- agex_band(r$cells, nm)
}
b7_tab <- bind_rows(b7)
write_csv(b7_tab, file.path(TAB, "B7_site_stratified_agex.csv"))
print(b7_tab, n = 40)

stamp("DONE. Tables in ", TAB)
