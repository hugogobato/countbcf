#!/usr/bin/env Rscript
# Low-dose city comparisons: posterior diagnostics and reference-group checks.
# Primary summaries are ratios of person-year-weighted counterfactual rates
# minus one, with quantile-specific Monte Carlo errors. E-values and dose
# rescalings are defined only for each city's rate ratio, not ERR differences.
# Dose-subcategory estimates are exploratory and do not identify a gradient.
# Reference constructions split minimally exposed survivors into outer-proximal
# (2.5-3 km) and distal (3-10 km) groups; distcat is the released grouping.
# In this file, observed reference-cell distances are 2473-2997 m proximally
# and 3219-7577 m distally. Outputs: T34-T40 in analysis/tables.

source(file.path(Sys.getenv("CBCF_ROOT", "."), "Applied_Study/code/common.R"))
suppressMessages(library(countbcf))
stopifnot(requireNamespace("posterior", quietly = TRUE))

QUICK <- nzchar(Sys.getenv("CBCF_QUICK"))
NB <- if (QUICK) 200L else NBURN_MAIN
NS <- if (QUICK) 100L else NSIM_MAIN
NCORE <- n_workers(6)

d  <- read_lss()
ds <- build_ds(d, "solid", "cola02w10")

#' Add the distance-ring variables Alina's design needs.
add_rings <- function(s) {
  s %>% mutate(
    distal    = as.integer(gdist >= 3000),                 # 0 = <3 km, 1 = >=3 km
    ring      = factor(ifelse(distal == 1, ">=3km", "<3km"), levels = c("<3km", ">=3km")),
    city_ring = as.integer(interaction(city_f, ring, drop = TRUE)))
}

sub_low <- add_rings(contrast_ds(ds, "low"))

# =====================================================================
# 0.  Descriptives of the reference arm: the object of Alina's concern
# =====================================================================
ref_desc <- sub_low %>% filter(z == 0) %>% group_by(city_f, ring) %>%
  summarise(cells = n(), subjects = sum(subjects), pyr = sum(pyr),
            cases = sum(outcome), crude_rate_1e4 = 1e4 * sum(outcome) / sum(pyr),
            gdist_min = min(gdist), gdist_med = median(gdist),
            gdist_max = max(gdist), .groups = "drop")
low_desc <- sub_low %>% filter(z == 1) %>% group_by(city_f) %>%
  summarise(ring = factor("low band", levels = c("<3km", ">=3km", "low band")),
            cells = n(), subjects = sum(subjects), pyr = sum(pyr),
            cases = sum(outcome), crude_rate_1e4 = 1e4 * sum(outcome) / sum(pyr),
            gdist_min = min(gdist), gdist_med = median(gdist),
            gdist_max = max(gdist), .groups = "drop")
write_csv(bind_rows(ref_desc, low_desc), file.path(TAB, "T40_city_reference_arm.csv"))
print(as.data.frame(bind_rows(ref_desc, low_desc)), digits = 5)

# =====================================================================
# 1.  Chain-wise city statistics on BOTH footings  (Igor a, b)
# =====================================================================
# Every quantity below is computed twice: once as the person-year weighted
# mean of the per-cell ERR ("per-cell", average of ratios) and once as the
# standardised ERR ("std", ratio of PY-weighted rates).  The separation
# test is CATE_g - ATE with the ATE taken over BOTH arms, exactly as in
# common.R::subgroup_table, so the two files can be cross-checked.

city_stats <- function(r) {
  cl <- r$cells; w <- cl$pyr
  tr  <- which(cl$z == 1)
  err <- exp(r$logrr) - 1
  gi  <- function(f) intersect(which(f), tr)
  hiT <- gi(cl$city_f == "Hiroshima"); naT <- gi(cl$city_f == "Nagasaki")
  hiA <- which(cl$city_f == "Hiroshima"); naA <- which(cl$city_f == "Nagasaki")
  list(
    # exposed-arm (ATT) footing, the one the city narrative uses
    err_hi   = draw_wmean(err, w, hiT),
    err_na   = draw_wmean(err, w, naT),
    serr_hi  = draw_serr(r$logmu, r$logrr, w, hiT),
    serr_na  = draw_serr(r$logmu, r$logrr, w, naT),
    att_err  = draw_wmean(err, w, tr),
    att_serr = draw_serr(r$logmu, r$logrr, w, tr),
    # cohort (ATE) footing, used by the separation test
    ate_err  = draw_wmean(err, w),
    ate_serr = draw_serr(r$logmu, r$logrr, w),
    cate_err_hi  = draw_wmean(err, w, hiA),
    cate_err_na  = draw_wmean(err, w, naA),
    cate_serr_hi = draw_serr(r$logmu, r$logrr, w, hiA),
    cate_serr_na = draw_serr(r$logmu, r$logrr, w, naA),
    # baseline (reference-arm) rate per 10^4 PY by city, both rings
    base_hi = draw_wmean(exp(r$logmu) * 1e4, w, which(cl$city_f == "Hiroshima" & cl$z == 0)),
    base_na = draw_wmean(exp(r$logmu) * 1e4, w, which(cl$city_f == "Nagasaki"  & cl$z == 0)))
}

stamp("chain-wise fits for the city statistic")
chains <- mclapply(1:2, function(ch)
  city_stats(fit_cbcf(sub_low, label = paste0("chain", ch), ps_spec = "main",
                      nburn = NB, nsim = NS, seed = SEED + 1000L * ch)),
  mc.cores = min(2L, NCORE), mc.preschedule = FALSE)
stopifnot(!any(vapply(chains, is.null, logical(1))))

derived <- function(s) c(s, list(
  d_naga_hiro_err   = s$err_na  - s$err_hi,
  d_naga_hiro_serr  = s$serr_na - s$serr_hi,
  sep_hi_err  = s$cate_err_hi  - s$ate_err,
  sep_na_err  = s$cate_err_na  - s$ate_err,
  sep_hi_serr = s$cate_serr_hi - s$ate_serr,
  sep_na_serr = s$cate_serr_na - s$ate_serr))
chains <- lapply(chains, derived)
saveRDS(list(chains = chains, nburn = NB, nsim = NS,
             seeds = SEED + 1000L * (1:2),
             estimand = "ratio of person-year-weighted counterfactual rates minus one"),
        file.path(OUTDIR, "city_statistic_draws.rds"))

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
# 2.  E-values and per-Gy rescaling for the city contrast  (Igor c)
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

# The diagnostic repair requires only these two chains. Historical distance,
# subcategory and site fits need not be rerun when this switch is selected.
if (identical(Sys.getenv("CBCF_CITY_DIAGNOSTICS_ONLY"), "1")) {
  dir.create(file.path(OUTDIR, "logs"), recursive = TRUE, showWarnings = FALSE)
  writeLines(capture.output(sessionInfo()),
             file.path(OUTDIR, "logs", "sessionInfo_08.txt"))
  stamp("08 city diagnostics done -> ", TAB)
  quit(save = "no", status = 0)
}

# =====================================================================
# 3.  Alina's distance-stratified designs + Igor's distance refit
# =====================================================================
score_city <- function(r, spec, note = "") {
  s <- city_stats(r)
  as_tibble(c(list(spec = spec, note = note, n_cells = r$n,
                   cases = sum(r$cells$outcome),
                   pihat_min = min(r$cells$pihat), pihat_max = max(r$cells$pihat)),
              qsum(s$att_err, "att_err"), qsum(s$att_serr, "att_serr"),
              qsum(s$err_hi, "err_hiroshima"), qsum(s$err_na, "err_nagasaki"),
              qsum(s$err_na - s$err_hi, "d_naga_hiro"),
              qsum(s$serr_hi, "serr_hiroshima"), qsum(s$serr_na, "serr_nagasaki"),
              qsum(s$serr_na - s$serr_hi, "d_serr_naga_hiro"),
              qsum(s$base_hi, "base_hiroshima"), qsum(s$base_na, "base_nagasaki")))
}

# Propensity specifications that DO see the distance ring.  Distance is the
# dose-assignment mechanism, so it stays out of the primary score; here it
# enters only through the <3 / >=3 km ring, which is the baseline-structure
# axis Alina identifies, not the continuous assignment variable.
ps_ring        <- z ~ city_f * sex_f + ns(agex, 4) + ns(age, 4) + city_f:ns(agex, 2) + city_f:ring
ps_gdist_full  <- z ~ city_f * sex_f + ns(agex, 4) + ns(age, 4) + city_f:ns(agex, 2) + city_f * ns(gdist, 4)

DESIGNS <- list(
  list(spec = "Baseline (full 0-<5 mGy reference)",
       note = "primary specification, reproduced here for comparison",
       filter = NULL, xc = c("agex", "age", "sex", "city"), ps = "main"),
  list(spec = "Reference = outer proximal only (2.5-3.0 km)",
       note = "RERF construction recommended by A. Brenner",
       filter = quote(z == 1 | distal == 0), xc = c("agex", "age", "sex", "city"), ps = "main"),
  list(spec = "Reference = distal only (3-10 km)",
       note = "complement of the RERF construction",
       filter = quote(z == 1 | distal == 1), xc = c("agex", "age", "sex", "city"), ps = "main"),
  list(spec = "Distal indicator in the prognostic forest",
       note = "full reference arm, ring in x_control",
       filter = NULL, xc = c("agex", "age", "sex", "city", "distal"), ps = "main"),
  list(spec = "City x distance in PS and forest",
       note = "city:ring in the score AND city_ring in x_control",
       filter = NULL, xc = c("agex", "age", "sex", "city", "city_ring"), ps = "ring"),
  list(spec = "City x distance in the PS only",
       note = "city:ring added to the primary logistic score",
       filter = NULL, xc = c("agex", "age", "sex", "city"), ps = "ring"),
  list(spec = "Continuous ground distance in PS and forest",
       note = "conditions on the assignment variable; positivity check",
       filter = NULL, xc = c("agex", "age", "sex", "city", "gdist"), ps = "gdist"))

run_design <- function(k) {
  g <- DESIGNS[[k]]
  s <- sub_low
  if (!is.null(g$filter)) s <- s[eval(g$filter, s), , drop = FALSE]
  ph <- if (g$ps == "main") fit_ps(s, "main") else {
    f <- if (g$ps == "ring") ps_ring else ps_gdist_full
    as.numeric(predict(glm(f, binomial, data = s), type = "response"))
  }
  r <- tryCatch(fit_cbcf(s, label = g$spec, pihat = ph, x_control_vars = g$xc,
                         nburn = NB, nsim = NS),
                error = function(e) { stamp("FAILED ", g$spec, ": ", conditionMessage(e)); NULL })
  if (is.null(r)) return(NULL)
  out <- score_city(r, g$spec, g$note); rm(r); gc(verbose = FALSE); out
}

stamp("distance-stratified designs: ", length(DESIGNS), " fits")
D <- bind_rows(mclapply(seq_along(DESIGNS), run_design,
                        mc.cores = NCORE, mc.preschedule = FALSE)) %>%
  mutate(separates_percell = (d_naga_hiro_lo > 0) | (d_naga_hiro_hi < 0),
         separates_std     = (d_serr_naga_hiro_lo > 0) | (d_serr_naga_hiro_hi < 0))
write_csv(D, file.path(TAB, "T36_city_distance_sensitivity.csv"))
print(as.data.frame(D %>% dplyr::select(spec, n_cells, pihat_min, pihat_max,
  serr_hiroshima, serr_nagasaki, d_serr_naga_hiro, d_serr_naga_hiro_lo,
  d_serr_naga_hiro_hi, separates_std, d_naga_hiro, separates_percell)), digits = 3)

# =====================================================================
# 4.  Gradient or step?  The low band split by dose  (Igor d)
# =====================================================================
# A radiation effect must grow with dose inside the band.  A population
# offset attached to being in the exposed ring at all shows up as a step at
# the 5 mGy boundary and no gradient above it.
SUBBANDS <- tribble(
  ~name,        ~lo,    ~hi,
  "5-<20 mGy",  0.005,  0.020,
  "20-<50 mGy", 0.020,  0.050,
  "50-<100 mGy",0.050,  0.100)

run_subband <- function(k) {
  b <- SUBBANDS[k, ]
  s <- sub_low %>% filter(z == 0 | (dose >= b$lo & dose < b$hi))
  r <- tryCatch(fit_cbcf(s, label = b$name, ps_spec = "main", nburn = NB, nsim = NS),
                error = function(e) { stamp("FAILED ", b$name, ": ", conditionMessage(e)); NULL })
  if (is.null(r)) return(NULL)
  out <- score_city(r, b$name, "low-band dose split") %>%
    mutate(dose_contrast_Gy = dose_contrast(s),
           dose_contrast_Gy_hi = dose_contrast(s, s$city_f == "Hiroshima"),
           dose_contrast_Gy_na = dose_contrast(s, s$city_f == "Nagasaki"))
  rm(r); gc(verbose = FALSE); out
}
stamp("low-band dose split: ", nrow(SUBBANDS), " fits")
G <- bind_rows(mclapply(seq_len(nrow(SUBBANDS)), run_subband,
                        mc.cores = NCORE, mc.preschedule = FALSE))
if (nrow(G)) G <- G %>% mutate(
  serr_per_Gy_hiroshima = serr_hiroshima / dose_contrast_Gy_hi,
  serr_per_Gy_nagasaki  = serr_nagasaki  / dose_contrast_Gy_na,
  separates_std = (d_serr_naga_hiro_lo > 0) | (d_serr_naga_hiro_hi < 0))
write_csv(G, file.path(TAB, "T37_lowband_dose_gradient.csv"))
print(as.data.frame(G %>% dplyr::select(spec, n_cells, dose_contrast_Gy,
  serr_hiroshima, serr_nagasaki, d_serr_naga_hiro, d_serr_naga_hiro_lo,
  d_serr_naga_hiro_hi, serr_per_Gy_hiroshima, serr_per_Gy_nagasaki)), digits = 3)

# =====================================================================
# 5.  Site composition on the low-band city contrast  (Igor f)
# =====================================================================
SITES <- tribble(
  ~name,                          ~outcome,
  "All solid (baseline)",         "solid",
  "Excluding breast and thyroid", "solid_no_bt",
  "Excluding non-melanoma skin",  "solid_no_skin",
  "Excluding stomach",            "solid_no_stom",
  "Non-sex-specific solid only",  "solid_ngs")

run_site <- function(k) {
  b <- SITES[k, ]
  s <- add_rings(contrast_ds(build_ds(d, b$outcome, "cola02w10"), "low"))
  r <- tryCatch(fit_cbcf(s, label = b$name, ps_spec = "main", nburn = NB, nsim = NS),
                error = function(e) { stamp("FAILED ", b$name, ": ", conditionMessage(e)); NULL })
  if (is.null(r)) return(NULL)
  out <- score_city(r, b$name, "site composition"); rm(r); gc(verbose = FALSE); out
}
stamp("site-composition perturbations: ", nrow(SITES), " fits")
S <- bind_rows(mclapply(seq_len(nrow(SITES)), run_site,
                        mc.cores = NCORE, mc.preschedule = FALSE))
if (nrow(S)) S <- S %>% mutate(
  separates_std = (d_serr_naga_hiro_lo > 0) | (d_serr_naga_hiro_hi < 0),
  separates_percell = (d_naga_hiro_lo > 0) | (d_naga_hiro_hi < 0))
write_csv(S, file.path(TAB, "T38_city_site_composition.csv"))
print(as.data.frame(S %>% dplyr::select(spec, cases, serr_hiroshima, serr_nagasaki,
  d_serr_naga_hiro, d_serr_naga_hiro_lo, d_serr_naga_hiro_hi, separates_std)), digits = 3)

# =====================================================================
# 6.  Parametric analogue of Alina's design  (cheap, GLM)
# =====================================================================
# Same Poisson model as 07_lowdose_city.R, refitted under three reference
# constructions, so the city x dose interaction can be read off a model
# that RERF would recognise.
ic <- d %>% filter(dcatn > 0, band %in% BANDS) %>%
  mutate(outcome = solid, sex_f = factor(sex, 1:2, c("Male", "Female")),
         city_f = factor(city, 1:2, c("Hiroshima", "Nagasaki")),
         ring = factor(ifelse(gdist >= 3000, ">=3km", "<3km"), levels = c("<3km", ">=3km")),
         band = droplevels(factor(as.character(band), levels = BANDS)))

fit_glm_city <- function(dat, extra = NULL, tag) {
  f0 <- outcome ~ band + sex_f * (ns(log(age), 4) + ns(agex, 3)) + city_f + offset(log(pyr))
  if (!is.null(extra)) f0 <- update(f0, extra)
  m0 <- glm(f0, poisson("log"), dat)
  m1 <- update(m0, . ~ . + band:city_f)
  cf <- coef(summary(m1)); lr <- anova(m0, m1, test = "LRT")
  tibble(model = tag, term = rownames(cf)) %>% filter(grepl("band.*city", term)) %>%
    mutate(RR_ratio = exp(cf[term, 1]),
           lo = exp(cf[term, 1] - 1.96 * cf[term, 2]),
           hi = exp(cf[term, 1] + 1.96 * cf[term, 2]),
           p = cf[term, 4], lrt_p = lr$`Pr(>Chi)`[2],
           n_cells = nrow(dat), cases = sum(dat$outcome))
}

P <- bind_rows(
  fit_glm_city(ic, NULL, "Primary (full 0-<5 mGy reference)"),
  fit_glm_city(ic, . ~ . + city_f:ring, "City-specific distal baseline term"),
  fit_glm_city(ic %>% filter(band != "reference" | ring == "<3km"),
               NULL, "Reference = outer proximal only"),
  fit_glm_city(ic %>% filter(band != "reference" | ring == ">=3km"),
               NULL, "Reference = distal only"))
write_csv(P, file.path(TAB, "T39_parametric_city_reference_groups.csv"))
print(as.data.frame(P), digits = 4)

# Reference-arm baseline rates, directly, adjusted for age and sex (Igor g)
base_glm <- glm(outcome ~ city_f * ring + sex_f * (ns(log(age), 4) + ns(agex, 3)) +
                  offset(log(pyr)), poisson("log"),
                ic %>% filter(band == "reference"))
nd <- ic %>% filter(band == "reference") %>% mutate(pyr = 1)
base_pred <- bind_rows(lapply(list(c("Hiroshima", "<3km"), c("Hiroshima", ">=3km"),
                                   c("Nagasaki", "<3km"), c("Nagasaki", ">=3km")), function(k) {
  nd2 <- nd %>% mutate(city_f = factor(k[1], levels = levels(ic$city_f)),
                       ring   = factor(k[2], levels = levels(ic$ring)))
  w <- ic %>% filter(band == "reference") %>% pull(pyr)
  tibble(city = k[1], ring = k[2],
         adj_rate_1e4 = 1e4 * weighted.mean(predict(base_glm, nd2, type = "response"), w))
}))
write_csv(base_pred, file.path(TAB, "T40b_reference_arm_adjusted_rates.csv"))
print(as.data.frame(base_pred), digits = 5)

writeLines(capture.output(sessionInfo()),
           file.path(OUTDIR, "logs", "sessionInfo_08.txt"))
stamp("08 done -> ", TAB)
