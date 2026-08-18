#!/usr/bin/env Rscript
# =====================================================================
# 03  Sensitivity analyses.
#
#  A. Regularisation of the MODERATING forest (the prior BCF deliberately
#     tightens to suppress spurious heterogeneity):
#       * number of moderating trees  -- in this implementation the leaf
#         prior on tau has variance a0^2 * m_tau / m_mu, so m_tau IS the
#         prior-scale knob (see note in the FINDINGS file);
#       * the CGM depth penalty (base_moderate, power_moderate);
#       * the overall leaf concentration a0.
#
#  B. Data perturbations:
#       * cancer-site composition of the all-solid aggregate;
#       * NIC survivors included in the control arm (with an NIC indicator
#         in the prognostic forest) vs excluded;
#       * reference-band and dose-cap definitions;
#       * organ-specific vs colon-dose banding for site endpoints;
#       * Poisson vs negative-binomial likelihood.
#
#  Every fit is scored on the same three quantities the main analysis
#  features: the ATT on the ERR scale, the female:male ERR difference, and
#  the age-at-exposure gradient of the ERR.
#
#  Each of those is reported on BOTH footings, because the manuscript quotes
#  the STANDARDISED ERR (ratio of person-year weighted counterfactual rates)
#  everywhere: `*_serr_*` columns are the standardised versions and are the
#  ones the tables use; the plain `d_female_male` / `d_young_old` columns are
#  the person-year weighted mean of the per-cell ERR (average of ratios) and
#  are kept only so the two footings can be cross-checked.
#
#  Labelling note (I. Shuryak, August 2026 review): ntree_moderate and the
#  (base, power) depth prior are MODERATING-forest knobs, whereas a0 is the
#  leaf concentration SHARED by both forests (R/countbcf.R:90), with the
#  implied treatment-effect prior scale a0 * sqrt(ntree_moderate /
#  ntree_control).  An a0 sweep therefore moves the prognostic surface too and
#  is labelled as such rather than as a moderating-forest-only sweep.
# =====================================================================

source(file.path(Sys.getenv("CBCF_ROOT", "."), "Applied_Study/code/common.R"))
suppressMessages(library(countbcf))

d <- read_lss()

# ---- scoring: the three featured quantities -------------------------
score_fit <- function(r, spec) {
  cl <- add_bands(r$cells); w <- cl$pyr; tr <- which(cl$z == 1)
  err <- exp(r$logrr) - 1
  ear <- ear_mat(r$logmu, r$logrr)
  att <- draw_wmean(err, w, tr)
  ear2 <- draw_wmean(ear, w, tr)
  gi <- function(f) intersect(which(f), tr)
  m <- draw_wmean(err, w, gi(cl$sex_f == "Male"))
  f <- draw_wmean(err, w, gi(cl$sex_f == "Female"))
  y <- draw_wmean(err, w, gi(cl$agex < 20))
  o <- draw_wmean(err, w, gi(cl$agex >= 40))
  hi <- draw_wmean(err, w, gi(cl$city_f == "Hiroshima"))
  na <- draw_wmean(err, w, gi(cl$city_f == "Nagasaki"))
  # ... and the same four contrasts on the standardised footing, which is what
  # the manuscript tabulates
  sg <- function(k) draw_serr(r$logmu, r$logrr, w, gi(k))
  sm <- sg(cl$sex_f == "Male");      sf <- sg(cl$sex_f == "Female")
  sy <- sg(cl$agex < 20);            so <- sg(cl$agex >= 40)
  shi <- sg(cl$city_f == "Hiroshima"); sna <- sg(cl$city_f == "Nagasaki")
  as_tibble(c(list(spec = spec, n_cells = nrow(cl), cases = sum(cl$outcome),
                   kappa = mean(r$kappa)),
              qsum(att, "att_err"), qsum(ear2, "att_ear"),
              qsum(draw_serr(r$logmu, r$logrr, w, tr), "att_serr"),
              qsum(sf - sm, "d_serr_female_male"),
              qsum(sy - so, "d_serr_young_old"),
              qsum(sna - shi, "d_serr_naga_hiro"),
              qsum(sf, "serr_female"), qsum(sm, "serr_male"),
              qsum(f - m, "d_female_male"), qsum(f / m, "ratio_female_male"),
              qsum(y - o, "d_young_old"), qsum(hi - na, "d_hiro_naga")))
}

run_spec <- function(sub, spec, ...) {
  r <- tryCatch(fit_cbcf(sub, label = spec, nburn = NBURN_SENS, nsim = NSIM_SENS, ...),
                error = function(e) { stamp("FAILED ", spec, ": ", conditionMessage(e)); NULL })
  if (is.null(r)) return(NULL)
  out <- score_fit(r, spec)
  rm(r); gc(verbose = FALSE)
  out
}

# The featured contrast for all sensitivity work: all-solid, high band.
ds_main  <- build_ds(d, "solid", "cola02w10")
sub_main <- contrast_ds(ds_main, "high")

# =====================================================================
# A. Moderating-forest regularisation
# =====================================================================
gridA <- bind_rows(
  tibble(spec = paste0("ntree_moderate=", c(10, 25, 50, 100, 200)),
         ntree_moderate = c(10, 25, 50, 100, 200),
         base_moderate = 0.25, power_moderate = 3, a0mult = 1),
  tibble(spec = c("tau tree prior (0.05,4) tight", "tau tree prior (0.25,3) default",
                  "tau tree prior (0.5,2) loose", "tau tree prior (0.95,2) mu-like"),
         ntree_moderate = 50,
         base_moderate = c(0.05, 0.25, 0.5, 0.95),
         power_moderate = c(4, 3, 2, 2), a0mult = 1),
  # a0 is SHARED by the mu and tau forests, so this sweep is not a
  # moderating-forest-only sweep; the label says so.
  tibble(spec = paste0("a0 x ", c(0.5, 1, 2), " (both forests)"),
         ntree_moderate = 50,
         base_moderate = 0.25, power_moderate = 3, a0mult = c(0.5, 1, 2)))

# default a0 (the value countbcf would pick) so the multiplier is meaningful
a0_default <- 0.5 * (log(max(quantile(sub_main$outcome, 0.95), 1)) -
                     mean(log(sub_main$pyr)))
a0_default <- if (!is.finite(a0_default) || a0_default <= 0) 1 else a0_default
stamp("default a0 for the featured contrast = ", signif(a0_default, 4))

stamp("A. moderating-forest regularisation: ", nrow(gridA), " fits")
resA <- mclapply(seq_len(nrow(gridA)), function(k) {
  g <- gridA[k, ]
  run_spec(sub_main, g$spec, ntree_moderate = g$ntree_moderate,
           base_moderate = g$base_moderate, power_moderate = g$power_moderate,
           a0 = a0_default * g$a0mult)
}, mc.cores = n_workers(6), mc.preschedule = FALSE)
A <- bind_rows(resA) %>% mutate(family = "A. tau regularisation", .before = 1)
# The (0.25,3) depth prior and a0 x 1 ARE the default specification, so three
# of the twelve rows are the same fit reported once per sweep; flag it here so
# the table caption can say "twelve rows, ten distinct fits".
A$is_default_repeat <- as.vector(duplicated(round(cbind(
  A$att_serr, A$d_serr_young_old, A$d_serr_female_male), 8)))
write_csv(A, file.path(TAB, "T20_sensitivity_tau_prior.csv"))
print(as.data.frame(A %>% dplyr::select(spec, att_serr, att_serr_lo, att_serr_hi,
  d_serr_female_male, d_serr_female_male_lo, d_serr_female_male_hi,
  d_serr_young_old, d_serr_young_old_lo, d_serr_young_old_hi)), digits = 3)
cat("rows that repeat the default specification: ", sum(A$is_default_repeat), "\n")

# =====================================================================
# B. Data perturbations
# =====================================================================
gridB <- tribble(
  ~spec,                                  ~outcome,        ~dose_var,   ~band_var,     ~include_nic, ~ref_bands,           ~cap,        ~model,
  "baseline (solid, colon bands, NB)",    "solid",         "cola02w10", "colon",       FALSE,        "reference",          "high",      "nb",
  "Poisson likelihood",                   "solid",         "cola02w10", "colon",       FALSE,        "reference",          "high",      "poisson",
  "NIC included in control arm",          "solid",         "cola02w10", "colon",       TRUE,         "reference",          "high",      "nb",
  "reference = 0-<100 mGy",               "solid",         "cola02w10", "colon",       FALSE,        "reference,low",      "high",      "nb",
  "high band includes >=2 Gy",            "solid",         "cola02w10", "colon",       FALSE,        "reference",          "high_plus", "nb",
  "solid excl. breast + thyroid",         "solid_no_bt",   "cola02w10", "colon",       FALSE,        "reference",          "high",      "nb",
  "solid excl. non-melanoma skin",        "solid_no_skin", "cola02w10", "colon",       FALSE,        "reference",          "high",      "nb",
  "solid excl. stomach",                  "solid_no_stom", "cola02w10", "colon",       FALSE,        "reference",          "high",      "nb",
  "non-gender-specific solid only",       "solid_ngs",     "cola02w10", "colon",       FALSE,        "reference",          "high",      "nb",
  "stomach, colon-dose bands",            "stomach",       "stoa02w10", "colon",       FALSE,        "reference",          "high",      "nb",
  "stomach, stomach-dose bands",          "stomach",       "stoa02w10", "stoa02w10",   FALSE,        "reference",          "high",      "nb",
  "lung, lung-dose bands",                "lung",          "luna02w10", "luna02w10",   FALSE,        "reference",          "high",      "nb",
  "colon, colon-dose bands",              "colon",         "cola02w10", "colon",       FALSE,        "reference",          "high",      "nb",
  "liver, liver-dose bands",              "liver",         "liva02w10", "liva02w10",   FALSE,        "reference",          "high",      "nb")

stamp("B. data perturbations: ", nrow(gridB), " fits")
resB <- mclapply(seq_len(nrow(gridB)), function(k) {
  g <- gridB[k, ]
  ds <- build_ds(d, g$outcome, g$dose_var, band_var = g$band_var,
                 include_nic = g$include_nic,
                 ref_bands = strsplit(g$ref_bands, ",")[[1]], cap_band = g$cap)
  if (g$cap == "high_plus")
    ds <- ds %>% mutate(arm = ifelse(arm == "excluded", "high", arm))
  sub <- contrast_ds(ds, "high")
  xcv <- c("agex", "age", "sex", "city")
  if (g$include_nic) xcv <- c(xcv, "nic")     # Brenner's NIC-indicator adjustment
  run_spec(sub, g$spec, x_control_vars = xcv, count_model = g$model)
}, mc.cores = n_workers(6), mc.preschedule = FALSE)
B <- bind_rows(resB) %>% mutate(family = "B. data perturbation", .before = 1)
write_csv(B, file.path(TAB, "T21_sensitivity_data.csv"))
print(as.data.frame(B %>% dplyr::select(spec, n_cells, cases, att_serr, att_serr_lo,
  att_serr_hi, d_serr_female_male, d_serr_female_male_lo, d_serr_female_male_hi,
  d_serr_young_old, d_serr_young_old_lo, d_serr_young_old_hi)), digits = 3)

# =====================================================================
# C. Low-dose sub-cohorts (the range Dropkin 2016 re-analysed with GAMs)
# =====================================================================
# Dropkin restricts to women and to 0-100 / 0-20 mGy sub-cohorts and reports
# a raised RR at low dose.  The same contrasts are run here so that the
# CountBCF low-dose estimates can be placed next to his.
stamp("C. low-dose sub-cohorts")

gridC <- tribble(
  ~spec,                                    ~sexes,  ~dcatn_max,
  "5-<100 mGy, women only",                 2,       7L,
  "5-<100 mGy, men only",                   1,       7L,
  "5-<20 mGy, both sexes",                  0,       3L,
  "5-<20 mGy, women only",                  2,       3L)

resC <- mclapply(seq_len(nrow(gridC)), function(k) {
  g <- gridC[k, ]
  sub <- ds_main %>%
    filter(arm == "reference" | (arm == "low" & dcatn <= g$dcatn_max)) %>%
    { if (g$sexes == 0) . else filter(., sex == g$sexes) } %>%
    mutate(z = as.integer(arm == "low")) %>% droplevels()
  xv <- c("agex", "age", "sex", "city")
  if (g$sexes != 0) xv <- setdiff(xv, "sex")
  run_spec(sub, g$spec, x_control_vars = xv, x_moderate_vars = xv,
           ps_spec = if (g$sexes == 0) "main" else "simple")
}, mc.cores = n_workers(4), mc.preschedule = FALSE)
C <- bind_rows(resC) %>% mutate(family = "C. low-dose sub-cohort", .before = 1)
write_csv(C, file.path(TAB, "T29_lowdose_subcohorts.csv"))
print(as.data.frame(C %>% dplyr::select(spec, n_cells, cases, att_err, att_err_lo,
  att_err_hi, att_ear, att_ear_lo, att_ear_hi)), digits = 3)

write_csv(bind_rows(A, B, C), file.path(TAB, "T22_sensitivity_all.csv"))
stamp("03 done -> ", TAB)
