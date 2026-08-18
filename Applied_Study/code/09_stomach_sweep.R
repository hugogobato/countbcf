#!/usr/bin/env Rscript
# =====================================================================
# 09  The regularisation sweep, repeated for stomach cancer.
#
#  Section 1.6.4 reports that CountBCF returns an essentially flat ERR
#  surface for stomach at high dose (posterior PY-weighted SD of the
#  per-cell ERR = 0.004) and argues that this is differential shrinkage
#  behaving correctly on a weak signal rather than over-shrinkage.  Until
#  now that argument rested on the published significance tests plus a
#  single CountBCF fit: the twelve-setting regularisation sweep of
#  Section 1.8 had only ever been run for the all-solid endpoint, and
#  Section 1.15 recorded the gap as an open item.
#
#  This script closes it.  The grid is byte-for-byte the grid of
#  03_sensitivity.R part A, so the stomach column is directly comparable
#  with the all-solid one:
#     * ntree_moderate in {10, 25, 50, 100, 200}   (moderating forest only)
#     * (base, power) in {(0.05,4), (0.25,3), (0.5,2), (0.95,2)}
#                                                  (moderating forest only)
#     * a0 x {0.5, 1, 2}                           (SHARED by both forests,
#       R/countbcf.R:90, implied treatment-effect prior scale
#       a0 * sqrt(ntree_moderate / ntree_control))
#
#  The question a sweep has to answer for a NULL heterogeneity result is
#  the mirror image of the one it answers for a positive one.  For
#  all-solid we ask whether a modifier that separates at the defaults
#  still separates when the moderating forest is tightened.  For stomach
#  we ask whether a modifier that does NOT separate at the defaults
#  appears once the moderating forest is loosened.  So the reported
#  quantities are the flatness of the surface (het_sd_err), the number of
#  the 49 subgroup separation tests that fire on the standardised scale,
#  the largest subgroup deviation from the ATE in either direction, and
#  the three contrasts the all-solid sweep features.
# =====================================================================

source(file.path(Sys.getenv("CBCF_ROOT", "."), "Applied_Study/code/common.R"))
suppressMessages(library(countbcf))

d <- read_lss()

# The 49 separation tests per fit are the same nine groupings used in
# 02_countbcf_main.R, so the counts below sit on the same denominator as
# T15b_multiplicity_accounting.csv.
BY <- list(c("sex_f"), c("agex_b"), c("age_b"), c("city_f"), c("tsx_b"),
           c("sex_f", "agex_b"), c("sex_f", "age_b"), c("city_f", "sex_f"),
           c("city_f", "agex_b"))

score_fit <- function(r, spec) {
  cl <- add_bands(r$cells); w <- cl$pyr; tr <- which(cl$z == 1)
  err <- exp(r$logrr) - 1
  gi <- function(f) intersect(which(f), tr)
  sg <- function(k) draw_serr(r$logmu, r$logrr, w, gi(k))
  sm  <- sg(cl$sex_f == "Male");       sf  <- sg(cl$sex_f == "Female")
  sy  <- sg(cl$agex < 20);             so  <- sg(cl$agex >= 40)
  shi <- sg(cl$city_f == "Hiroshima"); sna <- sg(cl$city_f == "Nagasaki")

  # the separation tests, on the same footing and with the same 10-case
  # minimum as the main analysis
  st <- bind_rows(lapply(BY, function(b)
    subgroup_table(cl, r$logmu, r$logrr, b, label = spec, treated_only = FALSE)))
  st_ok <- st %>% filter(enough)

  as_tibble(c(list(spec = spec, n_cells = nrow(cl), cases = sum(cl$outcome),
                   treated_cases = sum(cl$outcome[tr]),
                   kappa = mean(r$kappa),
                   n_tests = nrow(st_ok),
                   n_sep_serr = sum(st_ok$separates_serr),
                   n_sep_percell = sum(st_ok$separates_err),
                   n_sep_ear = sum(st_ok$separates_ear),
                   max_abs_d_serr = max(abs(st_ok$d_serr)),
                   widest_d_serr_ci = max(st_ok$d_serr_hi - st_ok$d_serr_lo)),
              qsum(draw_serr(r$logmu, r$logrr, w, tr), "att_serr"),
              qsum(draw_wmean(err, w, tr), "att_err"),
              qsum(draw_wsd(err, w, tr), "het_sd_err"),
              qsum(sf - sm,  "d_serr_female_male"),
              qsum(sy - so,  "d_serr_young_old"),
              qsum(sna - shi, "d_serr_naga_hiro")))
}

run_spec <- function(sub, spec, ...) {
  r <- tryCatch(fit_cbcf(sub, label = spec, nburn = NBURN_SENS, nsim = NSIM_SENS, ...),
                error = function(e) { stamp("FAILED ", spec, ": ", conditionMessage(e)); NULL })
  if (is.null(r)) return(NULL)
  out <- score_fit(r, spec)
  rm(r); gc(verbose = FALSE)
  out
}

# The featured stomach contrast, identical to the main fit in 02: stomach
# cases, stomach dose for the per-Gy conversion, colon-dose banding.
ds_stom  <- build_ds(d, "stomach", "stoa02w10", band_var = "colon", include_nic = FALSE)
sub_stom <- contrast_ds(ds_stom, "high")
stamp("stomach high-dose contrast: ", nrow(sub_stom), " cells, ",
      sum(sub_stom$outcome), " cases, ",
      sum(sub_stom$outcome[sub_stom$z == 1]), " exposed cases")

grid <- bind_rows(
  tibble(spec = paste0("ntree_moderate=", c(10, 25, 50, 100, 200)),
         ntree_moderate = c(10, 25, 50, 100, 200),
         base_moderate = 0.25, power_moderate = 3, a0mult = 1),
  tibble(spec = c("tau tree prior (0.05,4) tight", "tau tree prior (0.25,3) default",
                  "tau tree prior (0.5,2) loose", "tau tree prior (0.95,2) mu-like"),
         ntree_moderate = 50,
         base_moderate = c(0.05, 0.25, 0.5, 0.95),
         power_moderate = c(4, 3, 2, 2), a0mult = 1),
  tibble(spec = paste0("a0 x ", c(0.5, 1, 2), " (both forests)"),
         ntree_moderate = 50,
         base_moderate = 0.25, power_moderate = 3, a0mult = c(0.5, 1, 2)))

# default a0 computed on this contrast's own outcome, exactly as in 03
a0_default <- 0.5 * (log(max(quantile(sub_stom$outcome, 0.95), 1)) -
                     mean(log(sub_stom$pyr)))
a0_default <- if (!is.finite(a0_default) || a0_default <= 0) 1 else a0_default
stamp("default a0 for the stomach contrast = ", signif(a0_default, 4))

stamp("stomach regularisation sweep: ", nrow(grid), " fits")
res <- mclapply(seq_len(nrow(grid)), function(k) {
  g <- grid[k, ]
  run_spec(sub_stom, g$spec, ntree_moderate = g$ntree_moderate,
           base_moderate = g$base_moderate, power_moderate = g$power_moderate,
           a0 = a0_default * g$a0mult)
}, mc.cores = n_workers(6), mc.preschedule = FALSE)

S <- bind_rows(res) %>% mutate(family = "stomach tau regularisation", .before = 1)
# as in 03: the (0.25,3) depth prior and a0 x 1 ARE the default specification,
# so three of the twelve rows are the same fit reported once per sweep
S$is_default_repeat <- as.vector(duplicated(round(cbind(
  S$att_serr, S$het_sd_err, S$d_serr_young_old), 8)))
write_csv(S, file.path(TAB, "T41_stomach_regularisation.csv"))

print(as.data.frame(S %>% dplyr::select(spec, att_serr, att_serr_lo, att_serr_hi,
  het_sd_err, het_sd_err_lo, het_sd_err_hi, n_sep_serr, n_sep_percell,
  max_abs_d_serr, widest_d_serr_ci)), digits = 3)
cat("\n-- subgroup contrasts on the standardised footing --\n")
print(as.data.frame(S %>% dplyr::select(spec, d_serr_female_male,
  d_serr_female_male_lo, d_serr_female_male_hi, d_serr_young_old,
  d_serr_young_old_lo, d_serr_young_old_hi, d_serr_naga_hiro,
  d_serr_naga_hiro_lo, d_serr_naga_hiro_hi)), digits = 3)
cat("\nrows that repeat the default specification: ", sum(S$is_default_repeat), "\n")
cat("specifications in which ANY subgroup separates on the standardised scale: ",
    sum(S$n_sep_serr > 0), " of ", nrow(S), "\n")
cat("range of het_sd_err across the sweep: ",
    signif(min(S$het_sd_err), 3), " to ", signif(max(S$het_sd_err), 3), "\n")
cat("range of att_serr across the sweep:   ",
    signif(min(S$att_serr), 3), " to ", signif(max(S$att_serr), 3), "\n")

stamp("09 done -> ", TAB)
