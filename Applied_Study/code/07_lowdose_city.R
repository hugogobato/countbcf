#!/usr/bin/env Rscript
# =====================================================================
# 07  Is the Nagasaki-vs-Hiroshima ERR difference at 5-<100 mGy real?
#
# §6.4 of RESULTS.md found the only subgroup that separates from the ATE
# in the low-dose contrast: Nagasaki ERR +0.21 against Hiroshima -0.03.
# City is a modifier Preston's ERR model cannot express, so if it is real
# it is genuinely new -- but city is also the one covariate with a real
# arm imbalance (58% Hiroshima in the reference arm, 78% in the low band).
#
# The same battery that §8-§10 applied to the high-dose contrast is
# applied here: the moderating-forest prior scale, the propensity model,
# trimming, the NIC group, and dropping each city's own dose banding.
#
# Every row is scored on BOTH footings, but the STANDARDISED ERR (ratio of
# person-year weighted counterfactual rates) is the one the manuscript
# quotes: it is the quantity comparable with a Poisson rate ratio and with
# Preston's case-weighted summaries.  The person-year weighted mean of the
# per-cell ERR (average of ratios) runs about 2.4x larger here and is kept
# only for cross-referencing against T15/T17.
#
# Two labelling points, following I. Shuryak's August 2026 review:
#   * ntree_moderate and the (base, power) depth prior are MODERATING-forest
#     knobs; a0 is the leaf concentration SHARED by both forests
#     (R/countbcf.R:90), so an a0 sweep moves the prognostic surface too and
#     is labelled accordingly;
#   * the propensity score for this contrast spans roughly [0.50, 0.73], so
#     a fixed [0.10, 0.90] trim removes nothing and is a no-op rather than a
#     passed test.  The trimming row therefore trims to the central 90% of
#     the observed score distribution, which does bite.
# =====================================================================

source(file.path(Sys.getenv("CBCF_ROOT", "."), "Applied_Study/code/common.R"))
suppressMessages(library(countbcf))

d <- read_lss()
ds <- build_ds(d, "solid", "cola02w10")
sub_low <- contrast_ds(ds, "low")

score_city <- function(r, spec, family = "") {
  cl <- add_bands(r$cells); w <- cl$pyr; tr <- which(cl$z == 1)
  err <- exp(r$logrr) - 1
  gi <- function(f) intersect(which(f), tr)
  hi <- draw_wmean(err, w, gi(cl$city_f == "Hiroshima"))
  na <- draw_wmean(err, w, gi(cl$city_f == "Nagasaki"))
  shi <- draw_serr(r$logmu, r$logrr, w, gi(cl$city_f == "Hiroshima"))
  sna <- draw_serr(r$logmu, r$logrr, w, gi(cl$city_f == "Nagasaki"))
  as_tibble(c(list(spec = spec, family = family, n_cells = nrow(cl),
                   cases = sum(cl$outcome),
                   pihat_min = min(cl$pihat), pihat_max = max(cl$pihat)),
              qsum(draw_wmean(err, w, tr), "att_err"),
              qsum(draw_serr(r$logmu, r$logrr, w, tr), "att_serr"),
              qsum(shi, "serr_hiroshima"), qsum(sna, "serr_nagasaki"),
              qsum(sna - shi, "d_serr_naga_hiro"),
              qsum(hi, "err_hiroshima"), qsum(na, "err_nagasaki"),
              qsum(na - hi, "d_naga_hiro")))
}

grid <- bind_rows(
  tibble(spec = paste0("ntree_moderate=", c(10, 25, 50, 100, 200)),
         family = "moderating forest", ntree_moderate = c(10, 25, 50, 100, 200),
         base_moderate = 0.25, power_moderate = 3, a0mult = 1,
         ps = "main", trim_q = NA_real_, nic = FALSE),
  tibble(spec = c("tau tree prior (0.05,4) tight", "tau tree prior (0.5,2) loose",
                  "tau tree prior (0.95,2) mu-like"),
         family = "moderating forest", ntree_moderate = 50,
         base_moderate = c(0.05, 0.5, 0.95), power_moderate = c(4, 2, 2),
         a0mult = 1, ps = "main", trim_q = NA_real_, nic = FALSE),
  tibble(spec = paste0("a0 x ", c(0.5, 2), " (both forests)"),
         family = "shared leaf concentration", ntree_moderate = 50,
         base_moderate = 0.25, power_moderate = 3, a0mult = c(0.5, 2),
         ps = "main", trim_q = NA_real_, nic = FALSE),
  tibble(spec = paste0("PS = ", c("const", "simple", "rich", "forest")),
         family = "propensity score", ntree_moderate = 50,
         base_moderate = 0.25, power_moderate = 3, a0mult = 1,
         ps = c("const", "simple", "rich", "forest"), trim_q = NA_real_, nic = FALSE),
  tibble(spec = "trim to central 90% of pihat", family = "data",
         ntree_moderate = 50, base_moderate = 0.25, power_moderate = 3,
         a0mult = 1, ps = "main", trim_q = 0.05, nic = FALSE),
  tibble(spec = "NIC included in control arm", family = "data",
         ntree_moderate = 50, base_moderate = 0.25, power_moderate = 3,
         a0mult = 1, ps = "main", trim_q = NA_real_, nic = TRUE))

a0_default <- 0.5 * (log(max(quantile(sub_low$outcome, 0.95), 1)) -
                     mean(log(sub_low$pyr)))
a0_default <- if (!is.finite(a0_default) || a0_default <= 0) 1 else a0_default

stamp("low-dose city check: ", nrow(grid), " fits")
res <- mclapply(seq_len(nrow(grid)), function(k) {
  g <- grid[k, ]
  s <- sub_low; xcv <- c("agex", "age", "sex", "city")
  if (g$nic) {
    dsn <- build_ds(d, "solid", "cola02w10", include_nic = TRUE)
    s <- contrast_ds(dsn, "low"); xcv <- c(xcv, "nic")
  }
  ph <- fit_ps(s, g$ps)
  if (!is.na(g$trim_q)) {                     # a trim that actually bites
    qq   <- quantile(ph, c(g$trim_q, 1 - g$trim_q))
    keep <- ph >= qq[1] & ph <= qq[2]
    s <- s[keep, , drop = FALSE]; ph <- ph[keep]
  }
  r <- tryCatch(fit_cbcf(s, label = g$spec, pihat = ph,
                         x_control_vars = xcv, ntree_moderate = g$ntree_moderate,
                         base_moderate = g$base_moderate,
                         power_moderate = g$power_moderate,
                         a0 = a0_default * g$a0mult,
                         nburn = NBURN_SENS, nsim = NSIM_SENS),
                error = function(e) { stamp("FAILED ", g$spec, ": ", conditionMessage(e)); NULL })
  if (is.null(r)) return(NULL)
  out <- score_city(r, g$spec, g$family); rm(r); gc(verbose = FALSE); out
}, mc.cores = n_workers(6), mc.preschedule = FALSE)

C <- bind_rows(res) %>%
  mutate(separates_std     = (d_serr_naga_hiro_lo > 0) | (d_serr_naga_hiro_hi < 0),
         separates_percell = (d_naga_hiro_lo > 0) | (d_naga_hiro_hi < 0))
write_csv(C, file.path(TAB, "T32_lowdose_city_robustness.csv"))
print(as.data.frame(C %>% dplyr::select(spec, family, n_cells, att_serr,
  serr_hiroshima, serr_nagasaki, d_serr_naga_hiro, d_serr_naga_hiro_lo,
  d_serr_naga_hiro_hi, separates_std, d_naga_hiro, separates_percell)), digits = 3)

# Same question asked of the parametric model, which can express a city x dose
# interaction even though Preston's ERR model does not include one.
ic <- d %>% filter(dcatn > 0, band %in% BANDS) %>%
  mutate(outcome = solid, sex_f = factor(sex, 1:2, c("Male", "Female")),
         city_f = factor(city, 1:2, c("Hiroshima", "Nagasaki")),
         band = droplevels(factor(as.character(band), levels = BANDS)))
m0 <- glm(outcome ~ band + sex_f * (ns(log(age), 4) + ns(agex, 3)) + city_f +
            offset(log(pyr)), poisson("log"), ic)
m1 <- update(m0, . ~ . + band:city_f)
lrt <- anova(m0, m1, test = "LRT")
cf <- coef(summary(m1))
gl <- tibble(term = rownames(cf)) %>% filter(grepl("band.*city", term)) %>%
  mutate(estimate = cf[term, 1], se = cf[term, 2],
         RR_ratio = exp(estimate), lo = exp(estimate - 1.96 * se),
         hi = exp(estimate + 1.96 * se), p = cf[term, 4])
gl$lrt_p <- lrt$`Pr(>Chi)`[2]
write_csv(gl, file.path(TAB, "T33_parametric_city_x_dose.csv"))
print(as.data.frame(gl), digits = 4)
stamp("07 done -> ", TAB)
