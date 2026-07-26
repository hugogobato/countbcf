#!/usr/bin/env Rscript
# =====================================================================
# 02  Main CountBCF fits.
#
#  * 2 endpoints (all solid, stomach) x 3 dose contrasts (low / medium /
#    high vs the 0-<5 mGy in-city reference) x 2 chains.
#  * ATE / ATT on the excess-rate (EAR) and excess-relative-risk (ERR)
#    scales, with MCMC convergence diagnostics.
#  * subgroup CATEs reported WITH the (CATE - ATE) contrast, so a subgroup
#    counts as a finding only if its posterior credibly separates from the
#    global effect.
#  * the CountBCF ERR surface is put side by side with the refitted
#    Preston (2007) ERR model to isolate what is genuinely new.
# =====================================================================

source(file.path(Sys.getenv("CBCF_ROOT", "."), "Applied_Study_v2/code/common.R"))
suppressMessages(library(countbcf))

d <- read_lss()

SPECS <- tribble(
  ~endpoint, ~outcome,  ~dose_var,    ~nice,
  "solid",   "solid",   "cola02w10",  "All solid cancers (colon dose)",
  "stomach", "stomach", "stoa02w10",  "Stomach cancer (stomach dose)")
CONTRASTS <- c("low", "medium", "high")
CHAINS    <- c(1L, 2L)

jobs <- expand_grid(i = seq_len(nrow(SPECS)), treated = CONTRASTS, chain = CHAINS) %>%
  mutate(endpoint = SPECS$endpoint[i], outcome = SPECS$outcome[i],
         dose_var = SPECS$dose_var[i],
         label = paste(endpoint, treated, sep = ":"),
         job = paste(label, chain, sep = "/"))

stamp("main fits: ", nrow(jobs), " chains on ", n_workers(5), " workers")

run_job <- function(k) {
  j <- jobs[k, ]
  ds  <- build_ds(d, j$outcome, j$dose_var, band_var = "colon", include_nic = FALSE)
  sub <- contrast_ds(ds, j$treated)
  r <- fit_cbcf(sub, label = j$label, ps_spec = "main",
                nburn = NBURN_MAIN, nsim = NSIM_MAIN,
                seed = SEED + 1000L * j$chain, keep_draws = TRUE)
  r$chain <- j$chain; r$endpoint <- j$endpoint; r$treated <- j$treated
  r
}

t0 <- Sys.time()
fits <- mclapply(seq_len(nrow(jobs)), run_job,
                 mc.cores = n_workers(5), mc.preschedule = FALSE)
stamp("fits done in ", round(difftime(Sys.time(), t0, units = "mins"), 1), " min")
bad <- vapply(fits, function(f) inherits(f, "try-error") || is.null(f$cells), logical(1))
if (any(bad)) { print(fits[bad]); stop("some chains failed") }
names(fits) <- jobs$job

# =====================================================================
# 1. Convergence diagnostics, then pool chains
# =====================================================================
conv <- bind_rows(lapply(unique(jobs$label), function(lab) {
  ff <- fits[jobs$job[jobs$label == lab]]
  ch <- lapply(ff, function(f) draw_wmean(ear_mat(f$logmu, f$logrr), f$cells$pyr,
                                          which(f$cells$z == 1)))
  ch_err <- lapply(ff, function(f) draw_wmean(exp(f$logrr) - 1, f$cells$pyr,
                                              which(f$cells$z == 1)))
  ka <- lapply(ff, function(f) f$kappa)
  tibble(label = lab,
         rhat_att_ear = rhat_ess(ch)["rhat"],   ess_att_ear = rhat_ess(ch)["ess"],
         rhat_att_err = rhat_ess(ch_err)["rhat"], ess_att_err = rhat_ess(ch_err)["ess"],
         rhat_kappa = rhat_ess(ka)["rhat"],
         kappa_mean = mean(unlist(ka)),
         kappa_lo = quantile(unlist(ka), .025), kappa_hi = quantile(unlist(ka), .975))
}))
write_csv(conv, file.path(TAB, "T11_mcmc_convergence.csv"))
print(as.data.frame(conv), digits = 4)

# pool the two chains for every label
pool <- lapply(unique(jobs$label), function(lab) {
  ff <- fits[jobs$job[jobs$label == lab]]
  list(label = lab, endpoint = ff[[1]]$endpoint, treated = ff[[1]]$treated,
       cells = ff[[1]]$cells,
       logmu = do.call(rbind, lapply(ff, `[[`, "logmu")),
       logrr = do.call(rbind, lapply(ff, `[[`, "logrr")),
       kappa = unlist(lapply(ff, `[[`, "kappa")))
})
names(pool) <- unique(jobs$label)
rm(fits); gc(verbose = FALSE)

# =====================================================================
# 2. ATE / ATT table + convergence to the parametric yardstick
# =====================================================================
ate_tab <- bind_rows(lapply(pool, function(p)
  summarise_fit(p$cells, p$logmu, p$logrr, p$label) %>%
    mutate(endpoint = p$endpoint, treated = p$treated,
           kappa = mean(p$kappa), .before = 1)))
write_csv(ate_tab, file.path(TAB, "T12_countbcf_ate.csv"))
print(as.data.frame(ate_tab %>% dplyr::select(label, n_cells, cases, dose_contrast_Gy,
  base_rate, att_ear, att_ear_lo, att_ear_hi, att_serr, att_serr_lo, att_serr_hi,
  att_err, att_err_lo, att_err_hi, err_per_Gy, err_per_Gy_lo, err_per_Gy_hi,
  kappa)), digits = 4)

param <- read_csv(file.path(TAB, "T07_parametric_ate.csv"), show_col_types = FALSE)
converge <- bind_rows(
  param %>% filter(model %in% c("NB (adjusted)", "Poisson (adjusted)")) %>%
    transmute(endpoint, treated, method = model, ear, ear_lo, ear_hi, rr, rr_lo, rr_hi),
  ate_tab %>% transmute(endpoint, treated, method = "CountBCF (NB)",
                        ear = ate_ear, ear_lo = ate_ear_lo, ear_hi = ate_ear_hi,
                        rr = 1 + ate_serr, rr_lo = 1 + ate_serr_lo,
                        rr_hi = 1 + ate_serr_hi)) %>%
  arrange(endpoint, treated, method)
write_csv(converge, file.path(TAB, "T13_ate_convergence.csv"))
print(as.data.frame(converge), digits = 4)

# E-values for the CountBCF rate ratios
ev <- ate_tab %>%
  transmute(label, endpoint, treated, rr = 1 + att_serr,
            rr_lo = 1 + att_serr_lo, rr_hi = 1 + att_serr_hi,
            evalue_point = vapply(rr, evalue_pt, numeric(1)),
            evalue_ci = mapply(evalue_ci, rr_lo, rr_hi))
write_csv(ev, file.path(TAB, "T14_evalues.csv"))

# =====================================================================
# 3. Subgroup CATEs and the (CATE - ATE) separation test
# =====================================================================
stamp("subgroup tables")
BY <- list(c("sex_f"), c("agex_b"), c("age_b"), c("city_f"), c("tsx_b"),
           c("sex_f", "agex_b"), c("sex_f", "age_b"), c("city_f", "sex_f"),
           c("city_f", "agex_b"))

sg <- bind_rows(lapply(pool, function(p) {
  cl <- add_bands(p$cells)
  bind_rows(lapply(BY, function(b)
    subgroup_table(cl, p$logmu, p$logrr, b, label = p$label, treated_only = TRUE)))
}))
write_csv(sg, file.path(TAB, "T15_subgroup_cate.csv"))

stamp("subgroups separating from the ATE (95% CrI on CATE - ATE excludes 0):")
print(as.data.frame(sg %>% filter(enough, separates_err | separates_ear) %>%
  dplyr::select(label, by, subgroup, cases, err, err_lo, err_hi,
                d_err, d_err_lo, d_err_hi, separates_err,
                ear, d_ear, d_ear_lo, d_ear_hi, separates_ear)), digits = 3)

# per-cell CATE dump (compressed)
cells_out <- bind_rows(lapply(pool, function(p) {
  cl <- add_bands(p$cells)
  cl$err      <- colMeans(exp(p$logrr) - 1)
  cl$err_lo   <- apply(exp(p$logrr) - 1, 2, quantile, .025)
  cl$err_hi   <- apply(exp(p$logrr) - 1, 2, quantile, .975)
  em          <- ear_mat(p$logmu, p$logrr)
  cl$ear      <- colMeans(em)
  cl$ear_lo   <- apply(em, 2, quantile, .025)
  cl$ear_hi   <- apply(em, 2, quantile, .975)
  cl$base_rate <- colMeans(exp(p$logmu)) * 1e4
  cl$label    <- p$label
  cl
}))
write_csv(cells_out, file.path(TAB, "T16_countbcf_cate_cells.csv.gz"))

# how much of the per-cell heterogeneity is resolved by the posterior?
cell_res <- cells_out %>% filter(z == 1) %>% group_by(label) %>%
  summarise(n_treated_cells = n(),
            pct_err_gt0 = mean(err > 0), pct_err_CrI_excl0 = mean(err_lo > 0),
            pct_ear_CrI_excl0 = mean(ear_lo > 0),
            err_med = median(err), err_q10 = quantile(err, .10),
            err_q90 = quantile(err, .90),
            ear_med = median(ear), ear_q10 = quantile(ear, .10),
            ear_q90 = quantile(ear, .90), .groups = "drop")
write_csv(cell_res, file.path(TAB, "T16b_cell_level_heterogeneity.csv"))
print(as.data.frame(cell_res), digits = 3)

# =====================================================================
# 4. What is genuinely new relative to Preston's parametric ERR model?
# =====================================================================
stamp("CountBCF vs the refitted Preston ERR model")
# The Preston ERR surface is evaluated from the fitted parameters directly
# (beta_s * d * exp(alpha_e (e - 30) + alpha_a log(a/70))), so that every
# person-year cell gets the value implied by ITS own dose, sex, age at
# exposure and attained age -- no key matching, no within-band collapsing.
pf <- readRDS(file.path(OUTDIR, "preston_fit.rds"))$preston_par
preston_surface <- function(cl, kind = "ERR") {
  par <- pf[[kind]]; k <- pf$p0
  bs <- ifelse(cl$sex == 2, par[k + 2], par[k + 1])
  bs * cl$dose * exp(par[k + 3] * (cl$agex - 30) + par[k + 4] * log(cl$age / 70))
}

compare_preston <- function(p, by) {
  cl <- add_bands(p$cells)
  cl$preston_err <- preston_surface(cl, "ERR")
  cl$preston_ear <- preston_surface(cl, "EAR")
  tr <- which(cl$z == 1 & is.finite(cl$preston_err))
  if (!length(tr)) return(NULL)
  w <- cl$pyr
  err_draws <- exp(p$logrr) - 1
  ear_draws <- ear_mat(p$logmu, p$logrr)
  r0        <- exp(p$logmu)
  g <- do.call(interaction, c(as.list(cl[, by, drop = FALSE]),
                              list(drop = TRUE, sep = " | ")))
  bind_rows(lapply(levels(g), function(lv) {
    idx <- intersect(which(g == lv), tr)
    if (length(idx) < 2) return(NULL)
    cb  <- draw_wmean(err_draws, w, idx)
    cbe <- draw_wmean(ear_draws, w, idx)
    pk  <- weighted.mean(cl$preston_err[idx], w[idx])
    pke <- weighted.mean(cl$preston_ear[idx], w[idx])
    # baseline-rate standardised versions: both surfaces weighted by the same
    # CountBCF control-arm rate, so they are directly comparable to a GLM
    # rate ratio (and to Preston's own case-weighted summaries).
    wr  <- sweep(r0[, idx, drop = FALSE], 2, w[idx], "*")
    den <- rowSums(wr)
    cb_s <- rowSums(wr * err_draws[, idx, drop = FALSE]) / den
    pk_s <- as.numeric(wr %*% cl$preston_err[idx]) / den
    as_tibble(c(list(label = p$label, by = paste(by, collapse = " x "), subgroup = lv,
                     n_cells = length(idx), cases = sum(cl$outcome[idx]),
                     preston_err = pk, preston_ear = pke),
                qsum(cb, "cbcf_err"), qsum(cb - pk, "diff"),
                qsum(cbe, "cbcf_ear"), qsum(cbe - pke, "diff_ear"),
                qsum(cb_s, "cbcf_serr"), qsum(pk_s, "preston_serr"),
                qsum(cb_s - pk_s, "diff_serr")))
  }))
}

newness <- bind_rows(lapply(pool, function(p)
  bind_rows(lapply(list("all", c("sex_f"), c("agex_b"), c("age_b"), c("city_f"),
                        c("tsx_b"), c("sex_f", "agex_b")), function(b) {
    if (identical(b, "all")) {
      p2 <- p; p2$cells$.all <- factor("overall")
      compare_preston(p2, ".all")
    } else compare_preston(p, b)
  })))) %>%
  mutate(departs = (diff_lo > 0) | (diff_hi < 0),
         departs_serr = (diff_serr_lo > 0) | (diff_serr_hi < 0),
         departs_ear  = (diff_ear_lo > 0) | (diff_ear_hi < 0))
write_csv(newness, file.path(TAB, "T17_vs_preston_model.csv"))
stamp("subgroups whose CountBCF ERR departs from the Preston-model ERR:")
print(as.data.frame(newness %>% filter(departs, cases >= 20) %>%
  dplyr::select(label, by, subgroup, cases, preston_err, cbcf_err,
                cbcf_err_lo, cbcf_err_hi, diff, diff_lo, diff_hi,
                preston_serr, cbcf_serr, diff_serr_lo, diff_serr_hi,
                departs_serr)), digits = 3)

# cell-level agreement between the two surfaces (treated cells)
agree <- bind_rows(lapply(pool, function(p) {
  cl <- p$cells
  pk <- preston_surface(cl, "ERR")
  tr <- which(cl$z == 1 & is.finite(pk))
  pke <- preston_surface(cl, "EAR")
  cb <- colMeans(exp(p$logrr) - 1)[tr]
  cbe <- colMeans(ear_mat(p$logmu, p$logrr))[tr]
  w  <- cl$pyr[tr]; pk <- pk[tr]; pke <- pke[tr]
  tibble(label = p$label, n = length(tr),
         cor_err_pearson = cor(cb, pk), cor_err_spearman = cor(cb, pk, method = "spearman"),
         wmean_err_cbcf = weighted.mean(cb, w), wmean_err_preston = weighted.mean(pk, w),
         wrmse_err = sqrt(weighted.mean((cb - pk)^2, w)),
         cor_ear_pearson = cor(cbe, pke), cor_ear_spearman = cor(cbe, pke, method = "spearman"),
         wmean_ear_cbcf = weighted.mean(cbe, w), wmean_ear_preston = weighted.mean(pke, w),
         wrmse_ear = sqrt(weighted.mean((cbe - pke)^2, w)))
}))
write_csv(agree, file.path(TAB, "T18_surface_agreement.csv"))
print(as.data.frame(agree), digits = 4)

# =====================================================================
# 5. Dose-response shape implied by the three band contrasts
# =====================================================================
shape <- ate_tab %>%
  transmute(endpoint, treated, dose_contrast_Gy,
            err = att_err, err_lo = att_err_lo, err_hi = att_err_hi,
            err_per_Gy, err_per_Gy_lo, err_per_Gy_hi,
            serr = att_serr, serr_lo = att_serr_lo, serr_hi = att_serr_hi,
            ear = att_ear, ear_lo = att_ear_lo, ear_hi = att_ear_hi)
write_csv(shape, file.path(TAB, "T19_dose_response_shape.csv"))
print(as.data.frame(shape), digits = 4)

saveRDS(lapply(pool, function(p) list(label = p$label, endpoint = p$endpoint,
                                      treated = p$treated, cells = p$cells,
                                      err_mean = colMeans(exp(p$logrr) - 1),
                                      ear_mean = colMeans(ear_mat(p$logmu, p$logrr)))),
        file.path(OUTDIR, "main_cells.rds"))
stamp("02 done -> ", TAB)
