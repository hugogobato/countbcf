# =====================================================================
# CountBCF applied study v2 -- RERF Life Span Study solid-cancer
# incidence, FULL 1958-1998 person-year table (lssinc07.csv, 26,806 cells).
#
# Shared data construction, fitting helpers and posterior summaries.
# Sourced by 01_*.R .. 05_*.R.
#
# Data handling follows Preston et al. (2007, Radiat Res 168:1-64) and the
# instructions from A. Brenner (RERF):
#   * records with dcat == 1 (unknown DS02 dose) are dropped;
#   * dcat == 2 is split into NIC (distcat == 3, dose exactly 0) and
#     in-city survivors with 0-<5 mGy  -> dcatn coding, exactly as in the
#     Epicure log lss07sitemod_al.log;
#   * the causal contrasts use the four-level categorical dose variable
#     recommended for the applied-study section:
#         reference 0-<5 mGy | low 5-<100 | medium 100-<500 | high 500-<2000
#     with >=2000 mGy excluded;
#   * NIC survivors are NOT part of the causal contrast (they have no
#     positive probability of exposure); they enter a sensitivity fit as
#     extra control-arm cells with an NIC indicator in the prognostic
#     forest, the CountBCF analogue of the standard LSS baseline treatment
#     of the NIC group.
# =====================================================================

suppressMessages({
  library(dplyr); library(readr); library(tibble); library(tidyr)
  library(splines); library(MASS); library(parallel)
})

`%||%` <- function(a, b) if (is.null(a)) b else a

ROOT   <- Sys.getenv("CBCF_ROOT", unset = ".")
DATA   <- file.path(ROOT, "Applied_Study", "lssinc07.csv")
OUTDIR <- file.path(ROOT, "Applied_Study", "analysis")
TAB    <- file.path(OUTDIR, "tables")
PLT    <- file.path(OUTDIR, "plots")
for (p in c(OUTDIR, TAB, PLT)) dir.create(p, showWarnings = FALSE, recursive = TRUE)

stamp <- function(...) message(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ...)

# ---- MCMC budget ----------------------------------------------------
NBURN_MAIN <- 2000L; NSIM_MAIN <- 1000L
NBURN_SENS <- 1000L; NSIM_SENS <- 1000L
SEED       <- 20260724L

# Dose-band cutpoints on weighted (gamma + 10 x neutron) DS02 organ dose, Gy
CUT_REF <- 0.005; CUT_LOW <- 0.100; CUT_MED <- 0.500; CUT_HIGH <- 2.000
BANDS   <- c("reference", "low", "medium", "high")

# =====================================================================
# 1. Data
# =====================================================================

#' Read lssinc07.csv and attach the RERF-standard derived variables.
read_lss <- function(path = DATA) {
  d <- read_csv(path, show_col_types = FALSE)
  names(d) <- trimws(names(d))
  d %>%
    mutate(
      nic   = as.integer(distcat == 3),
      # Brenner's dcatn: 0 = unknown, 1 = NIC (dose 0), 2 = in-city 0-<5 mGy
      dcatn = case_when(dcat == 1 ~ 0L,
                        dcat == 2 & distcat == 3 ~ 1L,
                        dcat == 2 & distcat != 3 ~ 2L,
                        TRUE ~ as.integer(dcat)),
      # colon-dose four-level band (the recommended categorisation)
      band  = factor(case_when(dcatn == 0  ~ NA_character_,
                               dcatn == 1  ~ "NIC",
                               dcatn == 2  ~ "reference",
                               dcatn <= 7  ~ "low",
                               dcatn <= 14 ~ "medium",
                               dcatn <= 20 ~ "high",
                               TRUE        ~ "excluded"),
                     levels = c("NIC", BANDS, "excluded")),
      # survivor "cohort" followed through the table: the finest grouping
      # whose cells are repeated person-years of the SAME people
      cluster = factor(paste(city, sex, agxcat, dcatn, distcat, sep = "_")),
      # coarser grouping that CROSSES the dose bands, so a random intercept
      # on it does not compete with the treatment contrast
      cohort  = factor(paste(city, sex, agxcat, sep = "_")),
      # solid-cancer aggregates used in the site-composition sensitivity
      solid_no_bt   = pmax(solid - breast - thyroid, 0),
      solid_no_skin = pmax(solid - nmskin, 0),
      solid_ngs     = pmax(solid - breast - prost - femgenca - malgenca, 0),
      solid_no_stom = pmax(solid - stomach, 0))
}

#' Band a continuous organ dose on the same cutpoints (site-specific analyses).
organ_band <- function(dose, nic) {
  factor(case_when(nic == 1              ~ "NIC",
                   dose <  CUT_REF       ~ "reference",
                   dose <  CUT_LOW       ~ "low",
                   dose <  CUT_MED       ~ "medium",
                   dose <  CUT_HIGH      ~ "high",
                   TRUE                  ~ "excluded"),
         levels = c("NIC", BANDS, "excluded"))
}

#' Build the analysis frame for one endpoint.
#'
#' @param d         output of read_lss()
#' @param outcome   name of the case-count column
#' @param dose_var  organ-dose column used for the ERR-per-Gy conversion
#' @param band_var  "colon" (dcatn bands, default / primary) or an organ-dose
#'                  column name to re-band on
#' @param include_nic  keep the NIC cells (they get band == "NIC")
#' @param ref_bands    which bands make up the reference arm
build_ds <- function(d, outcome = "solid", dose_var = "cola02w10",
                     band_var = "colon", include_nic = FALSE,
                     ref_bands = "reference", cap_band = "high") {
  b <- if (identical(band_var, "colon")) d$band else organ_band(d[[band_var]], d$nic)
  keep_bands <- c(if (include_nic) "NIC", BANDS)
  if (identical(cap_band, "high_plus")) keep_bands <- c(keep_bands, "excluded")
  out <- d %>%
    mutate(band_use = b, outcome = .data[[outcome]], dose = .data[[dose_var]]) %>%
    filter(dcatn > 0, !is.na(band_use), band_use %in% keep_bands,
           !is.na(outcome), pyr > 0, dose >= 0) %>%
    mutate(band_use = droplevels(factor(as.character(band_use),
                                        levels = c("NIC", BANDS, "excluded"))),
           arm = ifelse(as.character(band_use) %in% ref_bands | band_use == "NIC",
                        "reference", as.character(band_use)),
           sex_f = factor(sex, levels = c(1, 2), labels = c("Male", "Female")),
           city_f = factor(city, levels = c(1, 2), labels = c("Hiroshima", "Nagasaki")),
           tsx = age - agex) %>%      # time since exposure
    dplyr::select(outcome, dose, pyr, upyr, subjects, band_use, arm, nic, cluster,
                  cohort, sex, city, sex_f, city_f, agex, age, year, tsx,
                  agxcat, agecat, time, dcatn, gdist)
  out
}

#' Subset to a single two-arm contrast: reference vs `treated`.
contrast_ds <- function(ds, treated) {
  ds %>%
    filter(arm %in% c("reference", treated)) %>%
    mutate(z = as.integer(arm == treated)) %>%
    droplevels()
}

# =====================================================================
# 2. Propensity score
# =====================================================================
# The exposure here is the survivor's location and shielding at 08:15 on
# 6/9 August 1945.  Distance from the hypocentre is deliberately EXCLUDED
# from the score: it is the (near-deterministic) dose mechanism, not a
# confounder, and conditioning on it destroys positivity.  What is left is
# the demographic composition of the arms, which is where the arms really
# do differ (Hiroshima:Nagasaki mix in particular).

PS_SPECS <- list(
  main    = z ~ city_f * sex_f + ns(agex, 4) + ns(age, 4) + city_f:ns(agex, 2),
  simple  = z ~ city_f + sex_f + agex + age,
  rich    = z ~ city_f * sex_f * (ns(agex, 5) + ns(age, 5))
)

# drop terms that involve a covariate which is constant in this subset
.drop_const <- function(f, sub, vars = c("city_f", "sex_f", "agex", "age")) {
  const <- vars[vapply(vars, function(v) length(unique(sub[[v]])) < 2, logical(1))]
  if (!length(const)) return(f)
  tl <- attr(terms(f), "term.labels")
  bad <- Reduce(`|`, lapply(const, function(v) grepl(v, tl, fixed = TRUE)))
  keep <- tl[!bad]
  reformulate(if (length(keep)) keep else "1", response = "z")
}

fit_ps <- function(sub, spec = "main") {
  if (spec == "const") return(rep(mean(sub$z), nrow(sub)))
  if (spec == "forest") {
    stopifnot(requireNamespace("ranger", quietly = TRUE))
    f <- .drop_const(z ~ city_f + sex_f + agex + age, sub)
    rf <- ranger::ranger(update(f, factor(z) ~ .), data = sub, probability = TRUE,
                         num.trees = 500, min.node.size = 25, seed = SEED,
                         num.threads = 1)
    return(as.numeric(rf$predictions[, "1"]))
  }
  m <- glm(.drop_const(PS_SPECS[[spec]], sub), binomial, data = sub)
  as.numeric(predict(m, type = "response"))
}

# =====================================================================
# 3. CountBCF wrapper
# =====================================================================

#' Fit CountBCF for one contrast and return compact posterior objects.
#'
#' Returns the per-cell posterior draw matrices for the two quantities the
#' study needs -- log RR (= tau, the ERR channel) and the per-cell excess
#' rate (EAR channel) -- plus the cell frame, so that every downstream
#' summary (subgroup CATE, CATE-minus-ATE contrast, ...) can be computed on
#' the draw scale rather than from point estimates.
#'
#' @param scale  excess cases per `scale` person-years (1e4 = Preston's unit)
fit_cbcf <- function(sub, label = "", ps_spec = "main", pihat = NULL,
                     x_control_vars = c("agex", "age", "sex", "city"),
                     x_moderate_vars = c("agex", "age", "sex", "city"),
                     count_model = "nb", nburn = NBURN_SENS, nsim = NSIM_SENS,
                     ntree_control = 250, ntree_moderate = 50,
                     base_moderate = 0.25, power_moderate = 3,
                     a0 = NA, include_pihat = "control",
                     cluster_re = FALSE, re_var = "cluster", seed = SEED, scale = 1e4,
                     keep_draws = TRUE, trim = NULL) {

  set.seed(seed)
  if (is.null(pihat)) pihat <- fit_ps(sub, ps_spec)
  if (!is.null(trim)) {
    keep <- pihat > trim & pihat < (1 - trim)
    sub <- sub[keep, , drop = FALSE]; pihat <- pihat[keep]
  }

  keepv <- function(v) v[vapply(v, function(k) length(unique(sub[[k]])) > 1, logical(1))]
  xcv <- keepv(x_control_vars); xmv <- keepv(x_moderate_vars)
  xc <- as.matrix(sub[, xcv, drop = FALSE])
  xm <- as.matrix(sub[, xmv, drop = FALSE])

  re_design <- matrix(1); re_q <- matrix(1)
  if (cluster_re) {
    # WARNING: this switch is INERT as a group-level effect.  countbcf()
    # samples the random effects but they do not enter the count likelihood
    # (the fitted log-mean is offset + mu_f + z * tau_f only); see
    # R/countbcf.R:109-113.  Any fit made with cluster_re = TRUE is therefore
    # numerically a no-random-effect fit with an extra sampled variance
    # component, and must not be labelled a random-intercept model.  To model
    # group structure, put the group indicator in x_control_vars instead.
    # countbcf() also sorts randeff_design itself (R/countbcf.R:295), so the
    # rows are passed through in the original observation order.
    g  <- droplevels(factor(sub[[re_var]]))
    D  <- model.matrix(~ g - 1)
    re_design <- D
    re_q <- matrix(1, nrow = ncol(D), ncol = 1)
  }

  fit <- countbcf::countbcf(
    y = sub$outcome, z = sub$z, x_control = xc, x_moderate = xm,
    pihat = pihat, offset = sub$pyr, count_model = count_model,
    include_pihat = include_pihat, nburn = nburn, nsim = nsim,
    ntree_control = ntree_control, ntree_moderate = ntree_moderate,
    base_moderate = base_moderate, power_moderate = power_moderate, a0 = a0,
    randeff_design = re_design, randeff_variance_component_design = re_q,
    update_interval = 1000)

  logmu <- fit$mu_f_post                        # nsim x n : log control rate
  logrr <- fit$tau_f_post                       # nsim x n : log RR = tau(x)
  kappa <- as.numeric(fit$kappa)
  re_sd <- if (cluster_re) as.numeric(fit$random_effects_sd) else NA_real_
  rm(fit); gc(verbose = FALSE)

  sub$pihat <- pihat
  sub$label <- label
  out <- list(cells = sub, label = label, kappa = kappa, re_sd = re_sd,
              scale = scale, n = nrow(sub))
  if (keep_draws) { out$logmu <- logmu; out$logrr <- logrr }
  out$summ <- summarise_fit(sub, logmu, logrr, label, scale)
  out
}

#' Per-cell excess absolute rate draws (cases per `scale` person-years).
ear_mat <- function(logmu, logrr, scale = 1e4)
  (exp(logmu + logrr) - exp(logmu)) * scale

# =====================================================================
# 4. Posterior summaries (all computed on the DRAW scale)
# =====================================================================

qsum <- function(x, nm) {
  x <- x[is.finite(x)]
  v <- if (!length(x)) rep(NA_real_, 6) else
    c(mean(x), quantile(x, c(0.025, 0.5, 0.975)), mean(x > 0), mean(x < 0))
  setNames(as.list(as.numeric(v)),
           paste0(nm, c("", "_lo", "_med", "_hi", "_pgt0", "_plt0")))
}

#' PY-weighted mean over a set of cells, per draw.
draw_wmean <- function(M, w, idx = NULL) {
  if (!is.null(idx)) { M <- M[, idx, drop = FALSE]; w <- w[idx] }
  if (!length(w) || sum(w) == 0) return(rep(NA_real_, nrow(M)))
  as.numeric(M %*% (w / sum(w)))
}

#' Person-year weighted SD across cells, per draw: how much heterogeneity
#' the posterior actually carries, as opposed to how much is visible in
#' point estimates.
draw_wsd <- function(M, w, idx = NULL) {
  if (!is.null(idx)) { M <- M[, idx, drop = FALSE]; w <- w[idx] }
  if (!length(w)) return(rep(NA_real_, nrow(M)))
  ww <- w / sum(w)
  m  <- as.numeric(M %*% ww); m2 <- as.numeric((M^2) %*% ww)
  sqrt(pmax(m2 - m^2, 0))
}

#' Standardised ERR over a set of cells, per draw: the ratio of the
#' person-year weighted treated and control rates, minus one.  This is the
#' quantity a parametric g-computation returns, and the one that must be
#' compared with the GLM rate ratio.  It is NOT the same as the person-year
#' weighted mean of the per-cell ERRs, which is the average of a ratio.
draw_serr <- function(logmu, logrr, w, idx = NULL) {
  if (!is.null(idx)) { logmu <- logmu[, idx, drop = FALSE]
                       logrr <- logrr[, idx, drop = FALSE]; w <- w[idx] }
  if (!length(w)) return(rep(NA_real_, nrow(logmu)))
  r0 <- exp(logmu); r1 <- exp(logmu + logrr)
  as.numeric((r1 %*% w) / (r0 %*% w)) - 1
}

#' Headline ATE / ATT summary for one fit.
summarise_fit <- function(cells, logmu, logrr, label, scale = 1e4) {
  w  <- cells$pyr
  tr <- which(cells$z == 1)
  ear <- ear_mat(logmu, logrr, scale)
  ate_ear <- draw_wmean(ear, w)
  att_ear <- draw_wmean(ear, w, tr)
  ate_err <- draw_wmean(exp(logrr) - 1, w)
  att_err <- draw_wmean(exp(logrr) - 1, w, tr)
  ate_serr <- draw_serr(logmu, logrr, w)
  att_serr <- draw_serr(logmu, logrr, w, tr)
  base <- draw_wmean(exp(logmu) * scale, w, tr)
  # mean dose contrast in Gy, PY-weighted, treated minus reference
  dbar <- weighted.mean(cells$dose[tr], w[tr]) -
          weighted.mean(cells$dose[cells$z == 0], w[cells$z == 0])
  as_tibble(c(list(label = label, n_cells = nrow(cells),
                   n_treated = length(tr),
                   cases = sum(cells$outcome), pyr = sum(cells$pyr),
                   treated_pyr = sum(w[tr]),
                   treated_cases = sum(cells$outcome[tr]),
                   dose_contrast_Gy = dbar),
              qsum(ate_ear, "ate_ear"), qsum(att_ear, "att_ear"),
              qsum(ate_err, "ate_err"), qsum(att_err, "att_err"),
              qsum(ate_serr, "ate_serr"), qsum(att_serr, "att_serr"),
              qsum(base, "base_rate"),
              qsum(draw_wsd(exp(logrr) - 1, w, tr), "het_sd_err"),
              qsum(draw_wsd(ear, w, tr), "het_sd_ear"),
              qsum(att_serr / dbar, "err_per_Gy")))
}

#' Subgroup CATEs with the (CATE - ATE) contrast Igor asked for.
#'
#' For every subgroup g we report the PY-weighted posterior CATE and, on the
#' same draws, CATE_g - ATE.  A subgroup is flagged `separates = TRUE` only
#' if the 95% credible interval of that DIFFERENCE excludes zero, i.e. the
#' subgroup is distinguishable from the global effect rather than merely
#' having a different point estimate.  The reference is the textbook
#' population ATE: both the subgroup CATE and the ATE are person-year weighted
#' over the whole cohort (all cells), so the ATE here is identical to the
#' headline `ate_*` in `summarise_fit`.  Set `treated_only = TRUE` only if you
#' deliberately want the treated-arm (ATT) framing instead.
subgroup_table <- function(cells, logmu, logrr, by, label = "", treated_only = FALSE,
                           min_cases = 10, scale = 1e4) {
  w   <- cells$pyr
  ear <- ear_mat(logmu, logrr, scale)
  idx_all <- if (treated_only) which(cells$z == 1) else seq_len(nrow(cells))
  ate_ear  <- draw_wmean(ear, w, idx_all)
  ate_err  <- draw_wmean(exp(logrr) - 1, w, idx_all)
  ate_serr <- draw_serr(logmu, logrr, w, idx_all)
  g <- do.call(interaction, c(as.list(cells[, by, drop = FALSE]),
                              list(drop = TRUE, sep = " | ")))
  res <- lapply(levels(g), function(lv) {
    idx <- intersect(which(g == lv), idx_all)
    if (!length(idx)) return(NULL)
    ncase <- sum(cells$outcome[idx])
    e  <- draw_wmean(ear, w, idx)
    r  <- draw_wmean(exp(logrr) - 1, w, idx)
    sr <- draw_serr(logmu, logrr, w, idx)
    de <- e - ate_ear; dr <- r - ate_err; dsr <- sr - ate_serr
    as_tibble(c(list(label = label, by = paste(by, collapse = " x "),
                     subgroup = lv, n_cells = length(idx),
                     cases = ncase, pyr = sum(w[idx]),
                     enough = ncase >= min_cases),
                qsum(e, "ear"), qsum(r, "err"), qsum(sr, "serr"),
                qsum(de, "d_ear"), qsum(dr, "d_err"), qsum(dsr, "d_serr")))
  })
  out <- bind_rows(res)
  if (!nrow(out)) return(out)
  out %>% mutate(
    separates_ear  = (d_ear_lo > 0) | (d_ear_hi < 0),
    separates_err  = (d_err_lo > 0) | (d_err_hi < 0),
    separates_serr = (d_serr_lo > 0) | (d_serr_hi < 0))
}

#' Split a numeric covariate into the bands Preston uses.
agex_band <- function(x) cut(x, c(-Inf, 10, 20, 30, 40, Inf),
                             labels = c("0-9", "10-19", "20-29", "30-39", "40+"),
                             right = FALSE)
age_band  <- function(x) cut(x, c(-Inf, 40, 60, 80, Inf),
                             labels = c("<40", "40-59", "60-79", "80+"),
                             right = FALSE)

add_bands <- function(cells) {
  cells %>% mutate(agex_b = agex_band(agex), age_b = age_band(age),
                   tsx_b = cut(tsx, c(-Inf, 20, 30, 40, Inf),
                               labels = c("<20", "20-29", "30-39", "40+"),
                               right = FALSE))
}

# =====================================================================
# 5. Misc
# =====================================================================

#' Gelman-Rubin R-hat and effective sample size from a list of chains.
rhat_ess <- function(chains) {
  m <- length(chains); n <- length(chains[[1]])
  X <- do.call(cbind, chains)
  W <- mean(apply(X, 2, var)); B <- n * var(colMeans(X))
  vhat <- (n - 1) / n * W + B / n
  rh <- sqrt(vhat / W)
  ess <- tryCatch(sum(coda::effectiveSize(coda::as.mcmc.list(
    lapply(chains, coda::as.mcmc)))), error = function(e) NA_real_)
  c(rhat = rh, ess = ess)
}

#' VanderWeele-Ding E-value for a rate ratio and its CI.
evalue_pt <- function(rr) { if (rr < 1) rr <- 1 / rr; rr + sqrt(rr * (rr - 1)) }
evalue_ci <- function(lo, hi) {
  if (lo <= 1 && hi >= 1) return(1)
  b <- if (hi < 1) hi else lo
  if (b < 1) b <- 1 / b
  b + sqrt(b * (b - 1))
}

n_workers <- function(want) max(1L, min(as.integer(want),
                                        as.integer(Sys.getenv("CBCF_CORES", "6"))))
