#!/usr/bin/env Rscript
# Exact posterior evaluation at exposure age 30 and attained age 70.
# This remains a binary 0.5-<2 Gy versus <5 mGy contrast. Dividing its
# ERR by the sex-specific observed mean dose contrast is a descriptive
# rescaling, not an estimate of a derivative or a new continuous-dose fit.
# Only the requested prediction rows are retained; full trees are not saved.
# Use one fresh R process per chain with a 6 GiB address-space limit and
# OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1.

source(file.path(Sys.getenv("CBCF_ROOT", "."), "Applied_Study/code/common.R"))
suppressMessages(library(countbcf))
stopifnot(requireNamespace("posterior", quietly = TRUE))
PROFILE_NBURN <- as.integer(Sys.getenv("CBCF_PROFILE_NBURN", "5000"))
PROFILE_NSIM <- as.integer(Sys.getenv("CBCF_PROFILE_NSIM", "1000"))
PROFILE_THIN <- as.integer(Sys.getenv("CBCF_PROFILE_THIN", "10"))
stopifnot(PROFILE_NBURN >= 0L, PROFILE_NSIM > 0L, PROFILE_THIN > 0L,
          PROFILE_NBURN %% PROFILE_THIN == 0L)


sub <- contrast_ds(build_ds(read_lss(), "solid", "cola02w10"), "high")
ps_model <- glm(.drop_const(PS_SPECS$main, sub), binomial, data = sub)
pihat <- as.numeric(predict(ps_model, type = "response"))
vars <- c("agex", "age", "sex", "city")
xc <- xm <- as.matrix(sub[, vars])
grid <- expand_grid(sex = 1:2, city = 1:2) %>%
  mutate(agex = 30, age = 70,
         sex_f = factor(sex, 1:2, c("Male", "Female")),
         city_f = factor(city, 1:2, c("Hiroshima", "Nagasaki")))
grid$pihat <- as.numeric(predict(ps_model, newdata = grid, type = "response"))
new_xm <- as.matrix(grid[, vars])
check <- unique(round(seq(1, nrow(sub), length.out = 4)))
predict_x <- rbind(new_xm, xc[check, , drop = FALSE])
predict_ps <- c(grid$pihat, pihat[check])
grid$city_weight <- vapply(seq_len(nrow(grid)), function(i) {
  z1 <- sub$z == 1 & sub$sex == grid$sex[i]
  sum(sub$pyr[z1 & sub$city == grid$city[i]]) / sum(sub$pyr[z1])
}, numeric(1))
grid$dose_contrast_Gy <- vapply(grid$sex, function(s) {
  j1 <- sub$sex == s & sub$z == 1
  j0 <- sub$sex == s & sub$z == 0
  weighted.mean(sub$dose[j1], sub$pyr[j1]) - weighted.mean(sub$dose[j0], sub$pyr[j0])
}, numeric(1))

# Checkpoint identity covers the data, propensity construction, sampler R/C++
# sources, installed native binary, R version and all fitting settings. This
# prevents silently reusing draws from a different source or input snapshot.
identity_files <- c("Applied_Study/lssinc07.csv", "Applied_Study/code/common.R",
                    "R/countbcf.R", "R/get_forest_fit.R", "Applied_Study/code/10_profile_comparison.R",
                    list.files(file.path(ROOT, "src"), pattern = "\\.(cpp|h)$", full.names = FALSE) %>%
                      file.path("src", .))
identity_md5 <- tools::md5sum(file.path(ROOT, identity_files))
names(identity_md5) <- identity_files
native_file <- file.path(find.package("countbcf"), "libs", paste0("countbcf", .Platform$dynlib.ext))
profile_identity <- list(source_and_data_md5 = identity_md5,
  installed_native_md5 = unname(tools::md5sum(native_file)),
  R_version = R.version.string, countbcf_version = as.character(packageVersion("countbcf")),
  settings = list(nburn = PROFILE_NBURN, nsim = PROFILE_NSIM, nthin = PROFILE_THIN,
    seeds = SEED + 1000L * (1:2), count_model = "nb", include_pihat = "control",
    ntree_control = 250, ntree_moderate = 50, base_moderate = .25,
    power_moderate = 3, a0 = NA_real_, variables = vars,
    dose_group = "0.5-<2 Gy versus in-city <5 mGy", grid = grid))
saveRDS(profile_identity, file.path(OUTDIR, "profile_age30_70_identity.rds"))
if (identical(Sys.getenv("CBCF_PROFILE_IDENTITY_ONLY"), "1")) quit(save = "no", status = 0)

profile_memory <- function(stage, require_headroom = FALSE) {
  if (!file.exists("/proc/meminfo")) return(invisible(NULL))
  available <- as.numeric(strsplit(trimws(sub("^MemAvailable:\\s*", "",
    grep("^MemAvailable:", readLines("/proc/meminfo"), value = TRUE))), "\\s+")[[1]][1])
  rss <- if (file.exists("/proc/self/status"))
    grep("^Vm(RSS|Peak):", readLines("/proc/self/status"), value = TRUE) else character()
  stamp(stage, ": host available ", round(available / 1024^2, 1), " GiB; ", paste(rss, collapse = "; "))
  # A capped 6 GiB process must leave at least 8 GiB for the desktop/other work.
  if (require_headroom && is.finite(available) && available < 14 * 1024^2)
    stop("Insufficient available memory for a profile fit while reserving 8 GiB for other work")
  invisible(NULL)
}

profile_chain <- function(ch) {
  checkpoint <- file.path(OUTDIR, paste0("profile_age30_70_chain", ch, ".rds"))
  if (file.exists(checkpoint) && !identical(Sys.getenv("CBCF_PROFILE_FORCE"), "1")) {
    old <- readRDS(checkpoint)
    stopifnot(identical(old$seed, SEED + 1000L * ch),
              identical(old$nburn, PROFILE_NBURN), identical(old$nsim, PROFILE_NSIM),
              isTRUE(all.equal(old$grid, grid)),
              identical(old$source_identity, profile_identity))
    stamp("using completed profile checkpoint for chain ", ch)
    return(old)
  }
  profile_memory(paste0("before profile chain ", ch), require_headroom = TRUE)
  set.seed(SEED + 1000L * ch)
  fit <- countbcf::countbcf(y = sub$outcome, z = sub$z,
    x_control = xc, x_moderate = xm, pihat = pihat, offset = sub$pyr,
    count_model = "nb", include_pihat = "control",
    nburn = PROFILE_NBURN, nsim = PROFILE_NSIM, nthin = PROFILE_THIN,
    ntree_control = 250, ntree_moderate = 50,
    base_moderate = .25, power_moderate = 3,
    a0 = NA_real_, return_trees = FALSE, update_interval = 1000,
    x_control_predict = predict_x, x_moderate_predict = predict_x,
    pihat_predict = predict_ps)
  profile_memory(paste0("sampler returned for chain ", ch))
  full_mu <- fit$mu_f_post; full_tau <- fit$tau_f_post
  kappa <- fit$kappa
  logmu <- fit$mu_f_predict[, seq_len(nrow(grid)), drop = FALSE]
  logrr <- fit$tau_f_predict[, seq_len(nrow(grid)), drop = FALSE]
  # New-row evaluation must exactly reproduce saved training coefficients.
  stopifnot(isTRUE(all.equal(unname(fit$mu_f_predict[, nrow(grid) + seq_along(check), drop = FALSE]),
                            unname(full_mu[, check, drop = FALSE]), tolerance = 1e-10)),
            isTRUE(all.equal(unname(fit$tau_f_predict[, nrow(grid) + seq_along(check), drop = FALSE]),
                            unname(full_tau[, check, drop = FALSE]), tolerance = 1e-10)))
  rm(fit); gc(verbose = FALSE)
  tr <- which(sub$z == 1)
  att_serr <- numeric(nrow(full_mu))
  for (first in seq.int(1L, nrow(full_mu), by = 20L)) {
    rows <- seq.int(first, min(first + 19L, nrow(full_mu)))
    r0 <- exp(full_mu[rows, tr, drop = FALSE])
    r1 <- exp(full_mu[rows, tr, drop = FALSE] + full_tau[rows, tr, drop = FALSE])
    att_serr[rows] <- as.numeric((r1 %*% sub$pyr[tr]) / (r0 %*% sub$pyr[tr])) - 1
  }
  rm(full_mu, full_tau, r0, r1); gc(verbose = FALSE)
  main <- as_tibble(c(list(label = "solid:high"), qsum(att_serr, "att_serr")))
  out <- list(chain = ch, seed = SEED + 1000L * ch,
       nburn = PROFILE_NBURN, nsim = PROFILE_NSIM, grid = grid,
       source_identity = profile_identity,
       logmu = logmu, logrr = logrr, kappa = kappa,
       main_att_serr = att_serr, main_summary = main)
  saveRDS(out, paste0(checkpoint, ".tmp"))
  stopifnot(file.rename(paste0(checkpoint, ".tmp"), checkpoint))
  stamp("profile chain ", ch, " checkpoint saved")
  out
}

# A single-chain invocation lets an external memory-limited runner start a
# fresh R process for each chain and keep an already completed checkpoint.
chain_request <- Sys.getenv("CBCF_PROFILE_CHAIN")
if (nzchar(chain_request)) {
  ch <- as.integer(chain_request)
  stopifnot(ch %in% 1:2)
  invisible(profile_chain(ch))
  quit(save = "no", status = 0)
}

stamp("profile comparison: sequential extended profile chains (no parallel workers)")
if (identical(Sys.getenv("CBCF_PROFILE_SUMMARISE_ONLY"), "1")) {
  stopifnot(all(file.exists(file.path(OUTDIR, paste0("profile_age30_70_chain", 1:2, ".rds")))))
}
chains <- lapply(1:2, function(ch) {
  out <- profile_chain(ch)
  gc(verbose = FALSE)
  out
})
if (any(vapply(chains, function(x) inherits(x, "try-error") || is.null(x$logrr), logical(1))))
  stop("A profile chain failed")
saveRDS(list(chains = chains, grid = grid, nburn = PROFILE_NBURN, nsim = PROFILE_NSIM,
             dose_group = "0.5-<2 Gy versus in-city <5 mGy",
             city_weighting = "sex-specific exposed-group person-years; counterfactual rates pooled before taking their ratio"),
        file.path(OUTDIR, "profile_age30_70_draws.rds"))

pf <- readRDS(file.path(OUTDIR, "preston_fit.rds"))$preston_par
pt <- read_csv(file.path(TAB, "T08_preston_err_ear_refit.csv"), show_col_types = FALSE)
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
dir.create(file.path(OUTDIR, "logs"), recursive = TRUE, showWarnings = FALSE)
writeLines(capture.output(sessionInfo()), file.path(OUTDIR, "logs", "sessionInfo_10.txt"))
print(as.data.frame(res), digits = 4)
print(as.data.frame(diag), digits = 4)
stamp("10 profile comparison done -> ", TAB)
