#!/usr/bin/env Rscript
# =====================================================================
# 01  Data validation, descriptives, and the PARAMETRIC comparators.
#
#  (a) reproduce Preston et al. (2007) Table 4 exactly from lssinc07.csv
#      -> proof that the data handling matches the canonical analysis;
#  (b) crude rates / overdispersion / zero-inflation diagnostics;
#  (c) covariate balance across the four dose bands;
#  (d) Poisson & NB GLMs with categorical dose  -> the ATE yardstick;
#  (e) Preston's own ERR and EAR models, refitted here by maximum
#      Poisson likelihood, so that the CountBCF CATE surface can be
#      compared against what the canonical parametric form can represent.
# =====================================================================

source(file.path(Sys.getenv("CBCF_ROOT", "."), "Applied_Study/code/common.R"))
set.seed(SEED)

d <- read_lss()

# =====================================================================
# (a) Preston Table 4 replication
# =====================================================================
stamp("(a) Preston 2007 Table 4 replication")

preston_bands <- function(x) {
  case_when(x == 1 ~ "NIC", x == 2 ~ "<0.005", x <= 7 ~ "0.005-0.1",
            x <= 11 ~ "0.1-0.2", x <= 14 ~ "0.2-0.5", x <= 16 ~ "0.5-1",
            x <= 20 ~ "1-2", TRUE ~ "2-4")
}
published <- tribble(
  ~stratum,     ~p_cases, ~p_pyr,   ~p_people,
  "NIC",           3994,   680744,   25247,
  "<0.005",        5603,   918199,   35545,
  "0.005-0.1",     4406,   729603,   27789,
  "0.1-0.2",        968,   145925,    5527,
  "0.2-0.5",       1144,   153846,    5935,
  "0.5-1",          688,    81251,    3173,
  "1-2",            460,    41412,    1647,
  "2-4",            185,    13711,     564,
  "Male",          7969,  1040278,   42902,
  "Female",        9479,  1724452,   62525,
  "Hiroshima",    12931,  1967601,   73430,
  "Nagasaki",      4517,   797129,   31997,
  "TOTAL",        17448,  2764730,  105427)

kd <- d %>% filter(dcatn > 0)
obs <- bind_rows(
  kd %>% mutate(stratum = preston_bands(dcatn)) %>% group_by(stratum) %>%
    summarise(cases = sum(solid), pyr = sum(pyr), people = sum(subjects), .groups = "drop"),
  kd %>% mutate(stratum = ifelse(sex == 1, "Male", "Female")) %>% group_by(stratum) %>%
    summarise(cases = sum(solid), pyr = sum(pyr), people = sum(subjects), .groups = "drop"),
  kd %>% mutate(stratum = ifelse(city == 1, "Hiroshima", "Nagasaki")) %>% group_by(stratum) %>%
    summarise(cases = sum(solid), pyr = sum(pyr), people = sum(subjects), .groups = "drop"),
  tibble(stratum = "TOTAL", cases = sum(kd$solid), pyr = sum(kd$pyr),
         people = sum(kd$subjects)))

valid <- published %>% left_join(obs, by = "stratum") %>%
  transmute(stratum, cases, published_cases = p_cases, cases_ok = cases == p_cases,
            pyr = round(pyr), published_pyr = p_pyr,
            pyr_reldiff = round((pyr - p_pyr) / p_pyr, 6),
            people, published_people = p_people, people_ok = people == p_people)
write_csv(valid, file.path(TAB, "T01_preston_table4_validation.csv"))
print(as.data.frame(valid))
stamp("  cases match: ", all(valid$cases_ok), " | people match: ", all(valid$people_ok),
      " | max |PY rel. diff|: ", signif(max(abs(valid$pyr_reldiff)), 3))

# solid-cancer aggregate composition (does `solid` include non-melanoma skin?)
disj <- c("oralca", "digestca", "respca", "femgenca", "malgenca", "urinca")
named <- rowSums(d[, disj]) + d$breast + d$thyroid + d$cnsca + d$bone +
         d$connect + d$thymus + d$othsol
# Preston Table 2 is built from 11 mutually exclusive site groups; compare
# the known-dose totals directly with the published counts.
grp <- c("oralca", "digestca", "respca", "skinbone", "breast", "femgenca",
         "malgenca", "urinca", "cnsca", "thyroid", "othsol")
pub <- c(oralca = 277, digestca = 10052, respca = 2001, skinbone = 400,
         breast = 1082, femgenca = 1457, malgenca = 420, urinca = 741,
         cnsca = 281, thyroid = 471, othsol = 266)
comp <- tibble(
  group = c(grp, "sum of the 11 groups", "solid (RERF all-solid total)",
            "non-melanoma skin (inside skinbone)", "melanoma (inside skinbone)",
            "stomach", "lung", "colon", "liver", "prostate"),
  cases_known_dose = c(vapply(grp, function(v) sum(kd[[v]]), numeric(1)),
                       sum(rowSums(kd[, grp])), sum(kd$solid), sum(kd$nmskin),
                       sum(kd$skinbone) - sum(kd$nmskin) - sum(kd$bone) - sum(kd$connect),
                       sum(kd$stomach), sum(kd$lung), sum(kd$colon), sum(kd$liver),
                       sum(kd$prost)),
  preston_table2 = c(pub[grp], sum(pub), 17448, 330, 17, 4730, 1759, 1516, 1494, 387))
comp$matches <- comp$cases_known_dose == comp$preston_table2
write_csv(comp, file.path(TAB, "T02_allsolid_composition.csv"))
print(as.data.frame(comp))
stamp("  all-solid includes non-melanoma skin: ",
      all(d$solid >= named + d$nmskin), " (cells violating: ",
      sum(d$solid < named + d$nmskin), ")")

# =====================================================================
# (b) analysis-set descriptives, overdispersion, zeros
# =====================================================================
stamp("(b) descriptives / dispersion")

ana <- d %>% filter(dcatn > 0, band %in% c("NIC", BANDS))
eda <- ana %>% group_by(band) %>%
  summarise(n_cells = n(), n_subjects = sum(subjects), person_years = sum(pyr),
            solid_cases = sum(solid), stomach_cases = sum(stomach),
            solid_rate_1e4 = sum(solid) / sum(pyr) * 1e4,
            stomach_rate_1e4 = sum(stomach) / sum(pyr) * 1e4,
            mean_colon_dose_Gy = sum(cola02w10 * pyr) / sum(pyr),
            max_colon_dose_Gy = max(cola02w10),
            pct_zero_cells_solid = mean(solid == 0),
            var_mean_solid = var(solid) / mean(solid),
            pct_zero_cells_stomach = mean(stomach == 0),
            var_mean_stomach = var(stomach) / mean(stomach),
            mean_agex = sum(agex * pyr) / sum(pyr),
            mean_age = sum(age * pyr) / sum(pyr),
            pct_female = sum(pyr * (sex == 2)) / sum(pyr),
            pct_hiroshima = sum(pyr * (city == 1)) / sum(pyr), .groups = "drop")
write_csv(eda, file.path(TAB, "T03_band_summary.csv"))
print(as.data.frame(eda), digits = 4)

# crude rates by age at exposure x band x sex  (Preston Table 7 analogue)
crude <- ana %>% filter(band != "NIC") %>%
  group_by(sex_lab = ifelse(sex == 1, "Male", "Female"), agex_b = agex_band(agex), band) %>%
  summarise(cases = sum(solid), pyr = sum(pyr),
            rate_1e4 = sum(solid) / sum(pyr) * 1e4, .groups = "drop") %>%
  arrange(sex_lab, agex_b, band)
crude_rr <- crude %>% group_by(sex_lab, agex_b) %>%
  mutate(crude_RR = rate_1e4 / rate_1e4[band == "reference"]) %>% ungroup()
write_csv(crude_rr, file.path(TAB, "T04_crude_rates_by_agex_sex_band.csv"))

# overdispersion diagnostics on the analysis set (in-city, known dose, <2 Gy)
ic <- ana %>% filter(band %in% BANDS)
F_BASE <- outcome ~ band + sex_f * (ns(log(age), 4) + ns(agex, 3)) + city_f
ic2 <- ic %>% mutate(outcome = solid, sex_f = factor(sex, 1:2, c("Male", "Female")),
                     city_f = factor(city, 1:2, c("Hiroshima", "Nagasaki")))
mp <- glm(F_BASE, poisson("log"), ic2, offset = log(pyr))
mnb <- glm.nb(update(F_BASE, . ~ . + offset(log(pyr))), ic2)
disp <- tibble(
  quantity = c("n cells", "cases", "person-years", "% zero cells",
               "Poisson Pearson X2 / df", "Poisson deviance / df",
               "NB theta (size)", "NB implied var/mean at mean fitted rate",
               "AIC Poisson", "AIC NB", "LRT NB vs Poisson (p)"),
  value = c(nrow(ic2), sum(ic2$outcome), sum(ic2$pyr), mean(ic2$outcome == 0),
            sum(residuals(mp, "pearson")^2) / mp$df.residual,
            mp$deviance / mp$df.residual, mnb$theta,
            1 + mean(fitted(mnb)) / mnb$theta, AIC(mp), AIC(mnb),
            pchisq(2 * (logLik(mnb) - logLik(mp)), df = 1, lower.tail = FALSE) / 2))
write_csv(disp, file.path(TAB, "T05_dispersion_diagnostics.csv"))
print(as.data.frame(disp), digits = 5)

# =====================================================================
# (c) covariate balance across arms
# =====================================================================
wm  <- function(x, w) sum(x * w) / sum(w)
wsd <- function(x, w) sqrt(sum(w * (x - wm(x, w))^2) / sum(w))
bal <- ic %>% group_by(band) %>%
  summarise(person_years = sum(pyr),
            agex_mean = wm(agex, pyr), agex_sd = wsd(agex, pyr),
            age_mean  = wm(age, pyr),  age_sd  = wsd(age, pyr),
            female    = wm(as.numeric(sex == 2), pyr),
            hiroshima = wm(as.numeric(city == 1), pyr),
            gdist_mean = wm(gdist, pyr), .groups = "drop")
# standardised differences vs the reference arm
ref <- bal %>% filter(band == "reference")
psd <- function(p, q) sqrt((p * (1 - p) + q * (1 - q)) / 2)
bal <- bal %>% mutate(
  smd_agex = (agex_mean - ref$agex_mean) / sqrt((agex_sd^2 + ref$agex_sd^2) / 2),
  smd_age  = (age_mean - ref$age_mean) / sqrt((age_sd^2 + ref$age_sd^2) / 2),
  smd_female    = (female - ref$female) / psd(female, ref$female),
  smd_hiroshima = (hiroshima - ref$hiroshima) / psd(hiroshima, ref$hiroshima))
write_csv(bal, file.path(TAB, "T06_covariate_balance.csv"))
print(as.data.frame(bal), digits = 4)

# =====================================================================
# (d) Poisson / NB GLMs with categorical dose -- the ATE yardstick
# =====================================================================
stamp("(d) parametric GLM ATE")

F_MIN <- outcome ~ band
F_ADJ <- outcome ~ band + sex_f * (ns(log(age), 4) + ns(agex, 3)) + city_f
F_INT <- outcome ~ band * sex_f + band:ns(agex, 3) + band:ns(log(age), 3) +
                   sex_f * (ns(log(age), 4) + ns(agex, 3)) + city_f

marginal_effect <- function(model, dat, treated, nsim = 4000, scale = 1e4) {
  b <- coef(model); V <- vcov(model)
  b <- b[!is.na(b)]
  keep <- intersect(names(b), rownames(V))       # drop aliased terms
  b <- b[keep]; V <- V[keep, keep, drop = FALSE]
  tt <- delete.response(terms(model))
  d0 <- dat; d0$band <- factor("reference", levels = levels(dat$band))
  d1 <- dat; d1$band <- factor(treated, levels = levels(dat$band))
  X0 <- model.matrix(tt, d0)[, names(b), drop = FALSE]
  X1 <- model.matrix(tt, d1)[, names(b), drop = FALSE]
  w  <- dat$pyr / sum(dat$pyr)
  B  <- MASS::mvrnorm(nsim, b, V)
  r0 <- exp(X0 %*% t(B)); r1 <- exp(X1 %*% t(B))      # n x nsim
  ear_d <- as.numeric(t(r1 - r0) %*% w) * scale
  rr_d  <- as.numeric((t(r1) %*% w) / (t(r0) %*% w))
  tibble(treated = treated,
         base_rate_1e4 = weighted.mean(exp(X0 %*% b), dat$pyr) * scale,
         ear = mean(ear_d), ear_lo = quantile(ear_d, .025), ear_hi = quantile(ear_d, .975),
         rr = mean(rr_d), rr_lo = quantile(rr_d, .025), rr_hi = quantile(rr_d, .975))
}

param_all <- list()
for (ep in c("solid", "stomach")) {
  dd <- ic %>% mutate(outcome = .data[[ep]],
                      sex_f = factor(sex, 1:2, c("Male", "Female")),
                      city_f = factor(city, 1:2, c("Hiroshima", "Nagasaki")),
                      band = droplevels(factor(as.character(band), levels = BANDS)))
  models <- list(
    `Poisson (minimal)`  = glm(F_MIN, poisson("log"), dd, offset = log(pyr)),
    `Poisson (adjusted)` = glm(F_ADJ, poisson("log"), dd, offset = log(pyr)),
    `NB (minimal)`       = glm.nb(update(F_MIN, . ~ . + offset(log(pyr))), dd),
    `NB (adjusted)`      = glm.nb(update(F_ADJ, . ~ . + offset(log(pyr))), dd),
    `NB (dose x modifiers)` = glm.nb(update(F_INT, . ~ . + offset(log(pyr))), dd))
  for (nm in names(models)) {
    m <- models[[nm]]
    param_all[[paste(ep, nm)]] <- bind_rows(lapply(
      c("low", "medium", "high"), function(lv) marginal_effect(m, dd, lv))) %>%
      mutate(endpoint = ep, model = nm, aic = AIC(m),
             overdisp = if (inherits(m, "negbin")) NA_real_ else
               sum(residuals(m, "pearson")^2) / m$df.residual, .before = 1)
  }
}
param_tab <- bind_rows(param_all)
write_csv(param_tab, file.path(TAB, "T07_parametric_ate.csv"))
print(as.data.frame(param_tab %>% filter(grepl("adjusted|minimal", model))), digits = 4)

# =====================================================================
# (e) Preston's ERR / EAR models refitted by maximum Poisson likelihood
# =====================================================================
# ERR model:  lambda = pyr * exp(eta0) * [1 + beta_s * d * exp(a_e(e-30) + a_a log(a/70))]
# EAR model:  lambda = pyr * [exp(eta0) + beta_s * d * exp(...) / 1e4]
# eta0 = gender-specific splines in log(attained age) and in age at exposure
#        (= birth cohort, since exposure date is common) + city.
stamp("(e) Preston ERR / EAR model refit")

pr <- ic %>% mutate(outcome = solid, dose = cola02w10,
                    sex_f = factor(sex, 1:2, c("Male", "Female")),
                    city_f = factor(city, 1:2, c("Hiroshima", "Nagasaki")))
X0 <- model.matrix(~ sex_f * (ns(log(age), 4) + ns(agex, 3)) + city_f, pr)
p0 <- ncol(X0)
yy <- pr$outcome; py <- pr$pyr; dz <- pr$dose
ee <- pr$agex; aa <- pr$age; fem <- as.numeric(pr$sex == 2)

make_nll <- function(kind) function(par) {
  eta0 <- as.numeric(X0 %*% par[1:p0])
  bm <- par[p0 + 1]; bf <- par[p0 + 2]; ae <- par[p0 + 3]; aa_p <- par[p0 + 4]
  bs <- ifelse(fem == 1, bf, bm)
  eps <- exp(ae * (ee - 30) + aa_p * log(aa / 70))
  if (kind == "ERR") {
    mult <- 1 + bs * dz * eps
    if (any(!is.finite(mult)) || any(mult <= 1e-8)) return(1e10)
    lam <- py * exp(eta0) * mult
  } else {
    lam <- py * (exp(eta0) + bs * dz * eps / 1e4)
    if (any(!is.finite(lam)) || any(lam <= 1e-12)) return(1e10)
  }
  -sum(dpois(yy, lam, log = TRUE))
}

fit_preston <- function(kind) {
  m_start <- glm(outcome ~ sex_f * (ns(log(age), 4) + ns(agex, 3)) + city_f,
                 poisson("log"), pr, offset = log(pyr))
  st <- c(coef(m_start), if (kind == "ERR") c(0.35, 0.58) else c(30, 50),
          -0.017, -1.6)
  o <- optim(st, make_nll(kind), method = "BFGS",
             control = list(maxit = 3000, reltol = 1e-11))
  H <- optimHess(o$par, make_nll(kind))
  se <- suppressWarnings(sqrt(diag(solve(H))))
  k <- p0
  tibble(model = kind,
         parameter = c("beta_male (at e=30, a=70)", "beta_female (at e=30, a=70)",
                       "alpha_agex (per year)", "alpha_age (log-age power)",
                       "% decline in excess per decade of age at exposure",
                       "female:male ratio", "optim convergence code (0 = success)", "logLik", "AIC"),
         estimate = c(o$par[k + 1], o$par[k + 2], o$par[k + 3], o$par[k + 4],
                      100 * (1 - exp(10 * o$par[k + 3])),
                      o$par[k + 2] / o$par[k + 1], o$convergence,
                      -o$value, 2 * (o$value + length(o$par))),
         se = c(se[k + 1], se[k + 2], se[k + 3], se[k + 4], NA, NA, NA, NA, NA)) %>%
    mutate(lo = estimate - 1.645 * se, hi = estimate + 1.645 * se) %>%   # 90% CI, as Preston
    list(tab = ., par = o$par)
}

pres <- list(ERR = fit_preston("ERR"), EAR = fit_preston("EAR"))
preston_tab <- bind_rows(pres$ERR$tab, pres$EAR$tab)
write_csv(preston_tab, file.path(TAB, "T08_preston_err_ear_refit.csv"))
print(as.data.frame(preston_tab), digits = 4)

# model-implied ERR / EAR surface per cell -> used in 02 to isolate what is
# genuinely NEW in the CountBCF CATE surface
par_err <- pres$ERR$par; par_ear <- pres$EAR$par
pr$preston_err <- ifelse(fem == 1, par_err[p0 + 2], par_err[p0 + 1]) * dz *
  exp(par_err[p0 + 3] * (ee - 30) + par_err[p0 + 4] * log(aa / 70))
pr$preston_ear <- ifelse(fem == 1, par_ear[p0 + 2], par_ear[p0 + 1]) * dz *
  exp(par_ear[p0 + 3] * (ee - 30) + par_ear[p0 + 4] * log(aa / 70))
write_csv(pr %>% dplyr::select(band, sex, city, agxcat, agecat, time, agex, age,
                               pyr, dose, outcome, preston_err, preston_ear),
          file.path(TAB, "T09_preston_model_cell_surface.csv.gz"))

# categorical ("non-parametric") dose response: one parameter per dcatn
dcat_fit <- glm(solid ~ factor(dcatn) + sex_f * (ns(log(age), 4) + ns(agex, 3)) + city_f,
                poisson("log"), pr, offset = log(pyr))
cf <- coef(summary(dcat_fit))
dcr <- tibble(term = rownames(cf)) %>%
  filter(grepl("factor\\(dcatn\\)", term)) %>%
  mutate(dcatn = as.integer(sub(".*\\)", "", term)),
         logRR = cf[term, 1], se = cf[term, 2],
         RR = exp(logRR), RR_lo = exp(logRR - 1.96 * se), RR_hi = exp(logRR + 1.96 * se)) %>%
  left_join(pr %>% group_by(dcatn) %>%
              summarise(mean_dose_Gy = sum(dose * pyr) / sum(pyr),
                        cases = sum(outcome), pyr = sum(pyr), .groups = "drop"),
            by = "dcatn") %>%
  mutate(ERR_per_Gy = (RR - 1) / mean_dose_Gy)
write_csv(dcr, file.path(TAB, "T10_categorical_dose_response.csv"))
print(as.data.frame(dcr), digits = 4)

saveRDS(list(preston_par = list(ERR = par_err, EAR = par_ear, p0 = p0)),
        file.path(OUTDIR, "preston_fit.rds"))
stamp("01 done -> ", TAB)
