#!/usr/bin/env Rscript
# =====================================================================
# 05  Publication figures, built from the tables written by 01-04.
# =====================================================================

source(file.path(Sys.getenv("CBCF_ROOT", "."), "Applied_Study_v2/code/common.R"))
suppressMessages({library(ggplot2); library(scales)})

theme_set(theme_bw(base_size = 12))
SEXCOL  <- c(Male = "#2C6E9B", Female = "#B5442A")
CITYCOL <- c(Hiroshima = "#3F6E4B", Nagasaki = "#8C5A9E")
BANDLAB <- c(low = "5-<100 mGy", medium = "100-<500 mGy", high = "500-<2000 mGy")
gg <- function(p, f, w = 7.5, h = 5) ggsave(file.path(PLT, f), p, width = w, height = h, dpi = 300)

rd <- function(f) read_csv(file.path(TAB, f), show_col_types = FALSE)

cells   <- rd("T16_countbcf_cate_cells.csv.gz") %>% add_bands()
ate     <- rd("T12_countbcf_ate.csv")
conv    <- rd("T13_ate_convergence.csv")
sg      <- rd("T15_subgroup_cate.csv")
newness <- rd("T17_vs_preston_model.csv")
dcr     <- rd("T10_categorical_dose_response.csv")
ovl     <- rd("T24_overlap_cells.csv.gz")
sensA   <- rd("T20_sensitivity_tau_prior.csv")
sensB   <- rd("T21_sensitivity_data.csv")

cells <- cells %>% mutate(
  endpoint = sub(":.*", "", label), treated = sub(".*:", "", label),
  band_lab = factor(BANDLAB[treated], levels = BANDLAB))

# ---------------------------------------------------------------------
# F1  ATE convergence: parametric GLM vs CountBCF
# ---------------------------------------------------------------------
f1d <- conv %>% mutate(treated = factor(treated, levels = c("low", "medium", "high")),
                       endpoint = factor(endpoint, c("solid", "stomach"),
                                         c("All solid cancers", "Stomach cancer")))
p1 <- ggplot(f1d, aes(treated, ear, colour = method, shape = method)) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey55") +
  geom_pointrange(aes(ymin = ear_lo, ymax = ear_hi), position = position_dodge(.55)) +
  facet_wrap(~endpoint, scales = "free_y") +
  scale_x_discrete(labels = BANDLAB) +
  labs(x = NULL, y = "Excess absolute rate  (cases / 10,000 PY)", colour = NULL, shape = NULL,
       title = "ATE agreement: parametric GLMs vs CountBCF",
       subtitle = "Marginal excess rate vs the 0-<5 mGy in-city reference, 95% intervals") +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 20, hjust = 1))
gg(p1, "F01_ate_convergence.png", 9, 5)

# ---------------------------------------------------------------------
# F2  categorical dose response (parametric) with the CountBCF band ERRs
# ---------------------------------------------------------------------
band_err <- ate %>% filter(endpoint == "solid") %>%
  transmute(mean_dose_Gy = dose_contrast_Gy, ERR = att_err,
            lo = att_err_lo, hi = att_err_hi, what = "CountBCF band contrast")
p2 <- ggplot(dcr, aes(mean_dose_Gy, RR - 1)) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey55") +
  geom_pointrange(aes(ymin = RR_lo - 1, ymax = RR_hi - 1), colour = "grey35",
                  size = .25, alpha = .8) +
  geom_smooth(method = "lm", formula = y ~ x - 1, se = TRUE, colour = "#2C6E9B", fill = "#2C6E9B22") +
  geom_pointrange(data = band_err, aes(mean_dose_Gy, ERR, ymin = lo, ymax = hi),
                  colour = "#B5442A", size = .7) +
  labs(x = "Person-year weighted colon dose (Gy)", y = "Excess relative risk",
       title = "Dose response for all solid cancers",
       subtitle = paste("Grey: 18 dose-category Poisson estimates with a linear fit through the origin.",
                        "Red: the three CountBCF band contrasts.", sep = "\n"))
gg(p2, "F02_dose_response.png", 7.5, 5)

# ---------------------------------------------------------------------
# F3  ERR by age at exposure and sex (the canonical gradient)
# ---------------------------------------------------------------------
p3 <- ggplot(cells %>% filter(z == 1, endpoint == "solid"),
             aes(agex, err, colour = sex_f)) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey55") +
  geom_point(alpha = .10, size = .5) +
  geom_smooth(method = "loess", span = .8, se = TRUE, formula = y ~ x) +
  facet_wrap(~band_lab) +
  scale_colour_manual(values = SEXCOL) +
  labs(x = "Age at exposure (years)", y = "Excess relative risk", colour = NULL,
       title = "CountBCF excess relative risk by age at exposure, all solid cancers",
       subtitle = "Posterior mean per person-year cell; no interaction terms were specified") +
  theme(legend.position = "bottom")
gg(p3, "F03_err_by_agex_sex.png", 10, 4.6)

# ---------------------------------------------------------------------
# F4  EAR by attained age and sex
# ---------------------------------------------------------------------
p4 <- ggplot(cells %>% filter(z == 1, endpoint == "solid"),
             aes(age, ear, colour = sex_f)) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey55") +
  geom_point(alpha = .10, size = .5) +
  geom_smooth(method = "loess", span = .8, se = TRUE, formula = y ~ x) +
  facet_wrap(~band_lab, scales = "free_y") +
  scale_colour_manual(values = SEXCOL) +
  labs(x = "Attained age (years)", y = "Excess cases / 10,000 PY", colour = NULL,
       title = "CountBCF excess absolute rate by attained age, all solid cancers",
       subtitle = "The relative effect fades with age while the absolute burden grows") +
  theme(legend.position = "bottom")
gg(p4, "F04_ear_by_age_sex.png", 10, 4.6)

# ---------------------------------------------------------------------
# F5  subgroup CATE minus ATE -- the separation test
# ---------------------------------------------------------------------
f5 <- sg %>% filter(enough, by %in% c("sex_f", "agex_b", "age_b", "city_f"),
                    grepl("^solid", label)) %>%
  mutate(treated = sub(".*:", "", label),
         band_lab = factor(BANDLAB[treated], levels = BANDLAB),
         by = recode(by, sex_f = "Sex", agex_b = "Age at exposure",
                     age_b = "Attained age", city_f = "City"),
         subgroup = factor(subgroup, levels = rev(unique(subgroup))))
p5 <- ggplot(f5, aes(subgroup, d_err, colour = separates_err)) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey55") +
  geom_pointrange(aes(ymin = d_err_lo, ymax = d_err_hi), size = .35) +
  coord_flip() +
  facet_grid(by ~ band_lab, scales = "free_y", space = "free_y") +
  scale_colour_manual(values = c(`TRUE` = "#B5442A", `FALSE` = "grey45"),
                      labels = c(`TRUE` = "separates from ATE", `FALSE` = "not distinguishable"),
                      name = NULL) +
  labs(x = NULL, y = "Subgroup ERR minus overall ERR (posterior, 95% CrI)",
       title = "Which subgroups are genuinely different from the average effect?",
       subtitle = "All solid cancers. A subgroup counts only if the CrI on the DIFFERENCE excludes zero") +
  theme(legend.position = "bottom")
gg(p5, "F05_subgroup_separation.png", 10, 7)

# ---------------------------------------------------------------------
# F6  CountBCF vs the refitted Preston ERR model
# ---------------------------------------------------------------------
f6 <- newness %>% filter(grepl("^solid", label), cases >= 20,
                         by %in% c("sex_f", "agex_b", "age_b", "city_f")) %>%
  mutate(treated = sub(".*:", "", label),
         band_lab = factor(BANDLAB[treated], levels = BANDLAB))
p6 <- ggplot(f6, aes(preston_err, cbcf_err, colour = departs)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey55") +
  geom_pointrange(aes(ymin = cbcf_err_lo, ymax = cbcf_err_hi), size = .3) +
  geom_text(aes(label = subgroup), size = 2.4, hjust = -0.12, show.legend = FALSE) +
  facet_wrap(~band_lab, scales = "free") +
  scale_colour_manual(values = c(`TRUE` = "#B5442A", `FALSE` = "grey45"),
                      labels = c(`TRUE` = "departs", `FALSE` = "agrees"), name = NULL) +
  expand_limits(x = 0, y = 0) +
  labs(x = "ERR implied by the refitted Preston (2007) model",
       y = "CountBCF ERR (95% CrI)",
       title = "What is new relative to the canonical parametric ERR model?",
       subtitle = "All solid cancers, subgroup means over treated person-year cells") +
  theme(legend.position = "bottom")
gg(p6, "F06_vs_preston_model.png", 10, 5)

# ---------------------------------------------------------------------
# F7  CATE surface over age at exposure x attained age
# ---------------------------------------------------------------------
hm <- cells %>% filter(z == 1, endpoint == "solid", treated == "high") %>%
  mutate(ax = cut(agex, seq(0, 80, 10)), ay = cut(age, seq(0, 100, 10))) %>%
  filter(!is.na(ax), !is.na(ay)) %>%
  group_by(sex_f, ax, ay) %>%
  summarise(err = weighted.mean(err, pyr), ear = weighted.mean(ear, pyr),
            pyr = sum(pyr), .groups = "drop") %>% filter(pyr > 200)
p7a <- ggplot(hm, aes(ax, ay, fill = err)) + geom_tile() + facet_wrap(~sex_f) +
  scale_fill_gradient2(low = "#2C6E9B", mid = "white", high = "#B5442A",
                       midpoint = 0, name = "ERR") +
  labs(x = "Age at exposure", y = "Attained age",
       title = "CountBCF ERR surface, all solid cancers, 500-<2000 mGy") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
gg(p7a, "F07_err_surface.png", 9, 5)
p7b <- ggplot(hm, aes(ax, ay, fill = ear)) + geom_tile() + facet_wrap(~sex_f) +
  scale_fill_gradient2(low = "#2C6E9B", mid = "white", high = "#B5442A",
                       midpoint = 0, name = "EAR /1e4 PY") +
  labs(x = "Age at exposure", y = "Attained age",
       title = "CountBCF excess absolute rate surface, all solid cancers, 500-<2000 mGy") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
gg(p7b, "F08_ear_surface.png", 9, 5)

# ---------------------------------------------------------------------
# F9  overlap
# ---------------------------------------------------------------------
ovl2 <- ovl %>% mutate(arm = factor(z, 0:1, c("reference (0-<5 mGy)", "exposed band")),
                       band_lab = factor(BANDLAB[contrast], levels = BANDLAB))
p9 <- ggplot(ovl2, aes(pihat, fill = arm, colour = arm)) +
  geom_density(alpha = .35, adjust = 1.1) +
  facet_wrap(~band_lab, scales = "free") +
  scale_fill_manual(values = c("#2C6E9B", "#B5442A")) +
  scale_colour_manual(values = c("#2C6E9B", "#B5442A")) +
  labs(x = "Estimated propensity score", y = "Density", fill = NULL, colour = NULL,
       title = "Overlap by arm (spline logistic propensity model)",
       subtitle = "No cell has an estimated propensity near 0 or 1: positivity is not strained") +
  theme(legend.position = "bottom")
gg(p9, "F09_overlap.png", 10, 4.2)

# ---------------------------------------------------------------------
# F10  sensitivity of the featured quantities
# ---------------------------------------------------------------------
mk_sens <- function(df, ttl) {
  long <- bind_rows(
    df %>% transmute(spec, q = "overall ERR", est = att_err, lo = att_err_lo, hi = att_err_hi),
    df %>% transmute(spec, q = "female minus male ERR", est = d_female_male,
                     lo = d_female_male_lo, hi = d_female_male_hi),
    df %>% transmute(spec, q = "ERR: exposed <20 minus 40+", est = d_young_old,
                     lo = d_young_old_lo, hi = d_young_old_hi)) %>%
    mutate(spec = factor(spec, levels = rev(df$spec)))
  ggplot(long, aes(spec, est)) +
    geom_hline(yintercept = 0, linetype = 2, colour = "grey55") +
    geom_pointrange(aes(ymin = lo, ymax = hi), size = .35, colour = "#2C6E9B") +
    coord_flip() + facet_wrap(~q, scales = "free_x") +
    labs(x = NULL, y = "Posterior mean and 95% CrI", title = ttl)
}
gg(mk_sens(sensA, "Sensitivity to the regularisation of the moderating forest"),
   "F10_sensitivity_tau_prior.png", 11, 5)
gg(mk_sens(sensB, "Sensitivity to the data construction"),
   "F11_sensitivity_data.png", 12, 6)

stamp("05 done -> ", PLT)
