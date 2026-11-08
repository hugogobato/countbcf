#!/usr/bin/env Rscript
# =====================================================================
# 05  Publication figures, built from the tables written by 01-04.
# =====================================================================

source(file.path(Sys.getenv("CBCF_ROOT", "."), "Applied_Study/code/common.R"))
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
  transmute(treated, ERR = att_serr,
            lo = att_serr_lo, hi = att_serr_hi,
            what = "CountBCF exposed-group band ERR") %>%
  left_join(cells %>% filter(endpoint == "solid", z == 1) %>%
              group_by(treated) %>%
              summarise(mean_dose_Gy = weighted.mean(dose, pyr), .groups = "drop"),
            by = "treated")
# Preserve the plotted values so this comparison can be checked without
# reverse-engineering a raster image. Both series use mean absolute dose.
write_csv(band_err, file.path(TAB, "T42_dose_response_figure_points.csv"))
p2 <- ggplot(dcr, aes(mean_dose_Gy, RR - 1)) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey55") +
  geom_pointrange(aes(ymin = RR_lo - 1, ymax = RR_hi - 1), colour = "grey35",
                  size = .25, alpha = .8) +
  geom_smooth(method = "lm", formula = y ~ x - 1, se = FALSE, colour = "#2C6E9B") +
  geom_pointrange(data = band_err, aes(mean_dose_Gy, ERR, ymin = lo, ymax = hi),
                  colour = "#B5442A", size = .7) +
  labs(x = "Person-year weighted colon dose (Gy)", y = "Excess relative rate (ERR)",
       title = "Dose response for all solid cancers",
       subtitle = paste("Grey: dose-category Poisson estimates; blue: descriptive linear trend.",
                        "Red: CountBCF exposed-group ratios of counterfactual rates minus one (95% CrI).", sep = "\n"))
gg(p2, "F02_dose_response.png", 7.5, 5)

# ---------------------------------------------------------------------
# F3  ERR by age at exposure and sex (descriptive posterior-mean smooths)
# Do not compute loess standard-error ribbons: these are fitted cell means,
# not independent observed effects, and loess SE calculation is memory-heavy.
# ---------------------------------------------------------------------
p3 <- ggplot(cells %>% filter(z == 1, endpoint == "solid"),
             aes(agex, err, colour = sex_f)) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey55") +
  geom_point(alpha = .10, size = .5) +
  geom_smooth(method = "loess", span = .8, se = FALSE, formula = y ~ x) +
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
  geom_smooth(method = "loess", span = .8, se = FALSE, formula = y ~ x) +
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
# standardised footing throughout, matching the tables
p5 <- ggplot(f5, aes(subgroup, d_serr, colour = separates_serr)) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey55") +
  geom_pointrange(aes(ymin = d_serr_lo, ymax = d_serr_hi), size = .35) +
  coord_flip() +
  facet_grid(by ~ band_lab, scales = "free_y", space = "free_y") +
  scale_colour_manual(values = c(`TRUE` = "#B5442A", `FALSE` = "grey45"),
                      labels = c(`TRUE` = "interval excludes zero", `FALSE` = "interval contains zero"),
                      name = NULL) +
  labs(x = NULL, y = "Subgroup band ERR minus cohort-averaged band ERR (95% CrI)",
       title = "Subgroup departures from the cohort-averaged effect",
       subtitle = "All solid cancers; person-year-weighted ratios of rates. Intervals describe subgroup-minus-cohort differences.") +
  theme(legend.position = "bottom")
gg(p5, "F05_subgroup_separation.png", 10, 7)

# ---------------------------------------------------------------------
# F6  CountBCF vs the refitted Preston ERR model
# ---------------------------------------------------------------------
f6 <- newness %>% filter(grepl("^solid", label), cases >= 20,
                         by %in% c("sex_f", "agex_b", "age_b", "city_f")) %>%
  mutate(treated = sub(".*:", "", label),
         band_lab = factor(BANDLAB[treated], levels = BANDLAB),
         plot_label = case_when(by == "agex_b" ~ paste("Exposed", subgroup),
                                by == "age_b" ~ paste("Attained", subgroup),
                                TRUE ~ subgroup)) %>%
  arrange(desc(departs_serr))
p6 <- ggplot(f6, aes(preston_serr, cbcf_serr, colour = departs_serr)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey55") +
  geom_pointrange(aes(ymin = cbcf_serr_lo, ymax = cbcf_serr_hi), size = .3) +
  geom_text(data = f6 %>% filter(by %in% c("sex_f", "city_f") | departs_serr),
            aes(label = plot_label), size = 2.7, hjust = -0.08,
            check_overlap = TRUE, show.legend = FALSE) +
  facet_wrap(~band_lab, scales = "free") +
  scale_colour_manual(values = c(`TRUE` = "#B5442A", `FALSE` = "grey45"),
                      labels = c(`TRUE` = "excludes fitted comparator", `FALSE` = "contains fitted comparator"), name = NULL) +
  expand_limits(x = 0, y = 0) +
  scale_x_continuous(expand = expansion(mult = c(.05, .45))) +
  labs(x = "Exposed-group band ERR implied by the refitted Preston (2007) model",
       y = "CountBCF exposed-group band ERR (95% CrI)",
       title = "CountBCF and fitted Preston subgroup effects",
       subtitle = "All solid cancers; both surfaces use the same posterior baseline-rate weights") +
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
ovl2 <- ovl %>% mutate(arm = factor(z, 0:1, c("reference (0-<5 mGy)", "higher-dose group")),
                       band_lab = factor(BANDLAB[contrast], levels = BANDLAB))
p9 <- ggplot(ovl2, aes(pihat, fill = arm, colour = arm)) +
  geom_density(alpha = .35, adjust = 1.1) +
  facet_wrap(~band_lab, scales = "free") +
  scale_fill_manual(values = c("#2C6E9B", "#B5442A")) +
  scale_colour_manual(values = c("#2C6E9B", "#B5442A")) +
  labs(x = "Estimated propensity score", y = "Density", fill = NULL, colour = NULL,
       title = "Demographic propensity by dose group (spline logistic model)",
       subtitle = "Scores overlap for the measured demographic covariates; distance is omitted") +
  theme(legend.position = "bottom")
gg(p9, "F09_overlap.png", 10, 4.2)

# ---------------------------------------------------------------------
# F10  sensitivity of the featured quantities
# ---------------------------------------------------------------------
# All three panels use the ratio-of-rates summaries, so the figure and
# Tables 10-11 report the same quantities.
mk_sens <- function(df, ttl) {
  long <- bind_rows(
    df %>% transmute(spec, q = "exposed-group band ERR", est = att_serr,
                     lo = att_serr_lo, hi = att_serr_hi),
    df %>% transmute(spec, q = "female minus male ERR", est = d_serr_female_male,
                     lo = d_serr_female_male_lo, hi = d_serr_female_male_hi),
    df %>% transmute(spec, q = "ERR: exposed <20 minus 40+", est = d_serr_young_old,
                     lo = d_serr_young_old_lo, hi = d_serr_young_old_hi)) %>%
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
