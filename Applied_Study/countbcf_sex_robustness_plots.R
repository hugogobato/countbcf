## ---------------------------------------------------------------------------
## countbcf_sex_robustness_plots.R
## Figures for the co-author revision round (Bucket B sensitivity analyses).
## Reads the B1-B7 summary tables produced by countbcf_sex_robustness.R and
## renders one PNG per experiment into countbcf_applied_analysis/plots/.
## Style matches countbcf_applied_analysis.R (ggplot2, theme_minimal, 300 dpi).
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(ggplot2)
})

BASE   <- "countbcf_applied_analysis"
TABLES <- file.path(BASE, "tables")
PLOTS  <- file.path(BASE, "plots")
dir.create(PLOTS, showWarnings = FALSE, recursive = TRUE)

rd   <- function(f) read.csv(file.path(TABLES, f), stringsAsFactors = FALSE)
save <- function(p, file, w = 8, h = 5.2) {
  ggsave(file.path(PLOTS, file), p, width = w, height = h, dpi = 300)
  cat("wrote", file.path(PLOTS, file), "\n")
}

theme_set(theme_minimal(base_size = 12) +
            theme(panel.grid.minor = element_blank(),
                  plot.title    = element_text(face = "bold", size = 13),
                  plot.subtitle = element_text(colour = "grey35", size = 10),
                  legend.position = "top"))

sex_cols <- c(Male = "#2c7fb8", Female = "#d95f0e")

## ===========================================================================
## B1 - what the RERF all-solid endpoint is made of
## ===========================================================================
b1 <- rd("B1_allsolid_construction.csv")
b1 <- b1[!grepl("^\\(|^solid \\(", b1$group), ]            # drop the two total rows
b1$status <- ifelse(grepl("^NO", b1$in_all_solid), "Excluded (RERF convention)",
              ifelse(b1$group %in% c("breast", "thyroid"),
                     "Included – breast / thyroid", "Included – other sites"))
b1$group  <- factor(b1$group, levels = b1$group[order(b1$cases)])
b1$status <- factor(b1$status, levels = c("Included – breast / thyroid",
                                          "Included – other sites",
                                          "Excluded (RERF convention)"))

p1 <- ggplot(b1, aes(group, cases, fill = status)) +
  geom_col(width = .72) +
  geom_text(aes(label = cases), hjust = -0.15, size = 3.1, colour = "grey25") +
  coord_flip(clip = "off") +
  scale_fill_manual(values = c("Included – breast / thyroid" = "#d95f0e",
                               "Included – other sites"      = "#41ab5d",
                               "Excluded (RERF convention)"       = "#bdbdbd")) +
  scale_y_continuous(expand = expansion(mult = c(0, .12))) +
  labs(title = "Composition of the RERF all-solid endpoint",
       subtitle = "Breast (316) and thyroid (140) ARE inside all-solid; non-melanoma skin is not.",
       x = NULL, y = "Incident cases (1958–1998)", fill = NULL) +
  theme(legend.position = "top")
save(p1, "B1_allsolid_construction.png", w = 8.5, h = 5.4)

## ===========================================================================
## B2 - breast & thyroid fit on their OWN organ dose (sex-specific ERR)
## ===========================================================================
g <- rd("B2_breast_thyroid_glm_err.csv")          # err + 95% CI (GLM)
c2 <- rd("B2_breast_thyroid_cbcf_sexerr.csv")     # CountBCF PY-weighted ERR
g$site  <- tools::toTitleCase(g$label)
g$Sex   <- tools::toTitleCase(g$sexlab)
c2$site <- tools::toTitleCase(c2$label)
c2$Sex  <- tools::toTitleCase(c2$sexlab)

p2 <- ggplot(g, aes(site, err, colour = Sex)) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey55") +
  geom_pointrange(aes(ymin = err_lo, ymax = err_hi),
                  position = position_dodge(width = .5), size = .7) +
  geom_point(data = c2, aes(site, err_pyw, colour = Sex),
             position = position_dodge(width = .5), shape = 4, size = 3, stroke = 1.1) +
  scale_colour_manual(values = sex_cols) +
  labs(title = "Breast & thyroid ERR on their own organ dose",
       subtitle = "Pointrange = GLM ERR ± 95% CI; × = CountBCF ERR. Excess large for both sexes.",
       x = NULL, y = "Excess relative risk (ERR)", colour = NULL)
save(p2, "B2_breast_thyroid_err.png", w = 7.5, h = 5)

## ===========================================================================
## B3 - all-solid WITH vs WITHOUT breast/thyroid: sex ERR & female/male ratio
## ===========================================================================
b3 <- rd("B3_solid_excl_bt_sexratio.csv")
defs <- c(solid_full = "All-solid\n(full)",
          solid_excl_breast_thy = "All-solid\nexcl. breast+thyroid")
err_long <- rbind(
  data.frame(def = defs[b3$label], Sex = "Male",   err = b3$err_male_pyw),
  data.frame(def = defs[b3$label], Sex = "Female", err = b3$err_female_pyw))
err_long$def <- factor(err_long$def, levels = defs)
ratio_lab <- data.frame(def = factor(defs[b3$label], levels = defs),
                        ratio = b3$ratio_fm_pyw,
                        y = pmax(b3$err_male_pyw, b3$err_female_pyw) + .07)

p3 <- ggplot(err_long, aes(def, err, fill = Sex)) +
  geom_col(position = position_dodge(width = .7), width = .62) +
  geom_text(data = ratio_lab,
            aes(def, y, label = sprintf("F/M ratio = %.2f", ratio)),
            inherit.aes = FALSE, size = 3.4, fontface = "bold", colour = "grey25") +
  scale_fill_manual(values = sex_cols) +
  scale_y_continuous(expand = expansion(mult = c(0, .15))) +
  labs(title = "Removing breast & thyroid barely moves the sex ratio",
       subtitle = "Female/male ERR ratio shifts only 0.91 → 0.92 — the sites are not what equalises the sexes.",
       x = NULL, y = "ERR (PY-weighted)", fill = NULL)
save(p3, "B3_solid_excl_bt_sexerr.png", w = 7.5, h = 5)

## ===========================================================================
## B4 - design-choice sensitivity: ERR by sex + female/male ratio per variant
## ===========================================================================
b4 <- rd("B4_design_sensitivity_sexratio.csv")
vlab <- c(baseline_cap2_b.1.5 = "Baseline\n(cap 2 Gy)",
          cap4Gy              = "Cap 4 Gy",
          bands_hi1Gy         = "High band\n≥1 Gy",
          bands_low.05        = "Low band\n0.05 Gy")
b4$variant <- factor(vlab[b4$label], levels = vlab)

err4 <- rbind(
  data.frame(variant = b4$variant, series = "Male",   value = b4$err_male_pyw,
             panel = "ERR (PY-weighted)"),
  data.frame(variant = b4$variant, series = "Female", value = b4$err_female_pyw,
             panel = "ERR (PY-weighted)"))
rat4 <- data.frame(variant = b4$variant, series = "Female / Male",
                   value = b4$ratio_fm_pyw, panel = "Female / Male ERR ratio")
b4l <- rbind(err4, rat4)
b4l$panel <- factor(b4l$panel, levels = c("ERR (PY-weighted)", "Female / Male ERR ratio"))
href <- data.frame(panel = factor("Female / Male ERR ratio",
                                  levels = levels(b4l$panel)), y = 1)

p4 <- ggplot(b4l, aes(variant, value, fill = series)) +
  geom_col(position = position_dodge(width = .72), width = .64) +
  geom_hline(data = href, aes(yintercept = y), linetype = 2, colour = "grey45") +
  facet_wrap(~panel, ncol = 1, scales = "free_y") +
  scale_fill_manual(values = c(Male = "#2c7fb8", Female = "#d95f0e",
                               "Female / Male" = "#756bb1")) +
  labs(title = "Female/male ERR ordering under alternative dose designs",
       subtitle = "Widening the cap or raising the high band flips the ratio above 1 (females higher).",
       x = NULL, y = NULL, fill = NULL)
save(p4, "B4_design_sensitivity.png", w = 8.5, h = 7)

## ===========================================================================
## B6 - E-values for the ATE relative rates
## ===========================================================================
b6 <- rd("B6_evalues.csv")
clab <- c(control_vs_moderate = "Control vs moderate",
          control_vs_high     = "Control vs high")
b6$contrast <- factor(clab[b6$contrast], levels = clab)
b6$endpoint <- tools::toTitleCase(b6$endpoint)

p6 <- ggplot(b6, aes(contrast, rr, colour = endpoint)) +
  geom_hline(yintercept = 1, linetype = 2, colour = "grey55") +
  geom_pointrange(aes(ymin = rr_lo, ymax = rr_hi),
                  position = position_dodge(width = .55), size = .7) +
  geom_text(aes(label = sprintf("E = %.2f (CI %.2f)", evalue_point, evalue_ci)),
            position = position_dodge(width = .55), vjust = -1.1, size = 3.1,
            show.legend = FALSE) +
  scale_colour_manual(values = c(Solid = "#2c7fb8", Stomach = "#d95f0e")) +
  scale_y_continuous(expand = expansion(mult = c(.05, .18))) +
  labs(title = "Relative rates and their E-values",
       subtitle = "E-value = confounder strength (RR scale) needed to explain away the estimate.",
       x = NULL, y = "Relative rate (high/moderate vs control)", colour = NULL)
save(p6, "B6_evalues.png", w = 7.8, h = 5)

## ===========================================================================
## B7 - site-stratified ERR by age at exposure (Little 2009 upturn check)
## ===========================================================================
b7 <- rd("B7_site_stratified_agex.csv")
ord <- c("<10", "10-20", "20-30", "30-40", "40-50", "50-60", "60+")
b7$band <- factor(b7$band, levels = ord)
b7$endpoint <- tools::toTitleCase(b7$endpoint)

p7 <- ggplot(b7, aes(band, err_pyw, colour = endpoint, group = endpoint)) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey55") +
  geom_line(linewidth = .9) +
  geom_point(size = 2.2) +
  scale_colour_brewer(palette = "Set1") +
  labs(title = "Site-stratified ERR by age at exposure",
       subtitle = "Within each site the trend is monotone – no within-site upturn (cf. Little 2009 aggregation artefact).",
       x = "Age at exposure (years)", y = "ERR (PY-weighted)", colour = NULL)
save(p7, "B7_site_stratified_agex.png", w = 8, h = 5)

cat("\nAll robustness figures written to", PLOTS, "\n")
