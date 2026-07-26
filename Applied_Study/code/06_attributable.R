#!/usr/bin/env Rscript
# =====================================================================
# 06  Excess cases and attributable fractions implied by CountBCF,
#     put next to Preston et al. (2007) Table 9.
#
#  excess cases in a band = ATT excess rate x treated person-years / 1e4
#  attributable fraction  = excess cases / observed cases in the band
# =====================================================================

source(file.path(Sys.getenv("CBCF_ROOT", "."), "Applied_Study_v2/code/common.R"))

cells <- read_csv(file.path(TAB, "T16_countbcf_cate_cells.csv.gz"), show_col_types = FALSE) %>%
  mutate(endpoint = sub(":.*", "", label), treated = sub(".*:", "", label))
ate <- read_csv(file.path(TAB, "T12_countbcf_ate.csv"), show_col_types = FALSE)

trt <- cells %>% filter(z == 1) %>% group_by(endpoint, treated) %>%
  summarise(treated_pyr = sum(pyr), treated_cases = sum(outcome), .groups = "drop")

# Preston Tables 9 (all solid) and 15 (stomach), aggregated to the four-level
# band definition used here.  Preston's stomach bands are defined on stomach
# dose; the case counts nevertheless line up closely with the colon-dose bands.
preston_t9 <- tribble(
  ~endpoint, ~treated, ~p_cases, ~p_excess, ~p_af,
  "solid",   "low",     4406,      81,      0.018,
  "solid",   "medium",  2112,     254,      0.120,
  "solid",   "high",    1148,     402,      0.350,
  "stomach", "low",     1226,      14,      0.012,
  "stomach", "medium",   559,      46,      0.082,
  "stomach", "high",     266,      69,      0.259)

af <- ate %>% dplyr::select(endpoint, treated, att_ear, att_ear_lo, att_ear_hi) %>%
  left_join(trt, by = c("endpoint", "treated")) %>%
  mutate(excess_cases = att_ear * treated_pyr / 1e4,
         excess_lo = att_ear_lo * treated_pyr / 1e4,
         excess_hi = att_ear_hi * treated_pyr / 1e4,
         attrib_frac = excess_cases / treated_cases,
         af_lo = excess_lo / treated_cases, af_hi = excess_hi / treated_cases) %>%
  left_join(preston_t9, by = c("endpoint", "treated")) %>%
  arrange(endpoint, factor(treated, levels = c("low", "medium", "high")))
write_csv(af, file.path(TAB, "T30_excess_cases_attributable_fraction.csv"))
print(as.data.frame(af), digits = 4)

tot <- af %>% group_by(endpoint) %>%
  summarise(excess_cases = sum(excess_cases), lo = sum(excess_lo), hi = sum(excess_hi),
            cases = sum(treated_cases), af = excess_cases / cases,
            preston_excess = sum(p_excess), preston_cases = sum(p_cases),
            .groups = "drop")
write_csv(tot, file.path(TAB, "T31_total_excess.csv"))
print(as.data.frame(tot), digits = 4)
stamp("06 done -> ", TAB)
