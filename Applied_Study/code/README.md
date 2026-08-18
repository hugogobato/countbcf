# CountBCF applied study, version 2 — code

Analysis of the **full** RERF Life Span Study solid-cancer incidence person-year
table (`../lssinc07.csv`, 26,806 cells, 1958–1998), replacing the partial file
used in `../../Applied_Study/`.

## Run it

```bash
cd <repo root>
CBCF_CORES=6 bash Applied_Study/code/run_all.sh
```

Roughly 40 minutes on a 20-thread machine with ~12 GB of free RAM.
`CBCF_CORES` caps the number of parallel MCMC workers; each worker holds two
`n_draws x n_cells` double matrices (≈150 MB for the largest contrast), so keep
`CBCF_CORES * 0.7 GB` comfortably below free memory. Outputs land in
`../analysis/tables` (CSV), `../analysis/plots` (PNG) and `../analysis/logs`.

## Files

| file | what it does |
|---|---|
| `common.R` | data construction (Brenner's `dcatn` recode, four-level dose bands, cluster/cohort ids), propensity-score specifications, the `countbcf` wrapper, and every posterior summary. Sourced by all other scripts. |
| `01_eda_parametric.R` | validates the data against Preston et al. (2007) Tables 2 and 4, descriptives and dispersion diagnostics, covariate balance, Poisson/NB GLM average effects, and a maximum-likelihood refit of Preston's own ERR and EAR risk models. |
| `02_countbcf_main.R` | the main CountBCF fits: 2 endpoints × 3 dose contrasts × 2 chains. ATE/ATT, convergence diagnostics, subgroup CATEs with the (CATE − ATE) separation test, and the comparison against the refitted Preston surface. |
| `03_sensitivity.R` | (A) regularisation of the moderating forest, (B) data perturbations (site composition, NIC, reference band, dose cap, likelihood, organ-dose banding), (C) low-dose sub-cohorts matching Dropkin (2016). |
| `04_overlap_cluster.R` | propensity-model comparison and overlap diagnostics, trimming sensitivity, subgroup overlap, CountBCF with random intercepts, and a cluster bootstrap of the parametric model. |
| `06_attributable.R` | excess cases and attributable fractions, placed next to Preston Tables 9 and 15. |
| `07_lowdose_city.R` | the robustness battery for the low-dose Hiroshima-vs-Nagasaki ERR difference: moderating-forest prior scale and depth prior, shared leaf concentration, four propensity specifications, a trim that actually bites, and the NIC group. Plus an explicit parametric dose × city interaction test. |
| `08_city_footing.R` | the follow-up on that difference. Chain-wise draws, R-hat, ESS and Monte-Carlo standard error for the standardised city statistic; E-values and per-Gy rescaling; the RERF distance-stratified reference constructions (outer proximal 2.5–3.0 km vs distal 3–10 km, city × distance in the score and in the prognostic forest); the low band split into 5–<20 / 20–<50 / 50–<100 mGy to separate a dose gradient from a step at the 5 mGy boundary; site-composition perturbations; and the parametric analogue under four reference-group constructions. |
| `09_stomach_sweep.R` | repeats the twelve-setting moderating-forest regularisation grid for the stomach high-dose contrast, recording the standardised and per-cell subgroup flags, surface flatness, and the age, sex, and city contrasts. |
| `05_figures.R` | all figures. Runs last because it reads the tables written by 01–04 and 06–09. |

## Estimands

Person-years are the offset, so `exp(mu_f)` is the cases-per-person-year rate in
the reference arm and `tau_f` is the log rate ratio for the exposed band. For
each person-year cell *i*:

- **EAR / absolute CATE** — `(exp(mu_f + tau_f) − exp(mu_f)) × 10^4`, excess
  cases per 10,000 person-years (Preston's unit);
- **ERR / relative CATE** — `exp(tau_f) − 1`;
- **standardised ERR** — `(Σ w r₁)/(Σ w r₀) − 1`, the ratio of person-year
  weighted rates. This, not the average of the per-cell ERRs, is what a
  parametric g-computation returns, so it is the quantity compared with the GLM
  rate ratio and with Preston's case-weighted summaries.

`ATE` averages over both arms of a contrast, `ATT` over the exposed arm only.
Every summary is computed on the posterior draws, so a subgroup contrast such as
`CATE_g − ATE` carries a genuine credible interval rather than a difference of
point estimates.

## Design decisions worth knowing

- Records with `dcat == 1` (unknown DS02 dose) are dropped; `dcat == 2` is split
  into not-in-city survivors and in-city survivors with 0–<5 mGy, exactly as in
  the Epicure log supplied by A. Brenner.
- The reference arm is the in-city 0–<5 mGy group. Not-in-city survivors are
  excluded from the primary contrasts — they had no positive probability of
  exposure, so including them would violate positivity — and enter a sensitivity
  fit as extra control-arm cells with an NIC indicator in the prognostic forest.
- Cells with weighted colon dose ≥ 2 Gy are excluded from the primary analysis
  and restored in a sensitivity fit.
- Ground distance from the hypocentre is deliberately kept out of the *primary*
  propensity model: it is the dose mechanism rather than a confounder, and
  conditioning on it destroys overlap. `08_city_footing.R` relaxes this in a
  dedicated family of fits, using the coarse `<3 km` / `>=3 km` ring, which is a
  baseline-composition variable rather than the continuous assignment variable.
  Note the structural fact behind that choice: the 5–<100 mGy band contains no
  cell beyond 2.71 km, so reference cells beyond 3.22 km have probability zero
  of low-band membership — the same positivity argument that excludes the NIC
  group, one ring closer in.
- `agex` and `age` determine calendar year, birth cohort and time since exposure
  (`year = age − agex + t0`, with the implied `t0` between 1945.594 and 1945.611
  across all cells, i.e. 6 and 9 August 1945 to within about six days), so only
  those two are given to the forests alongside sex and city.
- The **standardised ERR** is the footing every subgroup claim is made on.
  Scripts that score subgroups report both footings; the `*_serr_*` columns are
  the ones the manuscript uses.
- `countbcf`'s `randeff_*` arguments are **inert** — the random effects are
  sampled but never enter the count likelihood (`R/countbcf.R:109-113`). Group
  structure goes in `x_control` instead. Never label a `cluster_re = TRUE` fit a
  random-intercept model.
- `a0` is the leaf concentration **shared** by the prognostic and moderating
  forests (`R/countbcf.R:90`), with implied treatment-effect prior scale
  `a0 * sqrt(ntree_moderate / ntree_control)`. A moderating-forest-only sweep
  uses `ntree_moderate`, `base_moderate` or `power_moderate`.
