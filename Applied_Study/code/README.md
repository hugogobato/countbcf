# CountBCF applied-study analysis

The analysis uses the public RERF Life Span Study incidence table
(`Applied_Study/lssinc07.csv`, 26,806 cells, 1958–1998). Run commands from the
repository root. Required R packages include the dependencies in `DESCRIPTION`
and dplyr, readr, tibble, tidyr, splines, MASS, ggplot2, posterior, ranger, pROC,
sandwich, lmtest, scales, coda, and boot. The scripts report missing optional dependencies.

## Running within the memory budget

Run the pipeline, which installs the captured local package source in an
isolated temporary library before analysis:

```sh
CBCF_OUTDIR="$PWD/Applied_Study/analysis/clean_revision_run" CBCF_CORES=1 bash Applied_Study/code/run_all.sh
```

The default is one worker, with native-library threads limited to one. On Linux,
`run_all.sh` applies a 6 GiB address-space limit to each script. Leave at least
8 GiB for the desktop and other work. Large fits can take hours in total;
the former six-worker timing and memory estimates should not be used.
`CBCF_OUTDIR` can select a fresh output directory, and `CBCF_R_LIBRARY` can
select the isolated library. Set `CBCF_SKIP_PACKAGE_INSTALL=1` only when that
library already contains the intended source build. The runner preserves
normal dependency search paths and records the loaded package location. C++14 is required to
compile with current RcppArmadillo.

The age-30/70 comparison evaluates only eight prediction rows during sampling
(four target profiles and four training-row checks). It does not store full
forest histories. Its default profile budget is 5,000 burn-in iterations and
1,000 saved draws spaced ten iterations apart per chain. This budget retains
compact outputs, but the revision run still has between-chain disagreement
at the fixed profiles (pooled-sex R-hat 1.060 and 1.144, bulk ESS 24 and 9).
The table therefore contains provisional MCMC means and draw quantiles,
not validated posterior credible intervals. For a fresh process per chain, use:

```sh
CBCF_PROFILE_CHAIN=1 CBCF_CORES=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 prlimit --as=6442450944 -- Rscript Applied_Study/code/10_profile_comparison.R
CBCF_PROFILE_CHAIN=2 CBCF_CORES=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 prlimit --as=6442450944 -- Rscript Applied_Study/code/10_profile_comparison.R
CBCF_PROFILE_SUMMARISE_ONLY=1 Rscript Applied_Study/code/10_profile_comparison.R
```

Checkpoints contain source, installed-binary, data, seed, and settings identities;
reuse is refused when these differ. A successful checkpoint is written atomically.
Use a fresh `CBCF_OUTDIR` for a clean rerun after a source change; do not
discard the archived revision draws. Both profile summary scripts warn when
R-hat exceeds 1.01 or bulk ESS is below 400. Passing these checks alone would
not establish convergence or interval calibration.

## Scripts and outputs

| Script | Purpose |
|---|---|
| `common.R` | Data construction, propensity specifications, fitting and posterior summaries. |
| `01_eda_parametric.R` | Cohort validation, descriptives, Poisson/NB comparisons, Preston ERR/EAR refits. |
| `02_countbcf_main.R` | Two endpoints, three dose contrasts, two chains, average and subgroup summaries. |
| `03_sensitivity.R` | Prior and data-construction sensitivity analyses. |
| `04_overlap_cluster.R` | Propensity/overlap checks, trimming, a cohort-code mean-model sensitivity, and a parametric cluster bootstrap. The cohort code does not provide cluster-adjusted CountBCF intervals. |
| `06_attributable.R` | Excess cases and attributable fractions. |
| `07_lowdose_city.R` | Low-dose city contrast and sensitivity analyses. |
| `08_city_footing.R` | City posterior diagnostics, quantile MCSEs, reference-group/distance sensitivities, and E-values only for valid rate-ratio estimands. `CBCF_CITY_DIAGNOSTICS_ONLY=1` runs only the diagnostics. |
| `09_stomach_sweep.R` | Stomach prior-sensitivity and surface-flatness checks. |
| `10_profile_comparison.R` | Posterior evaluation at exposure age 30 and attained age 70, separately by sex and city, with a defined city mixture and diagnostics. |
| `11_rebuild_posterior_tables.R` | Recomputes the final city diagnostics/E-values and fixed-profile tables from archived compact draws, without fitting. |
| `05_figures.R` | Rebuilds figures from final tables; run last. |
| `export_reproducibility.py` | Exports only whitelisted final artifacts and executable source to `reproducibility/`, with hashes. |

Tables, plots, logs and model objects are written under `Applied_Study/analysis`
unless `CBCF_OUTDIR` is set. The public export excludes manuscripts, internal
notes, correspondence and preliminary outputs, without changing `.gitignore`.

## Estimands and interpretation

With person-years as the offset, `exp(mu_f)` is the reference rate per
person-year. Cell ERR is `exp(tau_f)-1`; EAR is
`(exp(mu_f+tau_f)-exp(mu_f))*10000`. The `*_serr_*` columns are the ratio of
person-year-weighted counterfactual rates minus one. This ratio differs from
the weighted mean of cell ERRs. Cohort summaries average over both exposure
groups; exposed-group summaries average only over exposed cells.

These are binary dose-group contrasts. Dividing an ERR by the observed mean
dose contrast gives a descriptive per-Gy rescaling, not a continuous-dose
slope. The age-30/70 table preserves this distinction and identifies CountBCF
draw quantiles versus Preston 90% confidence intervals. Inspect R-hat/ESS
before interpreting the former as credible intervals. A difference
of two ERRs is not a rate ratio; E-values and per-Gy slopes are left undefined
for that quantity. Quantile uncertainty is estimated with `posterior::mcse_quantile`,
not the MCSE of a posterior mean.

Unknown doses, NIC survivors, and doses at least 2 Gy are excluded from the
primary contrasts. Dedicated sensitivities change these restrictions. Distance
is strongly related to dose assignment, so distance-adjusted comparisons can
have severe overlap limitations. The cohort identifier is an integer-valued
forest covariate; it changes the mean model without modeling within-cohort
sampling dependence. The package's existing random-effect arguments do not
supply a cluster-adjusted count likelihood.

## Provenance

Historical tables come from several runs and lack complete execution provenance.
The revision adds fresh city diagnostics and profile posterior draws; it does
not retroactively establish the environments of the historical runs. The
runner captures executable source and data hashes before execution, records
package versions, and exports final summaries. See `reproducibility/README.md`
for archive reconstruction and the remaining final tagged-run release step.
