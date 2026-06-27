# countbcf

**Bayesian Causal Forests for Count and Zero-Inflated Count Outcomes.**

This repository contains two things:

1. **The `countbcf` R / C++ package** — a self-contained R package that
   implements Bayesian Causal Forests for count and zero-inflated count
   outcomes. Part I below documents installation, the model, usage, and
   the package internals.
2. **The replication code for the accompanying paper** — the full
   simulation pipeline, runners, analysis scripts, and outputs behind
   the working paper:

   > Souto, H. G., Wang, E. & Shuryak, I. (2026). *CountBCF: Bayesian Causal Forests for Count
   > and Zero-Inflated Count Outcomes.*

   Part II below documents the layout, the experimental design, and
   how to re-run the simulation studies end-to-end.

The package extends the Bayesian Causal Forest (BCF) of
Hahn, Murray and Carvalho (2020) to count and zero-inflated count
likelihoods. The MCMC backend, the GIG leaf prior, and the
zero-inflation bookkeeping are reused (and re-cast inside the
`(mu, tau)` BCF decomposition) from the upstream **`countbart`**
package of **Nathan B. Wikle and Corwin M. Zigler**
(<https://github.com/nbwikle/estimating-interference>), which in turn
builds on the log-linear BART model of Murray (2021). All credit for
that backend belongs to those authors; this package adds a causal
decomposition on top of it.

---

# Part I — the `countbcf` R package

## Installation

The package contains a sizeable C++ backend (Rcpp + RcppArmadillo,
OpenMP, Cereal). The most reliable way to install it on Google Colab
or any clean R session is via `devtools::install_github`:

```r
install.packages("remotes")
if (!require("devtools")) {
  install.packages("devtools")
}
devtools::install_github("hugogobato/countbcf")
library(countbcf)
```

The build pulls in `Rcpp`, `RcppArmadillo`, `Rcereal`, `GIGrvg`,
`fastDummies`, and `methods` automatically.

On Linux/macOS no extra setup is required. On Windows users need
Rtools matching their R version (`Rtools43` or newer is recommended).

### Google Colab quick start

```r
install.packages("remotes")
if (!require("devtools")) {
  install.packages("devtools")
}
devtools::install_github("hugogobato/countbcf")
library(countbcf)
```

---

## Model

For a count outcome `Y_i`, binary treatment `Z_i ∈ {0, 1}`, and
covariates `X_i` with an optional propensity score estimate
`pi_hat_i = E[Z | X_i]`, `countbcf` fits

**Non-zero-inflated count models** (`"poisson"`, `"nb"`):

```
log E[Y_i | X_i, Z_i] = mu_f(X_i, pi_hat_i) + Z_i * tau_f(X_i)
```

**Zero-inflated count models** (`"zipoisson"`, `"zinb"`): adds two
log-odds components for the structural-zero indicator `S_i`:

```
log lambda(X_i, Z_i)  =  mu_f (X_i, pi_hat_i) + Z_i * tau_f (X_i)
zi_logit(X_i, Z_i)    = (mu_f0(X_i, pi_hat_i) + Z_i * tau_f0(X_i))
                       - (mu_f1(X_i, pi_hat_i) + Z_i * tau_f1(X_i))
P(S_i = 1 | X_i, Z_i) =  sigmoid(zi_logit(X_i, Z_i))
Y_i | S_i = 1         =  0
Y_i | S_i = 0         ~  Poisson(lambda(X_i, Z_i))    [or NegBin(lambda, kappa)]
```

| Function group  | Prognostic forest          | Moderating forest           |
|-----------------|----------------------------|-----------------------------|
| Count (`f`)     | `mu_f(X, pi_hat)`          | `tau_f(X) * Z`              |
| Structural-zero (`f0`) | `mu_f0(X, pi_hat)`  | `tau_f0(X) * Z`             |
| Not structural-zero (`f1`) | `mu_f1(X, pi_hat)` | `tau_f1(X) * Z`            |

Totals: **2 forests** for non-ZI models (Poisson, NB) and **6 forests**
for ZI models (ZIP, ZINB). Each prognostic (mu) forest is a "vanilla"
BART; each moderating (tau) forest is entered linearly in `Z`. The leaf
prior is the generalized inverse Gaussian mixture of Murray (2021),
with smaller concentration on the tau forests (`a0/2`, `z_conc/2` by
default) to regularize heterogeneous treatment effects more strongly
than the prognostic component.

---

## Quick usage

A minimal end-to-end example with a Zero-Inflated Poisson likelihood:

```r
library(countbcf)

set.seed(1)
n  <- 1000
p  <- 5
X  <- matrix(rnorm(n * p), n, p)
pi <- plogis(0.5 * X[, 1] - 0.4 * X[, 2])
Z  <- rbinom(n, 1, pi)

log_lambda <- 1 + 0.5 * X[, 1] - 0.3 * X[, 2] + Z * (0.30 + 0.20 * X[, 1])
p_zi       <- plogis(-1 + 0.5 * X[, 2]    + Z * (-0.30 + 0.10 * X[, 3]))
Y          <- ifelse(rbinom(n, 1, p_zi) == 1, 0L, rpois(n, exp(log_lambda)))

fit <- countbcf(
  y           = Y,
  z           = Z,
  x_control   = X,
  x_moderate  = X,
  x_zero      = X,
  x_pos       = X,
  pihat       = pi,
  nburn       = 500,
  nsim        = 500,
  count_model = "zipoisson"
)
```

### Recovering CATE and ATE

`countbcf` returns the per-iteration raw forest coefficients in the
original input order. Combine them to recover potential outcomes and
treatment effects on the response scale:

```r
sigmoid <- function(z) 1 / (1 + exp(-z))

log_lambda_0 <- fit$mu_f_post
log_lambda_1 <- fit$mu_f_post + fit$tau_f_post
zi_logit_0   <- fit$mu_f0_post                   - fit$mu_f1_post
zi_logit_1   <- (fit$mu_f0_post + fit$tau_f0_post) -
                (fit$mu_f1_post + fit$tau_f1_post)

mu0          <- (1 - sigmoid(zi_logit_0)) * exp(log_lambda_0)   # nsim x n
mu1          <- (1 - sigmoid(zi_logit_1)) * exp(log_lambda_1)
cate_post    <- mu1 - mu0
ate_per_iter <- rowMeans(cate_post)

cate_hat <- colMeans(cate_post)                  # length n
ate_hat  <- mean(ate_per_iter)                   # posterior mean
ate_ci   <- quantile(ate_per_iter, c(0.025, 0.975))
```

For non-ZI models, only `mu_f_post` and `tau_f_post` are returned and
the CATE simplifies to
`exp(mu_f_post + tau_f_post) - exp(mu_f_post)`.

---

## Function signature

```r
countbcf(
  y, z, x_control,
  x_moderate = x_control, x_zero = x_control, x_pos = x_control,
  pihat = rep(0.5, length(y)),
  offset = NULL,
  nburn, nsim, nthin = 1, update_interval = 100,
  ntree_control  = 250, ntree_moderate  = 50,
  nztree_control = 100, nztree_moderate = 50,
  a0 = NA, a0_tau = NA,
  z_conc = 3.5 / sqrt(2), z_conc_tau = NA,
  base_control  = 0.95, power_control  = 2,
  base_moderate = 0.25, power_moderate = 3,
  kappa_a = 5, kappa_b = 3, kappa_prop_sd = 0.21,
  count_model   = "poisson",     # "poisson" | "nb" | "zipoisson" | "zinb"
  include_pihat = "control",     # "control" | "moderate" | "both" | "none"
  randeff_design = matrix(1),
  randeff_variance_component_design = matrix(1),
  randeff_scales = 1, randeff_df = 3,
  return_trees = FALSE,
  debug = FALSE
)
```

### Returned object (`countbcf_fit`)

| Element                                       | Shape          | Meaning                                                                                |
|-----------------------------------------------|----------------|----------------------------------------------------------------------------------------|
| `yhat_log`, `yhat`                            | `n × nsim`     | posterior of `log E[Y_i | X_i, Z_i = z_i^obs]` and its exponential                     |
| `order_vec`                                   | `length n`     | permutation used to sort `y` (zeros first); rows of `yhat*` follow this order          |
| `mu_f_post`, `tau_f_post`                     | `nsim × n`     | per-iter prognostic and moderating coefficients of the count component (original order)|
| `mu_f0_post`, `tau_f0_post`, `mu_f1_post`, `tau_f1_post` | `nsim × n` | per-iter coefficients of the ZI log-odds components (ZI models only)            |
| `mu_f_log`, `tau_f_log`, `mu_f0_log`, …       | `n × nsim`     | sorted-order in-sample fits (kept for backwards-compatibility with `countbart`)        |
| `kappa`, `kappa_acpt`                         | `length nsim`  | NB dispersion posterior and M-H acceptance rate (NB / ZINB only)                       |
| `control_fit$tree_samples`, …                 | object         | serialized `tree_samples` per forest (when `return_trees = TRUE`)                      |
| `random_effects`, `random_effects_sd`         | varies         | random-effects posteriors                                                              |
| `sigma`                                       | `length nsim`  | placeholder; unused in count models                                                    |

---

## Algorithm notes

The MCMC scheme follows `countbart` exactly; only the *forest update*
step is replaced to implement the BCF `(mu, tau)` decomposition. Within
the `bd → drmu → fit` block:

1. **mu forests** (vanilla, `omega = 1`): use the same sufficient
   statistics as `countbart`.
2. **tau forests** (non-vanilla, `omega = Z_trt`): because
   `omega == 0` for control units, those units contribute zero
   sufficient statistics to tau leaves — exactly as required by BCF.
   `(u_vec, r_tree)` are pre-multiplied by `omega` so that
   `allsuff_loglinear` does the right thing with `basis_dim == 1`.
3. **Per-group log f**: an auxiliary `log_f_per_group[g][k]` summing
   the in-sample fit of every forest in group `g` keeps the
   `r_tree` updates consistent across all six forests.
4. **kappa, latent variables, random effects**: identical to
   `countbart`.

Internally the data are sorted by `y` so that zeros come first
(required by the C++ kernel); `order_vec` is returned so you can map
sorted-order outputs (`yhat`, `*_log`) back to the original units. The
`*_post` matrices are already in the original order.

---

# Part II — Paper experiments and replication

This section documents the simulation pipeline that produces the
results reported in the paper
*CountBCF: Bayesian Causal Forests for Count and Zero-Inflated Count
Outcomes* (Souto, 2026). Everything in this section lives under the
`tests/` and `Papers/` folders of this repository and is independent
of installing the package as a user — it is intended for anyone who
wants to verify or extend the simulation study.

## Directory layout

```
tests/
├── _generate_runners.R       # emits one self-contained R script per (study, method, DGP, cell)
├── runners_count/            # ~56 emitted runners for the count-only sub-study
├── runners_zi/               # ~56 emitted runners for the zero-inflated sub-study
├── sim_count.R, sim_zi.R     # data-generating functions called by every runner
├── results_count/            # one .rds per runner: posterior samples + per-replicate truths
├── results_zi/
└── analysis/
    ├── analyze_results.R     # aggregates .rds files into per-cell metrics
    ├── summary_count.csv     # one row per cell, all metrics
    ├── summary_zi.csv
    ├── plots/                # PNGs used as figures in the paper
    └── tables/               # LaTeX tables used in the paper

Papers/
├── CountBCF_Merged.tex       # paper source
├── refs.bib                  # bibliography
└── CountBCF_Merged.pdf       # compiled paper (when present)
```

`tests/Experiment_Code_in_Notebooks/` contains Jupyter / Colab notebook
versions of the runners (one cell per `.R` runner) generated by
`tests/convert_to_notebooks.py`; use these on Google Colab to avoid the
local C++ toolchain. The paper itself includes figures from
`tests/analysis/plots/` via a relative `\graphicspath` directive, so
the document compiles directly against the analysis output without any
file copying.

## Factor grid and DGPs

`tests/_generate_runners.R` emits the one-at-a-time ablation grid
described in the paper: a reference cell at
`(N, p, ATE) = (250, 5, 1.25)` plus six ablation cells obtained by
moving exactly one factor off the reference (`N ∈ {100, 500}`,
`p ∈ {50, 250}`, `ATE ∈ {0.5, 2.5}`). The two sub-studies share a
single covariate distribution, propensity model, factor grid, and
Monte Carlo budget. The DGP factory is parameterised by four functions
per study:

| Object                | Role                                                    |
|-----------------------|---------------------------------------------------------|
| `F_MU_LINEAR / F_MU_NONLINEAR` | prognostic surface (log rate) for the count component |
| `F_TAU_LIN_CATE / F_TAU_NONLIN_CATE` | moderating surface (response-scale CATE) for the count component |
| `F_MU_ZI_LIN / F_MU_ZI_NONLIN` | prognostic surface (log-odds) for the structural-zero indicator (ZI study only) |
| `F_TAU_ZI_LIN_CATE / F_TAU_ZI_NL_CATE` | moderating surface (log-odds) for the structural-zero indicator (ZI study only) |

These functions are crossed with the four count likelihoods (Poisson,
NB, ZIP, ZINB) to give the eight DGPs reported in the paper. Each
runner draws `N_SIM = 100` independent replicates with `set.seed(s)`
for replicate `s`, so any single replicate is reproducible from its
`(dgp, N, p, ATE, sim_id)` tuple.

## Running the grid

```bash
Rscript tests/_generate_runners.R                # regenerate the runners
Rscript tests/runners_count/run_count_...R       # run one cell, ~5–15 min
Rscript tests/analysis/analyze_results.R         # aggregate all .rds → summary_*.csv + plots
```

Each runner writes a single `.rds` to `tests/results_count/` or
`tests/results_zi/` containing per-replicate point estimates,
posterior summaries, and the realised truths. `analyze_results.R`
walks the results folders, computes the per-cell metrics defined
below, and writes the summary CSVs and the PNG / LaTeX outputs used
in the paper.

## CSV column reference

`tests/analysis/summary_count.csv` and `summary_zi.csv` have one row
per `(method, DGP, N, p, ATE)` cell. The columns map onto the
mathematical notation used in the paper as follows.

| Column                  | Symbol in paper                  | Meaning |
|-------------------------|----------------------------------|---------|
| `n_sim`                 | `N_sim`                          | Monte Carlo replicates per cell (100). |
| `rmse_ate`              | RMSE(ATE)                        | RMSE of the posterior-mean ATE against the realised per-replicate ATE. |
| `ate_bias`              | bias(ATE)                        | Mean signed deviation of posterior-mean ATE from the realised ATE. |
| `ate_coverage_mean`     | `Cov_{95}(ATE)`                  | Empirical coverage of the 95% posterior interval for the ATE. |
| `ate_ci_width_mean`     | `W_{95}(ATE)`                    | Mean width of the 95% posterior interval for the ATE. |
| `pehe_mean`             | PEHE                             | Cross-replicate mean of the precision-in-estimating-heterogeneous-effects metric. |
| `pehe_sd`               | sd(PEHE)                         | Cross-replicate SD of PEHE (used for the error bars in the paper figures). |
| `cate_cov95_mean`       | `Cov_{95}(τ̂)`                   | Mean per-unit 95% credible-interval coverage for the CATE. |
| `cate_ci_width95_mean`  | `W_{95}(τ̂)`                     | Mean per-unit 95% credible-interval width for the CATE. |
| `elapsed_sec_mean`      | —                                | Mean wall-clock per replicate. |

Columns suffixed `_sd` are the cross-replicate SDs of the corresponding
`_mean` quantities. A small number of additional diagnostic columns
(`rmse_yhat`, `bias_yhat`, `mae_cate`, `cate_bias`, `cate_cor`,
`cate_cov50_mean`, `cate_sd_true`, `pct_zero`, `pct_struct_zero`) are
present in the CSV for diagnostic completeness but are not used in the
paper write-up; consult `tests/analysis/analyze_results.R` for the
definitions.

## Plot pipeline

The PNGs in `tests/analysis/plots/` are the ones included in the paper.
File names follow the pattern
`<study>_<sweep>_<metric>_<aggregator>.png`, e.g.
`count_N_cate_cov95_mean.png` is the per-unit CATE 95% coverage across
the `N` sweep on the count-only sub-study; `_normalized` variants are
per-DGP normalised to CountBCF's value at the smallest factor level
(used for the target-ATE sweep, where the absolute scale of PEHE and
RMSE(ATE) grows with the target ATE for trivial dimensional reasons).
The LaTeX `\graphicspath` directive in the paper points to
`../tests/analysis/plots/` so the figures resolve directly from the
analysis output without copying.

---

## Citation

If you use `countbcf` or its companion paper in academic work, please
cite the paper alongside the underlying methodology, including the
upstream `countbart` package:

```bibtex
@unpublished{souto2026countbcf,
  author = {Souto, Hugo Gobato},
  title  = {{CountBCF}: {B}ayesian Causal Forests for Count and
            Zero-Inflated Count Outcomes},
  year   = {2026}
}
```

---

## Acknowledgements

`countbcf` is a direct extension of the **`countbart`** R package by
**Nathan B. Wikle** and **Corwin M. Zigler**
(<https://github.com/nbwikle/estimating-interference>). The entire C++
sampling backend — log-linear BART, GIG leaf prior, zero-inflation
handling, the `tree_samples` serialization layer, and the
`bd/drmu/fit_loglinear` MCMC kernel — is theirs. CountBCF replaces only
the forest-update step to implement the BCF `(mu, tau)` decomposition.
Many thanks to those authors, and to Jared Murray for the log-linear
BART formulation that makes the count likelihood tractable in the
first place.

## Author

**Hugo Gobato Souto**
Dell Technologies
[hugo.souto@dell.com](mailto:hugo.souto@dell.com)
ORCID: <https://orcid.org/0000-0002-7039-0572>

## License

GPL (>= 3), matching the upstream `countbart` license.
