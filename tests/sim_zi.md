# Simulation Study 2 — Zero-Inflated DGPs (ZIP + ZINB)

This document specifies the second of the two simulation studies that will
support two separate papers built on `countbcf`.  It focuses on the
**zero-inflated component** of the model — settings where the DGP has both a
count regression part *and* an independent structural-zero process governed by
covariates.

The companion study (`sim_count.md`) is the analogous experiment for the
non-zero-inflated case (Poisson + NB only).

The benchmark is run as a set of **self-contained per-cell runners** under
`tests/runners_zi/`.  Each runner is a single `.R` file that depends only on
`library(countbart)` — no `source()` of any repo file — so it can be uploaded
verbatim to a Google Colab session and executed.  With 15 CPUs available, 15
runners execute in parallel; the full grid (56 files) takes ~4 waves.

The legacy library file `sim_zi.R` is retained as a reference implementation
of the DGPs and method wrappers but is **not** used by the per-cell runners.

---

## 1. Estimand

For unit $i$ with covariates $X_i \in \mathbb{R}^p$, treatment
$Z_i \in \{0,1\}$ and count outcome $Y_i \in \mathbb{Z}_{\ge 0}$, the potential
outcome means are
$$
\mu_z(x) = (1 - p_{\mathrm{zi}, z}(x)) \cdot \exp(\log\lambda_z(x)),
\quad z \in \{0,1\}.
$$
That is, $E[Y_i(z) \mid X_i = x]$ is the product of the probability the unit
is *not* a structural zero and the count-component mean.  The targets are the
response-scale CATE and ATE,
$$
\tau(x) = \mu_1(x) - \mu_0(x), \qquad \text{ATE} = \mathbb{E}_X[\tau(X)].
$$
All models are evaluated against these response-scale truths.

Because both the count rate **and** the structural-zero probability can move
with $Z$, the ATE here is a marginal effect that mixes two channels.

---

## 2. Data-generating processes

**Four** CATE-focused DGPs in a 2 × 2 design over **functional form**
(linear vs. nonlinear) and **count distribution** (Zero-Inflated Poisson vs.
Zero-Inflated Negative Binomial).  Covariates, propensity, and the qualitative
structural-zero mechanism are shared across DGPs so cross-DGP comparisons are
clean.

> **About CATE.** Each DGP below has strong, sign-changing, multi-covariate
> heterogeneous treatment effects on **both** the count rate and the
> structural-zero probability — both $\tau_h(X)$ and $\Delta\eta_{\mathrm{zi}}(X)$
> are nontrivial functions of $X$.  The previously-defined ATE-focused
> variants (mild single-covariate heterogeneity) have been dropped from this
> study; only the `*_cate` variants are run.

### 2.1 Covariates and treatment

Identical to the count-only study (`sim_count.md` §2.1):

- $X_i \in \mathbb{R}^p$, first $p-2$ columns iid $\mathcal{N}(0,1)$, last
  two columns Bernoulli($0.6$) and Bernoulli($0.4$);
- $\pi(x) = \sigma(0.5 x_1 - 0.4 x_2)$, $Z_i \sim \text{Bern}(\pi(X_i))$;
- $\hat\pi$ passed to the models equals the truth, isolating the outcome
  model.

Only $X_{\cdot, 1:5}$ carry signal; the rest are noise.

### 2.2 Count component

For all eight DGPs:
$$
\log \lambda_i \;=\; \mu_c(X_i) \;+\; Z_i \, \bigl[\,\tau_0 + \tau_h(X_i)\bigr].
$$

| DGP family | $\mu_c(x)$ | $\tau_h(x)$ — CATE-focused |
|------------|------------|----------------------------|
| **linear**    | $1.0 + 0.5 x_1 - 0.3 x_2 + 0.2 x_3$ | $0.60 x_1 - 0.45 x_2 + 0.30 x_3$ |
| **nonlinear** | $0.8 + 0.6 \sin(x_1) + 0.4 x_2 x_3 + 0.5 \mathbf{1}\{x_4 > 0\}$ | $0.80 \sin(x_1) + 0.50 \mathbf{1}\{x_2 > 0\} x_3 - 0.40 x_5$ |

### 2.3 Zero-inflation component

Baseline $P(\text{struct. zero}) \approx 25\text{–}30\%$.  Control logit
$\eta_{\mathrm{zi}}(x)$ is shared between ATE and CATE variants of the same
functional form; only the **treatment shift** $\Delta\eta_{\mathrm{zi}}(x)$
differs:

| DGP family | $\eta_{\mathrm{zi}}(x)$ (control logit) | $\Delta\eta_{\mathrm{zi}}(x)$ — CATE-focused |
|------------|-----------------------------------------|----------------------------------------------|
| **linear**    | $-1.0 + 0.5 x_2$                                  | $-0.50 x_1 + 0.30 x_3 - 0.40 x_4$ |
| **nonlinear** | $-1.0 + 0.7 \lvert x_2 \rvert - 0.5 \exp(-x_3^2)$ | $-0.60 x_1 x_4 + 0.40 x_2 \mathbf{1}\{x_5 > 0\}$ |

Hence $p_{\mathrm{zi}, z}(x) = \sigma(\eta_{\mathrm{zi}}(x) + z \Delta\eta_{\mathrm{zi}}(x))$.
The CATE-focused $\Delta\eta_{\mathrm{zi}}$ is sign-changing across units, so
for some units the treatment *increases* the structural-zero probability while
for others it decreases it — exactly the kind of heterogeneity a per-unit
method should be able to recover.

So the full set of DGPs is:

| Name | Functional form | Distribution |
|------|------------------|--------------|
| `linear_zip_cate` | linear | ZIP |
| `nonlinear_zip_cate` | nonlinear | ZIP |
| `linear_zinb_cate` | linear | ZINB($\kappa = 2$) |
| `nonlinear_zinb_cate` | nonlinear | ZINB($\kappa = 2$) |

### 2.4 Sampling distribution

For each unit $i$, draw a Bernoulli structural-zero indicator
$S_i \sim \text{Bern}(p_{\mathrm{zi}, Z_i}(X_i))$.  If $S_i = 1$ then
$Y_i = 0$; otherwise

| DGP | $Y_i \mid X_i, Z_i, S_i = 0$ |
|-----|------------------------------|
| **linear_zip**, **nonlinear_zip**     | $\text{Poisson}(\lambda_i)$ |
| **linear_zinb**, **nonlinear_zinb**   | $\text{NB}(\text{size} = \kappa, \text{mean} = \lambda_i)$, $\kappa = 2$ |

### 2.5 ATE calibration

The intercept $\tau_0$ on the count component is calibrated per cell so that
the realized response-scale ATE matches the target.  Concretely `sim_zi.R`
solves
$$
\mathbb{E}_X\!\Bigl[(1 - p_{\mathrm{zi},1}(X)) e^{\mu_c(X) + \tau_0 + \tau_h(X)}
                  - (1 - p_{\mathrm{zi},0}(X)) e^{\mu_c(X)}\Bigr]
\;=\; \text{ATE}_{\text{target}}
$$
with `uniroot()` on a held-out sample of 5,000 covariate vectors.  Only
$\tau_0$ is calibrated; $\Delta\eta_{\mathrm{zi}}$ (the ZI treatment effect) is
held fixed.

---

## 3. Experimental factors

Same one-at-a-time sweep as in the count-only study (sweeps are *not*
combinatorial — each ablation cell varies exactly one factor relative to the
reference).  Reference cell $(N, p, \text{ATE}) = (250, 5, 1.25)$:

| Factor | Levels (off-reference) | Held fixed |
|--------|--------|------------|
| **Sample size** $N$       | 100, 500                | $p=5$,   ATE = 1.25 |
| **Covariate count** $p$   | 50, 250                 | $N=250$, ATE = 1.25 |
| **Target ATE**            | 0.5, 2.5                | $N=250$, $p=5$ |

Together with the reference cell this gives **7 cells per (method, DGP)**.

For each cell, all four CATE-focused DGPs are run.  Within a cell × DGP,
**`N_SIM = 100`** Monte Carlo replicates are drawn and fit; replicate $s$
uses seed $s$ (`set.seed(sim_id)`).

---

## 4. Methods compared

Both methods produce a posterior over the response-scale ATE and the
unit-level outcome / CATE.

1. **CountBCF** — `countbcf()` with `count_model = "zipoisson"` on ZIP DGPs
   and `count_model = "zinb"` on ZINB DGPs.  Counterfactual response-scale
   means are computed from the six saved per-iteration coefficient matrices:
   $$
   \log\lambda_z = \mu_f + z \tau_f, \quad
   \eta_{\mathrm{zi}, z} = (\mu_{f_0} + z \tau_{f_0}) - (\mu_{f_1} + z \tau_{f_1}),
   $$
   $$
   \mu_z = (1 - \sigma(\eta_{\mathrm{zi}, z})) \exp(\log\lambda_z),
   \quad \widehat{\text{CATE}} = \mu_1 - \mu_0.
   $$
2. **BCF (Gaussian)** — `bcf_binary()` applied to the count outcome as if it
   were Gaussian.  Ignores zero-inflation entirely.  Tau posterior read off
   the moderate forest via `get_forest_fit()`.

Both methods use the same MCMC budget: `nburn = nsim = 1000`, `nthin = 1`.

> **CountBART (S-learner) — considered but dropped.** An S-learner built on
> `count_bart()` (augmenting the design matrix with $Z$ and $\hat\pi$, with
> `count_model = "zipoisson"` or `"zinb"` matching the DGP, and counterfactual
> predictions assembled from the three saved forests
> $\mu_z = (1 - \sigma(f_0(\tilde x) - f_1(\tilde x))) \exp(f(\tilde x))$) was
> originally part of the benchmark.  In the current implementation it
> consistently exhausts available RAM and crashes the run on the larger
> ablation cells (notably $p \in \{50, 250\}$ in the ZI study), so it has been
> removed from the comparison.  Until the memory footprint of `count_bart` is
> reduced, it is not feasible to include it here.

---

## 5. Metrics

Reported per (cell, dgp, sim, method) row of the pair-specific CSV.  Split
into outcome-fit, ATE-level, and CATE-level metrics:

**Outcome-fit (in-sample observed arm).**

| Metric | Definition |
|--------|------------|
| `rmse_yhat` | $\sqrt{\overline{(\hat y_i - \mu_{Z_i}(X_i))^2}}$ |
| `bias_yhat` | $\overline{\hat y_i - \mu_{Z_i}(X_i)}$ |

**ATE.**

| Metric | Definition |
|--------|------------|
| `ate_mean`, `ate_q025`, `ate_q975` | Posterior mean and 95 % CI for ATE |
| `ate_coverage` | 1 if the 95 % CI contains the realized ATE, 0 otherwise |
| `true_ate` | Realized ATE for that draw |

**CATE.**  All on the response scale.

| Metric | Definition |
|--------|------------|
| `pehe` | $\sqrt{\overline{(\widehat\tau(X_i) - \tau(X_i))^2}}$  — Precision in Estimation of Heterogeneous Effects |
| `mae_cate` | $\overline{\lvert \widehat\tau(X_i) - \tau(X_i) \rvert}$ |
| `cate_bias` | $\overline{\widehat\tau(X_i) - \tau(X_i)}$ |
| `cate_cor` | Pearson correlation between $\widehat\tau(X_i)$ and $\tau(X_i)$ — does the model rank units correctly? |
| `cate_cov95`, `cate_cov50` | Fraction of units whose 95 % / 50 % unit-level posterior CI contains $\tau(X_i)$ |
| `cate_ci_width95` | Mean width of the unit-level 95 % CIs (sharpness) |
| `cate_sd_true` | $\text{sd}(\tau(X_i))$ at the drawn dataset (DGP-level "how heterogeneous is the truth?") |

**Bookkeeping.**

| Metric | Definition |
|--------|------------|
| `pct_zero`, `pct_struct_zero` | Realized total-zero and structural-zero rates |
| `elapsed_sec` | Wallclock for the fit |
| `error` | `NA` on success; otherwise the caught error message |

Aggregating over the `N_SIM` replicates within a cell × DGP × method gives
the standard Monte Carlo summaries the paper will plot.  Pay particular
attention to `pehe` and `cate_cov95` on the `*_cate` DGPs — that is where
the methods are most differentiated.

---

## 6. Outputs

Each (method, DGP, cell) tuple writes its own CSV + log to the working
directory of the runner:

| File | Contents |
|------|----------|
| `sim_zi__<method>__<dgp>__N<N>_P<P>_ATE<ATE×100>.csv` | Long format, one row per sim. Streamed incrementally. |
| `sim_zi__<method>__<dgp>__N<N>_P<P>_ATE<ATE×100>.log` | Per-cell timestamped progress log. |

With 2 methods × 4 CATE DGPs × 7 cells = **56 cell-CSVs per study**.  Each
CSV holds 100 rows (one per Monte Carlo replicate).  Concatenate them with
`rbind` in R or `cat *.csv` in shell for downstream plotting.

---

## 7. Running the study

### 7.1 On Google Colab (recommended)

Each runner in `tests/runners_zi/` is **self-contained**: it depends only on
`library(countbart)` (and base R).  No `source()` of any repo file.  Upload
one runner per Colab session and `Rscript run_zi_<method>_<dgp>_N<N>_P<P>_ATE<ATE×100>.R`.
With 15 Colab CPUs the 56 files run in roughly 4 waves of 15 (plus a partial).

### 7.2 Locally with GNU parallel

```bash
ls tests/runners_zi/*.R | parallel -j 15 Rscript {}
```

Each file produces its own CSV + log in the working directory.

### 7.3 Regenerating the runners

The runners are emitted by `tests/_generate_runners.R`, which contains the
cell grid, DGP function bodies, and method wrappers as templated strings.
Edit it and re-run

```bash
Rscript tests/_generate_runners.R
```

to refresh both `runners_count/` and `runners_zi/` in lock-step.

---

## 8. Code listing

### 8.1 Per-cell runners (the actual benchmark)

Each `tests/runners_zi/run_zi_<method>_<dgp>_N<N>_P<P>_ATE<ATE×100>.R` file
inlines:

1. **Config block** — `N_SIM`, MCMC budget, the cell `(N, P, ATE)`, DGP
   metadata (`KAPPA`, `COUNT_MODEL`).
2. **DGP** — the DGP's `.f_mu_count`, `.f_tau_het`, `.f_mu_zi`, `.f_tau_zi`
   functions plus `.simulate_zi()`, `.calibrate_tau_const_zi()` (which
   accounts for the $(1 - p_{\mathrm{zi}})$ thinning of the count component),
   and `draw_dgp()`.
3. **Method wrapper** — `fit_method()` specialized to the chosen method.
4. **Metric helpers** — `summarize_fit()`, `empty_metrics()`.
5. **Driver** — `run_one_fit()` and the top-level 100-sim loop that streams
   one row per replicate to the output CSV.

### 8.2 Legacy library (reference only)

`tests/sim_zi.R` is kept for reference: it defines the same DGPs and method
wrappers as a library, plus a `run_pair()` driver that runs the full cell
sweep for one `(method, DGP)` pair.  It is **not** used by the per-cell
runners and is no longer the primary entry point.
