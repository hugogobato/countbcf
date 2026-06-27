# Simulation Study 1 — Count-Only DGPs (Poisson + Negative Binomial)

This document specifies the first of two simulation studies that will support
two separate papers built on `countbcf`.  It focuses on the **count component**
of the model — that is, settings where the data-generating process (DGP) is a
plain Poisson or Negative Binomial regression with **no structural-zero
component**.

The companion study (`sim_zi.md`) is the analogous experiment for the
**zero-inflated** part of the same model.

The benchmark is run as a set of **self-contained per-cell runners** under
`tests/runners_count/`.  Each runner is a single `.R` file that depends only on
`library(countbart)` — no `source()` of any repo file — so it can be uploaded
verbatim to a Google Colab session and executed.  With 15 CPUs available, 15
runners execute in parallel; the full grid (56 files) takes ~4 waves.

The legacy library file `sim_count.R` is retained as a reference implementation
of the DGPs and method wrappers but is **not** used by the per-cell runners.

---

## 1. Estimand

For unit $i$ with covariates $X_i \in \mathbb{R}^p$, treatment
$Z_i \in \{0,1\}$ and count outcome $Y_i \in \mathbb{Z}_{\ge 0}$, the potential
outcome means are
$$
\mu_0(x) = \mathbb{E}\!\left[ Y_i(0) \mid X_i = x \right], \qquad
\mu_1(x) = \mathbb{E}\!\left[ Y_i(1) \mid X_i = x \right].
$$
The targets are the conditional and average treatment effects on the
**response scale**:
$$
\tau(x) = \mu_1(x) - \mu_0(x), \qquad
\text{ATE} = \mathbb{E}_{X}\!\left[\tau(X)\right].
$$
All models are evaluated against these response-scale truths.

---

## 2. Data-generating processes

**Four** CATE-focused DGPs in a 2 × 2 design over **functional form**
(linear vs. nonlinear) and **count distribution** (Poisson vs. Negative
Binomial).  All share the same covariate distribution and propensity model so
cross-DGP comparisons are clean.

> **About CATE.** Each DGP below has strong, sign-changing, multi-covariate
> heterogeneous treatment effects — $\tau_h(X)$ is a nontrivial function of $X$
> so the per-unit CATE varies across the covariate space and per-unit
> estimation is a meaningful task.  The previously-defined ATE-focused
> variants (mild heterogeneity) have been dropped from this study; only the
> `*_cate` variants are run.

### 2.1 Covariates and treatment

- $X_i \in \mathbb{R}^p$ with the first $p-2$ columns drawn iid $\mathcal{N}(0,1)$;
- $X_{i, p-1} \sim \text{Bern}(0.6)$, $X_{i, p} \sim \text{Bern}(0.4)$ (so each
  design has a couple of binary columns);
- Propensity $\pi(x) = \sigma(0.5 x_1 - 0.4 x_2)$ where $\sigma$ is the logistic
  CDF; $Z_i \sim \text{Bern}(\pi(X_i))$. The estimated $\hat{\pi}$ passed to
  the models is set to $\pi(X_i)$ so that comparisons isolate the outcome
  model.

Only the first **five** covariates carry any signal; columns $X_{\cdot,6},\dots,X_{\cdot,p-2}$
are irrelevant noise (used to stress-test variable selection when $p$ grows).

### 2.2 Mean structure

For all four DGPs:
$$
\log \lambda_i \;=\; \mu_c(X_i) \;+\; Z_i \, \bigl[\,\tau_0 + \tau_h(X_i)\bigr].
$$

| DGP family | $\mu_c(x)$ | $\tau_h(x)$ — CATE-focused |
|------------|------------|----------------------------|
| **linear**    | $1.0 + 0.5 x_1 - 0.3 x_2 + 0.2 x_3$ | $0.60 x_1 - 0.45 x_2 + 0.30 x_3$ |
| **nonlinear** | $0.8 + 0.6\sin(x_1) + 0.4 x_2 x_3 + 0.5 \mathbf{1}\{x_4 > 0\}$ | $0.80 \sin(x_1) + 0.50 \mathbf{1}\{x_2 > 0\} x_3 - 0.40 x_5$ |

So the full set of DGPs is:

| Name | Functional form | Distribution |
|------|------------------|--------------|
| `linear_poisson_cate` | linear | Poisson |
| `nonlinear_poisson_cate` | nonlinear | Poisson |
| `linear_nb_cate` | linear | NB($\kappa = 2$) |
| `nonlinear_nb_cate` | nonlinear | NB($\kappa = 2$) |

The scalar $\tau_0$ is **calibrated** for each cell so that the realized ATE
matches the target value declared in the factor grid (see §3).  Concretely
`sim_count.R` solves
$$
\mathbb{E}_X\!\bigl[ e^{\mu_c(X)} (e^{\tau_0 + \tau_h(X)} - 1) \bigr] \;=\; \text{ATE}_{\text{target}}
$$
with `uniroot()` on a fresh held-out sample of 5,000 covariate vectors.

### 2.3 Sampling distribution

| DGP | $Y_i \mid X_i, Z_i$ |
|-----|---------------------|
| **linear_poisson**, **nonlinear_poisson** | $\text{Poisson}(\lambda_i)$ |
| **linear_nb**, **nonlinear_nb** | $\text{NB}(\text{size} = \kappa, \text{mean} = \lambda_i)$, $\kappa = 2$ |

The NB parameterization is $\text{Var}(Y) = \lambda + \lambda^2/\kappa$.  We
fix $\kappa = 2$ so the marginal variance is non-trivial without driving the
likelihood to a degenerate regime.

---

## 3. Experimental factors

The study sweeps three factors **one-at-a-time** (OAT) around a reference
configuration so that each plot has a single moving part.  Sweeps are *not*
combinatorial: each ablation cell varies exactly one factor relative to the
reference.  Reference values are $N_{\text{ref}} = 250$, $p_{\text{ref}} = 5$,
$\text{ATE}_{\text{ref}} = 1.25$.

| Factor | Levels (off-reference) | Held fixed |
|--------|-----------------------------------|----------------------------------|
| **Sample size** $N$  | 100, 500                | $p=5$,  ATE = 1.25 |
| **Covariate count** $p$ | 50, 250              | $N=250$, ATE = 1.25 |
| **Target ATE**       | 0.5, 2.5                | $N=250$, $p=5$ |

Together with the reference cell this gives **7 cells per (method, DGP)**:
one reference + 2 N ablations + 2 P ablations + 2 ATE ablations.

For each cell, all four CATE-focused DGPs are run.  Within a cell × DGP,
**`N_SIM = 100` Monte Carlo replicates** are drawn and fit.  The seed of
replicate $s$ is exactly $s$ — `set.seed(sim_id)`.  This makes individual
replicates trivially reproducible from their `(dgp, N, p, ATE, sim_id)` tuple.

`N_SIM` is a top-of-file constant in each runner.

---

## 4. Methods compared

Both methods produce a posterior over the ATE on the response scale and a
posterior mean of the unit-level outcome and CATE.

1. **CountBCF** — the proposed method.  Fit via `countbcf()` with
   `count_model = "poisson"` on Poisson DGPs and `count_model = "nb"` on NB
   DGPs.  Counterfactual means recovered from the raw forest coefficients:
   $\mu_0 = \exp(\mu_f)$, $\mu_1 = \exp(\mu_f + \tau_f)$, $\widehat{\text{CATE}} = \mu_1 - \mu_0$.
2. **BCF (Gaussian)** — `bcf_binary()` applied to the count outcome as if it
   were Gaussian.  The treatment-effect posterior is read off the moderate
   forest with `get_forest_fit()`.

Both methods use the same MCMC budget: `nburn = nsim = 1000`, `nthin = 1`.

> **CountBART (S-learner) — considered but dropped.** An S-learner built on
> `count_bart()` (augmenting the design matrix with $Z$ and $\hat\pi$ and
> evaluating the saved tree posterior at $\tilde x_i = (x_i, z, \hat\pi_i)$)
> was originally part of the benchmark.  In the current implementation it
> consistently exhausts available RAM and crashes the run on the larger
> ablation cells (notably $p \in \{50, 250\}$), so it has been removed from
> the comparison.  Until the memory footprint of `count_bart` is reduced,
> it is not feasible to include it here.

---

## 5. Metrics

Reported per (cell, dgp, sim, method) row of the pair-specific CSV.  Split
into ATE-level metrics and CATE-level metrics:

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
| `sim_count__<method>__<dgp>__N<N>_P<P>_ATE<ATE×100>.csv` | Long format, one row per sim. Streamed incrementally. |
| `sim_count__<method>__<dgp>__N<N>_P<P>_ATE<ATE×100>.log` | Per-cell timestamped progress log. |

With 2 methods × 4 CATE DGPs × 7 cells = **56 cell-CSVs per study**.  Each
CSV holds 100 rows (one per Monte Carlo replicate).  Concatenate them with
`rbind` in R or `cat *.csv` in shell for downstream plotting.

---

## 7. Running the study

### 7.1 On Google Colab (recommended)

Each runner in `tests/runners_count/` is **self-contained**: it depends only
on `library(countbart)` (and base R).  No `source()` of any repo file.  Upload
one runner per Colab session and `Rscript run_count_<method>_<dgp>_N<N>_P<P>_ATE<ATE×100>.R`.
With 15 Colab CPUs the 56 files run in roughly 4 waves of 15 (plus a partial).

### 7.2 Locally with GNU parallel

```bash
ls tests/runners_count/*.R | parallel -j 15 Rscript {}
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

Each `tests/runners_count/run_count_<method>_<dgp>_N<N>_P<P>_ATE<ATE×100>.R`
file inlines:

1. **Config block** — `N_SIM`, MCMC budget, the cell `(N, P, ATE)`, the DGP
   metadata (`KAPPA`, `COUNT_MODEL`).
2. **DGP** — the specific `.f_mu_count` and `.f_tau_het` functions for this
   DGP, plus `.simulate_count()`, `.calibrate_tau_const()`, and `draw_dgp()`.
3. **Method wrapper** — `fit_method()` specialized to the chosen method.
4. **Metric helpers** — `summarize_fit()`, `empty_metrics()`.
5. **Driver** — `run_one_fit()` and the top-level 100-sim loop that streams
   one row per replicate to the output CSV.

### 8.2 Legacy library (reference only)

`tests/sim_count.R` is kept for reference: it defines the same DGPs and method
wrappers as a library, plus a `run_pair()` driver that runs the full cell
sweep for one `(method, DGP)` pair.  It is **not** used by the per-cell
runners and is no longer the primary entry point.
