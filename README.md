# countbcf

**Bayesian Causal Forests for Count and Zero-Inflated Count Outcomes.**

`countbcf` is the R/C++ companion package to the working papers

> Souto, H. G. (2026a). *CountBCF: Bayesian Causal Forests for Count Outcomes.*
>
> Souto, H. G. (2026b). *Zero-Inflated CountBCF: Bayesian Causal Forests for Zero-Inflated Count Outcomes.*

It estimates conditional and average treatment effects for non-negative
integer responses by extending the Bayesian Causal Forest (BCF) of
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

## Installation

The package contains a sizeable C++ backend (Rcpp + RcppArmadillo,
OpenMP, Cereal). The most reliable way to install it on Google Colab
or any clean R session is via `devtools::install_github`:

```r
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
system("apt-get update -qq && apt-get install -y libcurl4-openssl-dev libssl-dev libxml2-dev")
install.packages(c("devtools", "Rcpp", "RcppArmadillo"))
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

## Citation

If you use `countbcf` in academic work, please cite **both** the
papers behind the package and the underlying methodology, including
the upstream `countbart` package:

```bibtex
@unpublished{souto2026countbcf,
  author = {Souto, Hugo Gobato},
  title  = {{CountBCF}: {B}ayesian Causal Forests for Count Outcomes},
  year   = {2026}
}

@unpublished{souto2026zicountbcf,
  author = {Souto, Hugo Gobato},
  title  = {{Zero-Inflated CountBCF}: {B}ayesian Causal Forests for
            Zero-Inflated Count Outcomes},
  year   = {2026}
}

@article{hahn2020bayesian,
  author  = {Hahn, P. Richard and Murray, Jared S. and Carvalho, Carlos M.},
  title   = {Bayesian Regression Tree Models for Causal Inference:
             Regularization, Confounding, and Heterogeneous Effects},
  journal = {Bayesian Analysis},
  volume  = {15},
  number  = {3},
  pages   = {965--1056},
  year    = {2020},
  doi     = {10.1214/19-BA1195}
}

@article{murray2021loglinear,
  author  = {Murray, Jared S.},
  title   = {Log-Linear {B}ayesian Additive Regression Trees for
             Multinomial Logistic and Count Regression Models},
  journal = {Journal of the American Statistical Association},
  volume  = {116},
  number  = {534},
  pages   = {756--769},
  year    = {2021},
  doi     = {10.1080/01621459.2020.1813587}
}

@article{wikle2023causal,
  author  = {Wikle, Nathan B. and Zigler, Corwin M.},
  title   = {Causal Health Impacts of Power Plant Emission Controls under
             Modeled and Uncertain Physical Process Interference},
  journal = {Annals of Applied Statistics},
  year    = {2023},
  note    = {Companion package \texttt{countbart}: \url{https://github.com/nbwikle/estimating-interference}}
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
