#' Fit log-linear BART for count regression models.
#'
#' @references Jared S. Murray (2021). Log-Linear Bayesian Additive Regression Trees for Multinomial Logistic and Count Regression Models.
#'  https://doi.org/10.1080/01621459.2020.1813587. 
#'
#' @details Fits a generalized version of the log-linear BART model (Murray 2021): For a response
#' variable y and covariates x,
#' \deqn{y_i ~ Pr(Y_i = y_i | \mu(x_i))}
#' where \eqn{Pr(Y_i = y_i | \mu(x_i))} is the PMF from a Poisson or negative binomial regression 
#' model (or their zero-inflated variants). Furthermore, \eqn{\mu(x_i) = \mu_{0i} f(x_i)}, for fixed
#' offset \eqn{\mu_{0i}} and BART prior over log-linear predictor function \eqn{log(f(x_i))}. 
#'
#' Some notes:
#' \itemize{
#'    \item x (and x_zero and x_pos) must be numeric matrices. See e.g. the makeModelMatrix function in the
#'    dbarts package for appropriately constructing a design matrix from a data.frame
#'    \item The prior for leaf parameter \eqn{\lambda} is a mixture of generalized inverse Gaussian (GIG)
#'    distributions, with hyperparameters \eqn{c} and \eqn{d}. By default, \eqn{c} and \eqn{d} are chosen
#'    so that \eqn{var(\lambda) = a^2_{0} / m}, where \eqn{a_0} is a concentration parameter and \eqn{m}
#'    denotes the number of trees in \eqn{f} (see Proposition 4.1 from Murray (2021) for details). For 
#'    the log-linear predictor function, \eqn{a_0} is chosen such that \eqn{Pr(f(x_i) \leq y^{*}) \approx 0.975}
#'    marginally, where \eqn{y^{*}} is the 0.95 quantile of the empirical distribution of y. For zero-inflated 
#'    odds components \eqn{f^{(0)}(x_i)} and \eqn{f^{(1)}(x_i)}, \eqn{a_0 = 3.5 / \sqrt{2}}.
#'    \item By default, \eqn{\kappa} (the negative binomial dispersion parameter) has a beta prime prior with 
#'    shape parameters \eqn{\alpha = 5} and \eqn{\beta = 3}.
#' }
#'
#' @param y Response variable (must be vector of non-negative integers).
#' @param x Design matrix for the log-linear predictor, f(x).
#' @param x_zero Design matrix for zero-tree, which gives the odds a zero is due to the zero-inflated point mass component. Defaults to x.
#' @param x_pos Design matrix for positive-tree, which gives the odds a zero is NOT due to the zero-inflated point mass component. Defaults to x.
#' @param offset Offset value; can be a scalar or vector. If 'NULL', an offset of 1 will be used.
#' @param nburn Number of burn-in MCMC iterations.
#' @param nsim Number of MCMC iterations to save after burn-in.
#' @param nthin Save every nthin'th MCMC iterate. The total number of MCMC iterations will be nsim*nthin + nburn.
#' @param update_interval Print status every update_interval MCMC iterations.
#' @param ntree Number of trees in log-linear predictor component.
#' @param a0 Concentration parameter. 'NA' gives default estimate from Murray (2022).
#' @param nztree Number of trees in zero-inflated point mass probability component.
#' @param z_conc Concentration parameter for zero-inflated log-odds prior.
#' @param sd_val SD(mu(x)) marginally at any covariate value (or its prior median if use_muscale=TRUE).
#' @param base Base for tree prior on trees (see details).
#' @param power Power for the tree prior on trees.
#' @param nu Degrees of freedom in the chisq prior on \eqn{sigma^2}.
#' @param lambda Scale parameter in the chisq prior on \eqn{sigma^2}.
#' @param sigq Calibration quantile for the chisq prior on \eqn{sigma^2}.
#' @param sighat Calibration estimate for the chisq prior on \eqn{sigma^2}.
#' @param kappa_a Shape parameter (alpha) for kappa prior (beta prime distribution); default is 5.
#' @param kappa_b Shape parameter (beta) for kappa prior (beta prime distribution); default is 3.
#' @param kappa_prop_sd Standard deviation of kappa's Metropolis-Hasting proposal distribution.
#' @param count_model Specifies the count model; takes values "poisson", "nb", "zipoisson", or "zinb".
#' @param debug XX
#' @param randeff_design XX
#' @param randeff_variance_component_design XX
#' @param randeff_scales XX
#' @param randeff_df XX
#'
#' @return A list with posterior estimates of \eqn{\log(E(Y_i | x_i))}, posterior samples of trees,
#' and posterior estimates of dispersion parameter \eqn{\kappa} (for negative binomial models).
#' 
#' @export
#' 
count_bart <- function(
  y,
  x,
  x_zero = x,
  x_pos = x,
  offset = NULL,
  nburn, 
  nsim, 
  nthin = 1, 
  update_interval = 100,                     
  ntree = 250,
  a0 = NA, # 
  nztree = 100,
  z_conc = 3.5 / sqrt(2),
  sd_val = 2 * sd(y),
  base = 0.95,
  power = 2,
  nu = 3, lambda = NULL, sigq = .9, sighat = NULL, # prior specification for continuous y... only used as 'ghost' inputs
  kappa_a = 5, kappa_b = 3, # shape parameters for kappa prior (beta prime distribution)
  kappa_prop_sd = 0.21,     # standard deviation of kappa M-H proposal distribution
  count_model = "poisson",  # must be either 'poisson', 'nb', 'zipoisson', or 'zinb'
  debug = FALSE,
  # used as 'ghost' inputs... but not relevant for count models...
  randeff_design = matrix(1),
  randeff_variance_component_design = matrix(1),
  randeff_scales = 1,
  randeff_df = 3,
  return_trees = TRUE
){
  
  ### check that inputs are valid

  # check if count_model is correctly specified
  if (!(count_model %in% c("poisson", "nb", "zipoisson", "zinb"))) {
    stop(
      "Undefined 'count_model'; must be one of: \n
         'poisson','nb' (negative binomial), 'zipoisson' (zero-inflated poisson), 'zinb' (zero-inflated negative binomial)"
    )
  } else {
    # convert to numeric value for 'countbart'
    model_type <- which(count_model == c("poisson", "nb", "zipoisson", "zinb"))
  }


  ### check that data structures are of correct size

  # determine if x is a list of design matrices or a single matrix
  if (is.vector(x)){
    # x is a list with multiple design matrices
    num_designs <- length(x)
  } else {
    num_designs <- 1
  }

  # check for correct inputs
  if (num_designs > 1) {
    
    # check that design matrices have the same row size
    xmat_rows <- unlist(lapply(x, nrow))
    if (length(unique(xmat_rows)) != 1){
      stop(
        "Data size mismatch. The row size should be the same across design matrices."
      )
    }

    # non-numeric values
    if (any(!is.finite(y))) stop("Non-numeric values in y")
    if (any(unlist(lapply(x, function(x_mat) {
      any(!is.finite(x_mat))
    })))) stop("Non-numeric values in x")


    # missing values
    if (any(is.na(y))) stop("Missing values in y")
    if (any(unlist(lapply(x, function(x_mat) {
      any(is.na(x_mat))
    })))) stop("Missing values in x")
    
    # check if input sizes are correct
    if ((model_type == 3) | (model_type == 4)) {
      # include x_zero and x_pos in data check
      if (
        !.ident(
          length(y), xmat_rows[1], nrow(x_zero), nrow(x_pos)
        )
      ) {
        stop(
          "Data size mismatch. The following should all be equal:
          length(y): ", length(y), "\n",
          "nrow(x): ", nrow(x), "\n",
          "nrow(x_zero): ", nrow(x_zero), "\n",
          "nrow(x_pos): ", nrow(x_pos), "\n"
        )
      }

      # check that x_zero and x_pos have numeric values
      if (any(!is.finite(x_zero))) stop("Non-numeric values in x_zero")
      if (any(!is.finite(x_pos))) stop("Non-numeric values in x_pos")

      # check for missing values
      if (any(is.na(x_zero))) stop("Missing values in x_zero")
      if (any(is.na(x_pos))) stop("Missing values in x_pos")

    } else {
      if (!.ident(length(y),xmat_rows[1])) {
        stop(
          "Data size mismatch. The following should all be equal:
          length(y): ", length(y), "\n",
          "nrow(x): ", nrow(x), "\n"
        )
      }
    }
    
  } else {

    # check if input sizes are correct
    if (
      !.ident(
        length(y),
        nrow(x),
        nrow(x_zero),
        nrow(x_pos)
      )
    ) {
      stop(
        "Data size mismatch. The following should all be equal:
        length(y): ", length(y), "\n",
        "nrow(x): ", nrow(x), "\n",
        "nrow(x_zero): ", nrow(x_zero), "\n",
        "nrow(x_pos): ", nrow(x_pos), "\n"
      )
    }

    # non-numeric values
    if (any(!is.finite(y))) stop("Non-numeric values in y")
    if (any(!is.finite(x))) stop("Non-numeric values in x")
    if (any(!is.finite(x_zero))) stop("Non-numeric values in x_zero")
    if (any(!is.finite(x_pos))) stop("Non-numeric values in x_pos")

    # missing values
    if (any(is.na(y))) stop("Missing values in y")
    if (any(is.na(x))) stop("Missing values in x")
    if (any(is.na(x_zero))) stop("Missing values in x_zero")
    if (any(is.na(x_pos))) stop("Missing values in x_pos")
  } 

  # check that y are non-negative counts
  if (any(y < 0)) stop("Negative values in y")
  if (!isTRUE(all(y == floor(y)))) stop("y must only contain integer values")

  # if(length(unique(y))<5) warning("y appears to be discrete")

  if (nburn < 0) stop("nburn must be positive")
  if (nsim < 0) stop("nsim must be positive")
  if (nthin < 0) stop("nthin must be positive")
  if (nthin > nsim + 1) stop("nthin must be < nsim")
  if (nburn < 100) warning("A low (<100) value for nburn was supplied")


  # number of trees
  if (length(ntree) != num_designs){
    ntree_vec <- rep(ntree[1], num_designs)
  } else {
    ntree_vec <- ntree
  }

  ### prepare inputs for countbart

  # offset vector
  if (any(is.null(offset))){
    # no offset specified; center the prior for the regression function at mean(y)
    mu0 <- rep(0, length(y))
  } else if (length(offset) == 1){
    # single offset value
    mu0 <- rep(log(offset), length(y))
  } else if (length(offset) == length(y)){
    # offset vector matches length of y
    mu0 <- log(offset)
  } else {
    # error
    stop("offset incorrectly specified: length(offset) does not match length(y)")
  }


  # sort y (from low to high)
  order_vec <- order(y)
  y_sort <- y[order_vec]
  mu0_sort <- mu0[order_vec]

  # convert x to matrices
  if (num_designs > 1){
    # X is a list
    X <- list()
    X_sort <- list()
    for (k in 1:num_designs) {
      X[[k]] <- matrix(x[[k]], ncol = ncol(x[[k]]))
      X_sort[[k]] <- X[[k]][order_vec, ]
    }
  } else {
    # X is a matrix
    X <- matrix(x, ncol = ncol(x))
    X_sort <- X[order_vec, ]
  }

  # determine leaf hyperparameters

  if (is.na(a0)){
    y_star <- quantile(y_sort, probs = 0.95)
    a0 <- 0.5 * (log(y_star) - mean(mu0_sort))
  }

  log_leaf_var <- a0^2 / sum(ntree_vec)
  c_start <- sum(ntree_vec) / a0^2 + 0.5
  opt_res <- optim(c_start, fn = function(x) {
    (trigamma(x) - log_leaf_var)^2
  }, method = "Brent", lower = 0, upper = 1000000000)
  leaf_c <- opt_res$par
  leaf_d <- exp(digamma(leaf_c))
  
  if ((model_type == 3) | (model_type == 4)){
    # zero-inflated model
    X_zero = matrix(x_zero, ncol = ncol(x_zero))
    X_pos = matrix(x_pos, ncol = ncol(x_pos))

    X_zero_s <- X_zero[order_vec, ]
    X_pos_s <- X_pos[order_vec, ]

    z_leaf_var <- z_conc^2 / nztree
    z_c0 <- nztree / z_conc^2 + 0.5
    z_opt_res <- optim(z_c0, fn = function(x) {
      (trigamma(x) - z_leaf_var)^2
    }, method = "Brent", lower = 0, upper = 1000000000)
    z_c <- z_opt_res$par
    z_d <- exp(digamma(z_c))
  } else {
    z_c <- 0
    z_d <- 0
  }
  
  # used to specify BART design matrices...
  muy = mean(y_sort)
  yscale = scale(y_sort)
  sdy = sd(y_sort)
  n = length(y_sort)

  # old code... used as 'ghost' inputs to 'make_bart_spec'
  if(is.null(lambda)) {
    if(is.null(sighat)) {
      if (num_designs > 1){
        lmf = lm(yscale~as.matrix(X_sort[[1]]))
      } else {
        lmf = lm(yscale~as.matrix(X_sort))
      }
      sighat = summary(lmf)$sigma*sdy #sd(y) #summary(lmf)$sigma
    }
    qchi = qchisq(1.0-sigq,nu)
    lambda = ((sighat/sdy)^2*qchi)/nu
  }
  con_sd = ifelse(abs(2*sdy - sd_val)<1e-6, 2, sd_val/sdy)
  Sigma0_con = matrix(con_sd*con_sd/sum(ntree_vec), nrow=1)


  ### prepare design structures for 'countbart'
  if (num_designs > 1){
    X_list <- X_sort

  } else {
    X_list <- list(X_sort)
  }

  basis_matrix_list <- list()
  for (k in 1:num_designs){
    basis_matrix_list[[k]] <- matrix(rep(1, n), ncol = 1)
  }

  if ((model_type == 3) | model_type == 4) {
    # add details for zero-inflated models
    X_list[[num_designs + 1]] <- X_zero_s
    X_list[[num_designs + 2]] <- X_pos_s

    basis_matrix_list[[num_designs + 1]] <- matrix(rep(1, n), ncol = 1)
    basis_matrix_list[[num_designs + 2]] <- matrix(rep(1, n), ncol = 1)

    ntree_vec <- c(ntree_vec, nztree, nztree)
  }

  # call to 'make_bart_designs' (in multibart_objects.R)
  designs <- countbcf:::make_bart_designs(X_list, basis_matrix_list)

  # call to 'make_bart_spec' (in multibart_objects.R)

  specs <- list()
  for (k in 1:length(designs)){
    specs[[k]] <- countbcf:::make_bart_spec(
       design = designs[[k]],
       ntree = ntree_vec[k],
       Sigma0 = Sigma0_con,
       scale_df = 3,
       vanilla = TRUE,
       alpha = base,
       beta = power
     )
  }
  
  # arguments to 'countbart'
  const_args <-
    list(
      y = y_sort,
      offset = mu0_sort,
      bart_specs = specs,
      bart_designs = designs,
      random_des = randeff_design,
      random_var_ix = randeff_variance_component_design,
      random_var = matrix(rep(0.00000001, ncol(randeff_variance_component_design)), ncol = 1),
      random_var_df = randeff_df,
      randeff_scales = randeff_scales / sdy,
      burn = nburn, nd = nsim, thin = nthin,
      lambda = lambda, nu = nu,
      kappa_a = kappa_a, kappa_b = kappa_b,
      kappa_prop_sd = kappa_prop_sd,
      leaf_c = leaf_c, leaf_d = leaf_d,
      z_c = z_c, z_d = z_d,
      status_interval = update_interval,
      text_trace = TRUE,
      count_model = model_type,
      return_trees = as.logical(return_trees)
    )
  
  # call 'countbart'
  fitcounts <- do.call(countbcf:::countbart, const_args)

  ### output
  out <- list(
    yhat = fitcounts$yhat_post,
    order_vec = order_vec
  )

  ## per-forest in-sample posterior means on the log scale (nsim x n,
  ## original observation order). Available regardless of return_trees.
  inv <- order(order_vec)
  n_sorted <- length(order_vec)
  .extract_coef <- function(cc) {
    nsim <- length(cc) / n_sorted
    arr <- matrix(as.numeric(cc), nrow = n_sorted, ncol = nsim)
    t(arr[inv, , drop = FALSE])
  }
  out$f_post <- .extract_coef(fitcounts$coefs[[1]])
  if ((model_type == 3) | (model_type == 4)) {
    out$f0_post <- .extract_coef(fitcounts$coefs[[num_designs + 1]])
    out$f1_post <- .extract_coef(fitcounts$coefs[[num_designs + 2]])
  }

  if (return_trees) {
    for (k in 1:num_designs){
      fit <- list(
        tree_samples = fitcounts$tree_trace[[k]],
        str = fitcounts$tree_trace[[k]]$save_string(),
        c = leaf_c, d = leaf_d,
        scale = sdy, shift = muy
      )
      nm <- if (num_designs > 1) paste0("tree_fit", k) else "tree_fit"
      out[[nm]] <- fit
    }

    if ((model_type == 3) | model_type == 4) {
      out$f0_fit <- list(
        tree_samples = fitcounts$tree_trace[[num_designs + 1]],
        str = fitcounts$tree_trace[[num_designs + 1]]$save_string(),
        c = z_c, d = z_d, scale = sdy, shift = muy
      )
      out$f1_fit <- list(
        tree_samples = fitcounts$tree_trace[[num_designs + 2]],
        str = fitcounts$tree_trace[[num_designs + 2]]$save_string(),
        c = z_c, d = z_d, scale = sdy, shift = muy
      )
    }
  }

  if ((model_type == 2) | (model_type == 4)){
    # negative binomial dispersion parameter
    out$kappa = fitcounts$kappa
    out$kappa_acpt = fitcounts$kappa_acceptance
    out$kappa_a = kappa_a
    out$kappa_b = kappa_b
  }
  
  if(debug) {
    out$const_args = const_args
    out$raw_fit = fitcounts
  }
  
  class(out) <- "count_fit"
  return(out)
}


