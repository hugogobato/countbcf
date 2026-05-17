#' Fit Bayesian Causal Forests
#'
#' @references Hahn, Murray, and Carvalho(2017). Bayesian regression tree models for causal inference: regularization, confounding, and heterogeneous effects.
#'  https://arxiv.org/abs/1706.09523. (Call citation("bcf") from the
#' command line for citation information in Bibtex format.)
#'
#' @details Fits a generalized version of the Bayesian Causal Forest model (Hahn et. al. 2018): For a response
#' variable y, treatment z, and covariates x,
#' \deqn{y_i = \mu(x_i, \hat z_i) + \tau(x_i, \pi_i)\omega(z_i) + \epsilon_i}
#' where \eqn{\z_i} is an (optional) estimate of \eqn{E(Z_i | X_i=x_i)} and
#' \eqn{\epsilon_i \sim N(0,\sigma^2)}
#'
#' Some notes:
#' \itemize{
#'    \item x_control and x_moderate must be numeric matrices. See e.g. the makeModelMatrix function in the
#'    dbarts package for appropriately constructing a design matrix from a data.frame
#'    \item sd_control and sd_moderate are the prior SD(mu(x)) and SD(tau(x)) at a given value of x (respectively). If
#'    use_muscale = FALSE, then this is the parameter \eqn{\sigma_\mu} from the original BART paper, where the leaf parameters
#'    have prior distribution \eqn{N(0, \sigma_\mu/m)}, where m is the number of trees.
#'    If use_muscale=TRUE then sd_control is the prior median of a half Cauchy prior for SD(mu(x)). If use_tauscale = TRUE,
#'    then sd_moderate is the prior median of a half Normal prior for SD(tau(x)).
#'    \item By default the prior on \eqn{\sigma^2} is calibrated as in Chipman, George and McCulloch (2008).
#'
#'
#' }
#'
#' @param y Response variable
#' @param z Treatment assignments
#' @param x_control Design matrix for the "prognostic" function mu(x)
#' @param x_moderate Design matrix for the covariate-dependent treatment effects tau(x)
#' @param zhat Length n estimates of E(Z|X)
#' @param nburn Number of burn-in MCMC iterations
#' @param nsim Number of MCMC iterations to save after burn-in
#' @param nthin Save every nthin'th MCMC iterate. The total number of MCMC iterations will be nsim*nthin + nburn.
#' @param update_interval Print status every update_interval MCMC iterations
#' @param ntree_control Number of trees in mu(x)
#' @param sd_control SD(mu(x)) marginally at any covariate value (or its prior median if use_muscale=TRUE)
#' @param base_control Base for tree prior on mu(x) trees (see details)
#' @param power_control Power for the tree prior on mu(x) trees
#' @param ntree_moderate Number of trees in tau(x)
#' @param sd_moderate SD(tau(x)) marginally at any covariate value (or its prior median if use_tauscale=TRUE)
#' @param base_moderate Base for tree prior on tau(x) trees (see details)
#' @param power_moderate Power for the tree prior on tau(x) trees (see details)
#' @param nu Degrees of freedom in the chisq prior on \eqn{sigma^2}
#' @param lambda Scale parameter in the chisq prior on \eqn{sigma^2}
#' @param sigq Calibration quantile for the chisq prior on \eqn{sigma^2}
#' @param sighat Calibration estimate for the chisq prior on \eqn{sigma^2}
#' @param include_zhat Takes values "control", "moderate", "both" or "none". Whether to
#' include zhat in mu(x) ("control"), tau(x) ("moderate"), both or none. Values of "control"
#' or "both" are HIGHLY recommended with observational data.
#' @param use_muscale Use a half-Cauchy hyperprior on the scale of mu.
#' @param randeff_design XX
#' @param randeff_variance_component_design XX
#' @param randeff_scales XX
#' @param randeff_df XX
#' @param debug XX
#' @param use_tauscale Use a half-Normal prior on the scale of tau.
#'
#' @return A list with elements
#' \item{tau}{\code{nsim} by \code{n} matrix of posterior samples of individual treatment effects}
#' \item{mu}{\code{nsim} by \code{n} matrix of posterior samples of individual treatment effects}
#' \item{sigma}{Length \code{nsim} vector of posterior samples of sigma}
#' 
#' @export
#' 
bcf_continuous_linear <- function(y, z, #Omega_con=matrix(rep(1,length(y)), ncol=1), Omega_mod=Omega_con,
                            x_control, x_moderate = x_control,
                            #Omega_con_out = Omega_con[1:2,], Omega_mod_out = Omega_mod[1:2],
                            #control_fits = FALSE, moderate_fits = FALSE, # Only returns coefficients for out of sample x's
                            #x_control_out = x_control[1:2,], x_moderate_out = x_moderate[1:2,],
                            zhat = rep(0.5, length(y)), #zhat_out = rep(0.5, nrow(x_control_out)),
                            randeff_design = matrix(1),
                            randeff_variance_component_design = matrix(1),
                            randeff_scales = 1,
                            randeff_df = 3,
                            nburn, nsim, nthin = 1, update_interval = 100,
                            ntree_control = 250,
                            sd_control = 2*sd(y),
                            base_control = 0.95,
                            power_control = 2,
                            ntree_moderate = 50,
                            sd_moderate = 0.25*sd(y)/sd(z),
                            base_moderate = 0.25,
                            power_moderate = 3,
                            nu = 3, lambda = NULL, sigq = .9, sighat = NULL,
                            include_zhat = "control", 
                            use_muscale=TRUE, use_tauscale=TRUE,
                            debug = FALSE
) {
  
  #cat"into main\n")
  
  # TODO: checks for random effects matrices
  
  zhat = as.matrix(zhat)
  if( !.ident(length(y),
    #          nrow(Omega_con),
    #          nrow(Omega_mod),
              nrow(x_control),
              nrow(x_moderate),
              nrow(zhat)
  )
  ) {
    
    stop("Data size mismatch. The following should all be equal:
         length(y): ", length(y), "\n",
         #"nrow(Omega_con): ", nrow(Omega_con), "\n",
         #"nrow(Omega_mod): ", nrow(Omega_mod), "\n",
         "nrow(x_control): ", nrow(x_control), "\n",
         "nrow(x_moderate): ", nrow(x_moderate), "\n",
         "nrow(zhat): ", nrow(zhat),"\n"
    )
  }
  # 
  # if( !.ident(nrow(Omega_con_out),
  #             nrow(Omega_mod_out),
  #             nrow(x_control_out),
  #             nrow(x_moderate_out),
  #             nrow(zhat_out)
  # )
  # ) {
  #   
  #   stop("Data size mismatch. The following should all be equal:
  #        nrow(Omega_con_out): ", nrow(Omega_con_out), "\n",
  #        "nrow(Omega_mod_out): ", nrow(Omega_mod_out), "\n",
  #        "nrow(x_control_out): ", nrow(x_control_out), "\n",
  #        "nrow(x_moderate_out): ", nrow(x_moderate_out), "\n",
  #        "nrow(zhat_out): ", nrow(zhat_out),"\n"
  #   )
  # }
  # 
  if(any(is.na(y))) stop("Missing values in y")
  if(any(is.na(z))) stop("Missing values in z")
  # if(any(is.na(Omega_mod))) stop("Missing values in Omega_mod")
  if(any(is.na(x_control))) stop("Missing values in x_control")
  if(any(is.na(x_moderate))) stop("Missing values in x_moderate")
  if(any(is.na(zhat))) stop("Missing values in zhat")
  
  if(any(!is.finite(y))) stop("Non-numeric values in y")
  if(any(!is.finite(z))) stop("Non-numeric values in z")
  # if(any(!is.finite(Omega_mod))) stop("Non-numeric values in Omega_mod")
  if(any(!is.finite(x_control))) stop("Non-numeric values in x_control")
  if(any(!is.finite(x_moderate))) stop("Non-numeric values in x_moderate")
  if(any(!is.finite(zhat))) stop("Non-numeric values in zhat")
  
  # if(any(is.na(Omega_con_out))) stop("Missing values in Omega_con_out")
  # if(any(is.na(Omega_mod_out))) stop("Missing values in Omega_mod_out")
  # if(any(is.na(x_control_out))) stop("Missing values in x_control_out")
  # if(any(is.na(x_moderate_out))) stop("Missing values in x_moderate_out")
  # if(any(is.na(zhat_out))) stop("Missing values in zhat_out")
  # 
  # if(any(!is.finite(Omega_con_out))) stop("Non-numeric values in Omega_con_out")
  # if(any(!is.finite(Omega_mod_out))) stop("Non-numeric values in Omega_mod_out")
  # if(any(!is.finite(x_control_out))) stop("Non-numeric values in x_control_out")
  # if(any(!is.finite(x_moderate_out))) stop("Non-numeric values in x_moderate_out")
  # if(any(!is.finite(zhat_out))) stop("Non-numeric values in zhat_out")
  
  #if(!all(sort(unique(z)) == c(0,1))) stop("z must be a vector of 0's and 1's, with at least one of each")
  
  if(length(unique(y))<5) warning("y appears to be discrete")
  
  if(nburn<0) stop("nburn must be positive")
  if(nsim<0) stop("nsim must be positive")
  if(nthin<0) stop("nthin must be positive")
  if(nthin>nsim+1) stop("nthin must be < nsim")
  if(nburn<100) warning("A low (<100) value for nburn was supplied")
  
  ### TODO range check on parameters
  
  ###
  x_c = matrix(x_control, ncol=ncol(x_control))
  x_m = matrix(x_moderate, ncol=ncol(x_moderate))
  # x_c_o = matrix(x_control_out, ncol=ncol(x_control_out))
  # x_m_o = matrix(x_moderate_out, ncol=ncol(x_moderate_out))
  
  if(include_zhat!="none"){
    if(length(unique(zhat)) == 1){
      warning("All values of zhat are equal. zhat will not be included among covariates")
      include_zhat="none"
    }
  }
  
  if(include_zhat=="both" | include_zhat=="control") {
    x_c = cbind(x_control, zhat)
  }
  if(include_zhat=="both" | include_zhat=="moderate") {
    x_m = cbind(x_moderate, zhat)
  }
  
  yscale = scale(y)
  sdy = sd(y)
  muy = mean(y)
  n = length(y)
  
  if(is.null(lambda)) {
    if(is.null(sighat)) {
      lmf = lm(yscale~as.matrix(x_c))
      sighat = summary(lmf)$sigma*sdy #sd(y) #summary(lmf)$sigma
    }
    qchi = qchisq(1.0-sigq,nu)
    lambda = ((sighat/sdy)^2*qchi)/nu
  }
  
  con_sd = ifelse(abs(2*sdy - sd_control)<1e-6, 2, sd_control/sdy)
  mod_sd = ifelse(abs(sdy - sd_moderate)<1e-6, 1, sd_moderate/sdy)/ifelse(use_tauscale,0.674,1) # if HN make sd_moderate the prior median
  

  Sigma0_con = matrix(con_sd*con_sd/(ntree_control), nrow=1)

  Sigma0_mod = matrix(0, nrow=1, ncol = 1)
  diag(Sigma0_mod) = diag(Sigma0_mod) + mod_sd*mod_sd/(ntree_moderate)
  
  X_list = list(x_c, x_m)
  basis_matrix_list = list(matrix(rep(1, n), ncol=1), matrix(z, ncol=1))
  
  designs = countbcf:::make_bart_designs(X_list, basis_matrix_list)
  specs = list(countbcf:::make_bart_spec(design = designs[[1]], 
                                          ntree = ntree_control, 
                                          Sigma0 = Sigma0_con, 
                                          scale_df = 3, vanilla=TRUE,
                                          alpha = base_control,
                                          beta = power_control),
               countbcf:::make_bart_spec(design=designs[[2]], 
                                          ntree=ntree_moderate, 
                                          Sigma0 = Sigma0_mod, 
                                          scale_df = -1, 
                                          alpha=base_moderate, 
                                          beta=power_moderate))
  
  
  # const_args = 
  #   list(y_ = yscale, Omega_con = t(Omega_con), Omega_mod = t(Omega_mod),
  #                     Omega_con_est = t(Omega_con_out), Omega_mod_est = t(Omega_mod_out),
  #                     t(x_c), t(x_m), t(x_c_o), t(x_m_o),
  #                     cutpoint_list_c, cutpoint_list_m,
  #                     random_des = randeff_design,
  #                     random_var = matrix(rep(0.00000001, ncol(randeff_variance_component_design)), ncol=1),
  #                     random_var_ix = randeff_variance_component_design,
  #                     random_var_df = randeff_df, randeff_scales = randeff_scales/sdy,
  #                     burn = nburn, nd = nsim, thin = nthin,
  #                     ntree_mod = ntree_moderate, ntree_con = ntree_control, 
  #                     lambda = lambda, nu = nu,
  #                     Sigma0_con = Sigma0_con, Sigma0_mod =Sigma0_mod,
  #                     #Sigma0_con=2*matrix(1/250), Sigma0_mod=1*matrix(0.02),
  #                     con_alpha = base_control, con_beta = power_control,
  #                     mod_alpha = base_moderate, mod_beta = power_moderate,
  #                     treef_name = tempdir(), prior_sample = FALSE,
  #                     use_con_scale = use_muscale, use_mod_scale = use_tauscale,
  #                     con_scale_df = 1, mod_scale_df = -1,
  #                     status_interval = update_interval,
  #                     vanilla = vanilla, 
  #                     dart = dart, var_sizes_con = var_sizes_con, var_sizes_mod = var_sizes_mod
  #   )
  # 
  # 
  #cat"here\n\n")
  const_args = 
    list(y = yscale, 
         bart_specs = specs,
         bart_designs = designs,
         random_des = randeff_design,
         random_var_ix = randeff_variance_component_design,
         random_var = matrix(rep(0.00000001, ncol(randeff_variance_component_design)), ncol=1),
         random_var_df = randeff_df, 
         randeff_scales = randeff_scales/sdy,
         burn = nburn, nd = nsim, thin = nthin,
         lambda = lambda, nu = nu,
         status_interval = update_interval,
         text_trace=TRUE
    )
  
  #cat"here\n\n")

  fitbcf = do.call(countbcf:::multibart, const_args)
  
  control_fit = list(tree_samples = fitbcf$tree_trace[[1]],
                     str = fitbcf$tree_trace[[1]]$save_string(),
                     scale=sdy, shift = muy)
  moderate_fit = list(tree_samples = fitbcf$tree_trace[[2]],
                     str = fitbcf$tree_trace[[2]]$save_string(),
                     scale=sdy, shift = 0)
  
  
  #dim(fitbcf$coefs[[1]]) = c(1, n, nsim)
  #dim(fitbcf$coefs[[2]]) = c(1, n, nsim)
  
  out = list(sigma = sdy*fitbcf$sigma,
       yhat = muy + sdy*fitbcf$yhat_post,
       #mu_post = t(muy + sdy*fitbcf$coefs[[1]][,,,drop=TRUE]),
       #tau_post = t(sdy*fitbcf$coefs[[2]][,,,drop=TRUE]),
       #       mu  = m_post,
       #tau = tau_post,
       #mu_scale = fitbcf$etas[,1]*sdy,
       #tau_scale = fitbcf$etas[,2]*sdy,
       #coefs_mod = fitbcf$coefs_mod*sdy,
       #coefs_con = fitbcf$coefs_con*sdy,
       #eta_con = fitbcf$eta_con_post,
       #eta_mod = fitbcf$eta_mod_post,
       #_["coefs_mod_est"] = scoefs_mod_est,
       #_["coefs_con_est"] = scoefs_con_est,
       #coefs_mod_est = fitbcf$coefs_mod_est*sdy,
       #coefs_con_est = fitbcf$coefs_con_est*sdy,
       random_effects = fitbcf$gamma*sdy,
       random_effects_sd = fitbcf$random_sd_post*sdy,
       control_fit = control_fit,
       moderate_fit = moderate_fit
       #splitting_prob_con = fitbcf$var_prob_con,
       #splitting_prob_mod = fitbcf$var_prob_mod
  )
  
  if(debug) {
    out$const_args = const_args
    out$raw_fit = fitbcf
  }
  
  class(out) = "bcf_fit"
  
  return(out)
}
