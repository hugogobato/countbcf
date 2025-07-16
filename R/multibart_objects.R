make_bart_designs = function(X_list, basis_matrix_list, groups_list = NA) {
  num_designs = length(X_list)
  if(num_designs != length(basis_matrix_list)) {
    stop("X_list and basis_matrix_list must be the same size")
  }
  
  if (all(!is.na(groups_list))) {
    group = sapply(groups_list, function(x) length(x) > 1)
    out = lapply(1:num_designs, function(i) {
      make_bart_design(
        X = X_list[[i]], basis_matrix = basis_matrix_list[[i]],
        group = group[i], groups = groups_list[[i]], index = i - 1
      )
    })
  } else {
    out = lapply(1:num_designs, function(i) make_bart_design(X = X_list[[i]], basis_matrix = basis_matrix_list[[i]], index = i - 1))
  }
  
  return(out)
}

make_bart_design = function(
  X,
  basis_matrix,
  group = FALSE,
  groups = 0,
  index = -1
  ) {
  
  cutpoint_list = lapply(1:ncol(X), function(i) countbart:::.cp_quantile(X[,i]))
  
  list(X=t(X),
       Omega = t(basis_matrix),
       info = cutpoint_list,
       index = index,
       group = group,
       groups = groups
  )
}

make_bart_spec = function(
  design,
  ntree,
  Sigma0,
  mu0 = NA,
  scale_df = -1,
  vanilla = FALSE,
  alpha = 0.95, 
  beta = 2,
  dart = FALSE,
  update_leaf_scale = TRUE
  ) {
  
  if(is.na(mu0)) mu0 = rep(0, ncol(Sigma0))
  if(ncol(design$Omega)==1 & all( (design$Omega - 1)<1e-8)) vanilla = TRUE
  
  list(
    design_index = design$index,
    ntree = ntree,
    scale_df = scale_df,
    mu0 = mu0,
    Sigma0 = Sigma0,
    vanilla = vanilla,
    alpha = alpha,
    beta = beta,
    dart = dart,
    sample_eta = ifelse(update_leaf_scale, 1, -1)
  )
}