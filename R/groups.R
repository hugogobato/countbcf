
#' Title
#'
#' @param dummies 
#'
#' @return
#' @export
#'
#' @examples
get_group_indices = function(dummies) {
  apply(dummies, 2, function(g) Position(function(x) abs(x - 0)>1e-6, g))
}

#' Construct (redundant) dummies from a factor
#'
#' @param grouping_var A factor defining the groups
#' @param original_colnames If true, use the original factor levels as column names. Otherwise "Group" will be pre-pended to the column names of the dummy matrix
#'
#' @return A matrix of dummy variables
#' @import fastDummies
#' @export
#'
#' @examples
get_dummies = function(grouping_var, original_colnames = TRUE) {
  gpind = fastDummies::dummy_cols(data.frame(group = factor(grouping_var)),
                                  remove_selected_columns = TRUE)
  gpind = as.matrix(gpind)
  if(original_colnames==TRUE) {
    colnames(gpind) = make.names(substring(colnames(gpind), first=7))
  }
  
  return(gpind)
}


#' Construct metadata for BCF with random intercepts and slopes (on the treatment variable)
#'
#' @param groups A factor variable defining the groups
#' @param treatment Treatment assignments/levels
#' @param original_colnames If true, use the original factor levels as column names for the matrix of dummies for the groups. Otherwise "Group" will be pre-pended to the column names of the dummy matrix
#'
#' @return
#' @export
#'
#' @examples
random_intercepts_slopes = function(groups, treatment, original_colnames=TRUE) {
  #dummies = makeModelMatrixFromDataFrame(as.data.frame(factor(groups)))
  #dummies = makeModelMatrixFromDataFrame(as.data.frame(factor(groups)))
  
  dummies = get_dummies(groups, original_colnames=original_colnames)
  
  random_des = cbind(dummies, dummies*trt) #cbind(Xschool, Xschool*dat$treatment)
  Q = matrix(0, nrow=ncol(random_des), ncol=2)
  Q[1:ncol(random_des)/2,1] = 1
  Q[(1+ncol(random_des)/2):nrow(Q),2] = 1
  
  intercept_ix = get_group_indices(dummies)
  trt_ix = get_group_indices(dummies*trt)
  
  return(list(randeff_design = random_des,
              randeff_variance_component_design = as.matrix(Q),
              group_dummies = dummies,
              intercept_ix = intercept_ix,
              treatment_ix = trt_ix
              
  ))
  
}

#' Construct metadata for BCF with random intercepts 
#'
#' @param groups A factor variable defining the groups
#' @param original_colnames If true, use the original factor levels as column names for the matrix of dummies for the groups. Otherwise "Group" will be pre-pended to the column names of the dummy matrix
#'
#' @return
#' @export
#'
#' @examples
random_intercepts = function(groups, original_colnames=TRUE) {
  
  dummies = get_dummies(groups, original_colnames=soriginal_colnames)
  
  random_des = cbind(dummies)#cbind(Xschool, Xschool*dat$treatment)
  Q = matrix(1, nrow=ncol(random_des), ncol=1)
  
  dummies = get_dummies(groups, original_colnames)
  
  intercept_ix = get_group_indices(dummies)
  
  return(list(randeff_design = random_des,
              randeff_variance_component_design = as.matrix(Q),
              group_dummies = dummies,
              intercept_ix = intercept_ix
  ))
  
}

#' Get posteriors for random effects/parameters in a varying intercepts/slopes model
#'
#' @param bcf_fit 
#' @param randeff_setup 
#'
#' @return
#' @export
#'
#' @examples
get_random_intercepts_slopes_posteriors = function(bcf_fit, randeff_setup) {
  ngroups = length(randeff_setup$intercept_ix)
  
  intercept_posterior = bcf_fit$random_effects[,1:ngroups]
  #cat(dim(intercept_posterior))
  colnames(intercept_posterior) = colnames(randeff_setup$group_dummies)
  
  
  treatment_posterior = bcf_fit$random_effects[,(ngroups+1):(2*ngroups)]
  #cat(dim(treatment_posterior))
  colnames(treatment_posterior) = colnames(randeff_setup$group_dummies)
  
  intercept_sd = bcf_fit$random_effects_sd[,1]
  treatment_sd = bcf_fit$random_effects_sd[,2]
  
  return(list(intercept_posterior=intercept_posterior,
              treatment_posterior=treatment_posterior,
              intercept_sd=intercept_sd,
              treatment_sd=treatment_sd))
}
