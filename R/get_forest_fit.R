#' Extract fits/coefficients from multibart forests
#'
#' @param mb_forest A multibart forest object
#' @param X Covariate matrix
#' @param Z Basis functions corresponding to X
#'
#' @return If Z is null, a (basis dimension) x n x (number of MCMC iterations) array
#' of coefficients. If Z is not null an n x (number of MCMC iterations) matrix
#' of fitted values
#' @export
get_forest_fit = function(mb_forest, X, Z=NULL) {
  test_ptr = try({m = mb_forest$tree_samples$ntree}, silent=TRUE)
  if(class(test_ptr)=="try-error") {
    mb_forest$tree_samples = new(tree_samples)
    mb_forest$tree_samples$load_string(mb_forest$str)
  }
  if(is.null(Z)) {
    out  = mb_forest$tree_samples$coefs(as.matrix(t(X)))
    if(mb_forest$tree_samples$basis_dim==1) {
      return( t(drop(out)*mb_forest$scale))
    } else {
      return(out*mb_forest$scale)
    }
  } else {
    out  = mb_forest$tree_samples$fits(as.matrix(t(X)), as.matrix(t(Z)))
    return( t(out*mb_forest$scale + mb_forest$shift))
  }
}
