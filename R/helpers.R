.ident <- function(...){
  # courtesy https://stackoverflow.com/questions/19966515/how-do-i-test-if-three-variables-are-equal-r
  args <- c(...)
  if( length( args ) > 2L ){
    #  recursively call ident()
    out <- c( identical( args[1] , args[2] ) , .ident(args[-1]))
  }else{
    out <- identical( args[1] , args[2] )
  }
  return( all( out ) )
}

### Permute the rows of a random-effects design matrix so that they line up
### with y/z/X/offset after the internal sort of y (zeros first). Rows of
### randeff_design correspond to observations on the C++ side (it forms
### random_des.t() * r and random_des * gamma), so failing to permute it
### silently attributes each observation to the wrong group.
.sort_randeff_design <- function(randeff_design,
                                 randeff_variance_component_design,
                                 order_vec) {
  W <- as.matrix(randeff_design)
  Q <- as.matrix(randeff_variance_component_design)
  n <- length(order_vec)

  ### the 1x1 defaults disable random effects on the C++ side; nothing to sort
  if (length(W) == 1L && length(Q) == 1L) return(W)

  if (nrow(W) != n) {
    stop(sprintf(paste0("randeff_design has %d row(s) but y has length %d. ",
                        "Rows of randeff_design must correspond to observations."),
                 nrow(W), n))
  }
  if (length(Q) == 1L) {
    stop(paste0("randeff_design was supplied but randeff_variance_component_design ",
                "is the 1x1 default, which disables random effects. Supply a ",
                "variance-component design with one row per column of ",
                "randeff_design and one column per variance component."))
  }
  if (ncol(W) != nrow(Q)) {
    stop(sprintf(paste0("Dimension mismatch: randeff_design has %d column(s) but ",
                        "randeff_variance_component_design has %d row(s); these must agree."),
                 ncol(W), nrow(Q)))
  }

  W[order_vec, , drop = FALSE]
}

.cp_quantile = function(x, num=10000, cat_levels=8){
  nobs = length(x)
  nuniq = length(unique(x))
  
  if(nuniq==1) {
    ret = x[1]
    warning("A supplied covariate contains a single distinct value.")
  } else if(nuniq < cat_levels) {
    xx = sort(unique(x))
    ret = xx[-length(xx)] + diff(xx)/2
  } else {
    q = approxfun(sort(x),quantile(x,p = 0:(nobs-1)/nobs))
    ind = seq(min(x),max(x),length.out=num)
    ret = q(ind)
  }
  
  return(ret)
}
