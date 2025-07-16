#' make_forest
#'
#' @param tree_samples 
#' @param tree_strings 
#' @param shift 
#' @param scale 
#'
#' @return
#' @export
make_forest = function(tree_samples, tree_strings, shift=0, scale=1) {
  
  forest = list(tree_samples = tree_samples,
                     str = tree_strings,
                     shift=shift,
                     scale=scale)
  class(forest) = "mb_forest"
  
  return(forest)
}