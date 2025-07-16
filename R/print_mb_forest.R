#' Print a multibart forest
#'
#' @param x A multibart forest
#'
#' @return
#' @export
print.mb_forest = function(x) {
  cat("Multibart forest with ", x$tree_samples$ntree, " trees")
}