#' countbart
#'
#' Description of your package
#'
#' @docType package
#' @import Rcpp
#' @importFrom stats approxfun lm qchisq quantile sd
#' @importFrom Rcpp evalCpp
#' @useDynLib countbart
#' @name countbart
#' @export tree_samples 
NULL

loadModule("tree_samples", TRUE)
