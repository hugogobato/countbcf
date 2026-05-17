#' countbcf: Bayesian Causal Forests for Count and Zero-Inflated Count Outcomes
#'
#' Bayesian Causal Forests for count and zero-inflated count outcomes
#' (CountBCF and Zero-Inflated CountBCF). Builds on the 'countbart' source
#' code of Wikle and Zigler (2023) and the log-linear BART backbone of
#' Murray (2021), combined with the BCF (mu, tau) decomposition of
#' Hahn, Murray and Carvalho (2020).
#'
#' @docType package
#' @import Rcpp
#' @importFrom stats approxfun lm qchisq quantile sd
#' @importFrom Rcpp evalCpp
#' @useDynLib countbcf
#' @name countbcf-package
#' @export tree_samples
NULL

loadModule("tree_samples", TRUE)
