#' badgers: Bayesian Graph-Guided Errors-in-Variables Regression
#'
#' Fit graph-guided Bayesian regression models with replicated predictors.
#' See [badgers()] for grid tuning, [badgers_mcmc()] for a fixed setting,
#' and [dg()] for simulation.
#'
#' @useDynLib badgers, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @importFrom stats rnorm rpois runif
"_PACKAGE"
