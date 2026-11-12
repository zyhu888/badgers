## Shared scalar checks keep invalid controls out of compiled code.
check_scalar <- function(x, name, lower = -Inf, upper = Inf,
                         integer = FALSE, open_lower = FALSE,
                         open_upper = FALSE) {
  valid <- is.numeric(x) && !is.complex(x) && length(x) == 1L &&
    is.null(dim(x)) && is.finite(x)
  if (valid) {
    valid <- x >= lower && x <= upper &&
      (!open_lower || x > lower) && (!open_upper || x < upper) &&
      (!integer || x == floor(x))
  }
  if (!valid) {
    bounds <- sprintf("%s%s, %s%s", if (open_lower) "(" else "[",
                      lower, upper, if (open_upper) ")" else "]")
    stop(sprintf("%s must be a finite %sscalar in %s.", name,
                 if (integer) "integer-valued " else "numeric ", bounds),
         call. = FALSE)
  }
  invisible(NULL)
}

check_flag <- function(x, name) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    stop(sprintf("%s must be TRUE or FALSE.", name), call. = FALSE)
  }
}

check_mcmc_controls <- function(n_iter, burnin, k, a_sigma, b_sigma,
                                a_iota, b_iota, a_rho, b_rho,
                                eta, epsilon, c_graph, lite, verbose) {
  check_scalar(n_iter, "n_iter", 1, .Machine$integer.max, integer = TRUE)
  check_scalar(burnin, "burnin", 0, n_iter, integer = TRUE, open_upper = TRUE)
  check_scalar(k, "k", 1, .Machine$integer.max, integer = TRUE)
  positive <- list(a_sigma = a_sigma, b_sigma = b_sigma, a_iota = a_iota,
                   b_iota = b_iota, a_rho = a_rho, b_rho = b_rho,
                   eta = eta, epsilon = epsilon)
  for (nm in names(positive)) check_scalar(positive[[nm]], nm, 0, open_lower = TRUE)
  check_scalar(c_graph, "c_graph", 0, 6, open_lower = TRUE, open_upper = TRUE)
  check_flag(lite, "lite")
  check_flag(verbose, "verbose")
}

check_grid <- function(x, name, positive = FALSE) {
  if (!is.numeric(x) || is.complex(x) || !is.null(dim(x)) ||
      !length(x) || any(!is.finite(x)) || (positive && any(x <= 0))) {
    stop(sprintf("%s must be a nonempty finite numeric vector%s.", name,
                 if (positive) " with strictly positive values" else ""),
         call. = FALSE)
  }
}

## Input checks shared by badgers_mcmc() and badgers().
##
## Wtrain, Wtune and Wtest are num_rep x n x p arrays. Predictors may have
## different numbers of replicates q_j: the first q_j rows of W[, , j] hold the
## replicates and the remaining rows are NA, with the same pattern for every
## observation and every data set.
check_badgers_inputs <- function(ytrain, Wtrain, ytune, Wtune, ytest, Wtest,
                                 adj_graph) {

  sets <- list(train = list(y = ytrain, W = Wtrain),
               tune = list(y = ytune, W = Wtune),
               test = list(y = ytest, W = Wtest))

  for (nm in names(sets)) {
    W <- sets[[nm]]$W
    y <- sets[[nm]]$y
    if (!is.numeric(W) || is.complex(W) || length(dim(W)) != 3L ||
        any(dim(W) < 1L)) {
      stop(sprintf("W%s must be a numeric num_rep x n x p array.", nm))
    }
    if (!is.numeric(y) || is.complex(y) || !is.null(dim(y)) ||
        length(y) != dim(W)[2]) {
      stop(sprintf("y%s must have length equal to dim(W%s)[2].", nm, nm))
    }
    if (any(!is.finite(y))) stop(sprintf("y%s must contain only finite values.", nm))
    if (any(is.infinite(W)) || any(is.nan(W))) {
      stop(sprintf("W%s must contain finite values or trailing NA padding.", nm))
    }
  }

  R <- dim(Wtrain)[1]
  p <- dim(Wtrain)[3]
  for (nm in c("tune", "test")) {
    d <- dim(sets[[nm]]$W)
    if (d[1] != R || d[3] != p) {
      stop(sprintf(
        "W%s must have the same number of replicate rows and predictors as Wtrain.",
        nm))
    }
  }

  if (!is.matrix(adj_graph) || !(is.numeric(adj_graph) || is.logical(adj_graph)) ||
      is.complex(adj_graph) || any(dim(adj_graph) != p)) {
    stop("adj_graph must be a p x p matrix, with p = dim(Wtrain)[3].")
  }

  if (anyNA(adj_graph) || any(!adj_graph %in% c(0, 1)) ||
      !all(adj_graph == t(adj_graph))) {
    stop("adj_graph must be symmetric with finite 0/1 entries.", call. = FALSE)
  }

  ## Replicate layout: left-packed NA padding, identical across observations.
  miss1 <- matrix(is.na(Wtrain[, 1, , drop = FALSE]), nrow = R, ncol = p)
  if (any(miss1[1, ])) {
    stop("Every predictor needs at least one replicate (first row of W must not be NA).")
  }
  if (any(apply(miss1, 2, function(m) any(diff(as.integer(m)) < 0)))) {
    stop("Replicates must be left-packed: for each predictor, NA padding may only follow the observed replicates.")
  }
  for (nm in names(sets)) {
    W <- sets[[nm]]$W
    n <- dim(W)[2]
    expected <- aperm(array(miss1, dim = c(R, p, n)), c(1, 3, 2))
    if (!all(is.na(W) == expected)) {
      stop(sprintf(
        "The NA pattern of W%s must be the same for every observation%s: each predictor needs the same number of replicates throughout.",
        nm, if (nm == "train") "" else " and match Wtrain"))
    }
  }

  invisible(TRUE)
}


#' Run the BADGERS MCMC Sampler
#'
#' @description Runs the fixed-\code{(mu, nu)} BADGERS Markov chain Monte Carlo sampler.
#'
#' @param ytrain Numeric vector of training responses, of length \code{n}.
#'
#' @param Wtrain Three-dimensional array of replicate measurements of the training predictors, with dimension \code{num_rep} by \code{n} by \code{p}. Predictors may have different numbers of replicates \eqn{q_j}: the first \eqn{q_j} rows of \code{Wtrain[, , j]} hold the replicates and the remaining rows are \code{NA}, with the same pattern for every observation.
#'
#' @param ytune Numeric vector of validation (tuning) responses.
#'
#' @param Wtune Replicate array of the validation predictors, with dimension \code{num_rep} by \code{ntune} by \code{p} and the same \code{NA} pattern as \code{Wtrain}.
#'
#' @param ytest Numeric vector of test responses.
#'
#' @param Wtest Replicate array of the test predictors, with dimension \code{num_rep} by \code{ntest} by \code{p} and the same \code{NA} pattern as \code{Wtrain}.
#'
#' @param adj_graph Symmetric finite numeric or logical \code{p} by \code{p} 0/1 adjacency matrix. The diagonal is ignored.
#'
#' @param mu Finite common mean of the log-shrinkage parameters \eqn{\alpha_j = \log \lambda_j}, where \eqn{\beta_j} has a Laplace prior with scale \eqn{\sigma / \lambda_j}.
#'
#' @param nu Strictly positive finite scale in \eqn{\alpha \sim N(\mu 1, \nu \Phi^{-1})}.
#'
#' @param n_iter Positive integer number of MCMC iterations.
#'
#' @param burnin Nonnegative integer strictly smaller than \code{n_iter}. Summaries use iterations \code{burnin + 1} through \code{n_iter}.
#'
#' @param k Positive integer adaptation interval. Adaptation continues throughout the chain; \code{k > n_iter} disables adaptation.
#'
#' @param a_sigma,b_sigma Shape and scale of the inverse-Gamma prior on the residual variance \eqn{\sigma^2}.
#'
#' @param a_iota,b_iota Shape \code{a_iota} of the Gamma prior on the shrinkage parameters \eqn{\iota_{jk}} of the latent-covariate precision matrix, and its rate \code{b_iota} for the diagonal entries. The rates of the off-diagonal entries are determined by the working graph and \code{c_graph}.
#'
#' @param a_rho,b_rho Shape and scale of the inverse-Gamma prior on the measurement error variance \eqn{\rho^2}.
#'
#' @param eta,epsilon Shape control \eqn{\eta} and scale control \eqn{\zeta} (called \code{epsilon} here) of the constrained Wishart prior on the precision matrix of the log-shrinkage parameters.
#'
#' @param c_graph Strength of the graph information, a scalar strictly between 0 and 6. The off-diagonal rates of the \eqn{\iota_{jk}} prior are \eqn{\max(10^{-6}, d_{jk}^{-1} 10^{-6 + c})}, where \eqn{d_{jk}} is the shortest-path distance between predictors \eqn{j} and \eqn{k} in \code{adj_graph}, and \eqn{10^{-6}} if they are not connected.
#'
#' @param lite Logical; if \code{TRUE}, full parameter traces are not stored.
#'
#' @param verbose Logical; if \code{TRUE}, print the iteration number every \code{k} iterations. The default is silent.
#'
#' @returns A list with the following components:
#' \describe{
#' \item{beta_mean}{Length-\code{p} post-burn-in posterior coefficient means.}
#' \item{CI95, CI90, CI80}{Length-\code{p} thresholded posterior means; see Details.}
#' \item{MSPE_tune, MSPE_test}{Mean squared errors of posterior mean latent linear predictions on validation and test observations.}
#' \item{BETA}{\code{p} by \code{n_iter} coefficient draws, including burn-in, even when \code{lite = TRUE}.}
#' \item{mu, nu, c_graph}{The supplied hyperparameter values.}
#' \item{acs_alpha}{\code{p} by \code{n_iter} Metropolis acceptance indicators.}
#' \item{r_alpha}{\code{p} by \code{floor(n_iter / k)} acceptance fractions per adaptation block.}
#' \item{alpha_sd_prop}{Length-\code{p} final proposal standard deviations.}
#' \item{time}{Elapsed time in minutes.}
#' \item{Xtrain, Xtune, Xtest}{With \code{lite = FALSE}, sample-size by \code{p} by \code{n_iter} latent-covariate arrays. Otherwise \code{NULL}.}
#' \item{IOTA, GAMMA}{With \code{lite = FALSE}, \code{p} by \code{p} by \code{n_iter} penalty and latent precision arrays. Otherwise \code{NULL}.}
#' \item{SIGMA2, RHO2}{With \code{lite = FALSE}, length-\code{n_iter} variance draws. Otherwise \code{NULL}.}
#' \item{ALPHA, TAU}{With \code{lite = FALSE}, \code{p} by \code{n_iter} log-shrinkage and local scale draws (\code{TAU} stores scales, not squared scales). Otherwise \code{NULL}.}
#' }
#'
#' @details The MSPE values are computed from the posterior mean of
#' \eqn{x_i^\top \beta} over the post-burn-in draws, using the latent covariates
#' sampled for the validation and test observations. The sampler is silent
#' unless \code{verbose = TRUE}.
#'
#' The inputs are checked before sampling: the arrays must have matching
#' dimensions, and the \code{NA} padding of the replicates must be left-packed and
#' identical for all observations.
#'
#' For large grids, \code{badgers_mcmc()} can be run for each \code{(mu, nu)} pair
#' in parallel. A complete template that works on Windows, macOS and Linux is
#' installed with the package and can be found with
#' \code{system.file("examples", "parallel-runs.R", package = "badgers")}; it is
#' also described in the parallel-runs vignette.
#'
#' Responses and observed replicate entries must be finite. Each array dimension
#' must be positive. Only trailing \code{NA} padding is allowed; \code{NaN},
#' \code{Inf}, and missing responses are rejected. Prior shape/rate/scale controls
#' \code{a_sigma}, \code{b_sigma}, \code{a_iota}, \code{b_iota}, \code{a_rho},
#' \code{b_rho}, \code{eta}, and \code{epsilon} must be finite and strictly positive.
#' An inverse-Gamma(shape, scale) variable is defined as the reciprocal of a
#' Gamma(shape, rate = scale) variable.
#'
#' There is no intercept or automatic centering/scaling. The latent-covariate
#' prior has mean zero. Preprocessing should be estimated on training data and
#' applied consistently to held-out data. Validation and test responses are used
#' only for scoring, not for latent-covariate sampling. Select tuning parameters
#' using validation MSPE and reserve test MSPE for final evaluation.
#'
#' The \code{CI95}, \code{CI90}, and \code{CI80} fields are thresholded posterior
#' means, not interval bounds or frequentist confidence intervals. They use
#' equal-tailed empirical quantiles (R's \code{quantile(type = 7)}). No standard
#' errors or convergence diagnostics are computed. Short examples illustrate the
#' interface; assess mixing, effective sample sizes, and sensitivity with longer,
#' independent chains before inference. Iteration counts alone do not establish
#' convergence.
#'
#' @references Hu, Z., Qu, K., Long, Q., Yan, J., and Li, W. A Graph-Guided Bayesian Regression Model for High-Dimensional Data with Replicated Covariates and Measurement Error. Manuscript.
#'
#' @examples
#' library(badgers)
#'
#' sim <- dg(
#'   seed = 123,
#'   p = 8,
#'   g = 2,
#'   gsm = 2,
#'   ntrain = 12,
#'   ntune = 12,
#'   ntest = 12,
#'   sigma2 = 2,
#'   rho2 = 2,
#'   num_rep = 3,
#'   qmin = 2,
#'   qmax = 3,
#'   edii = 0.6,
#'   ediu = 0,
#'   eduu = 0.3,
#'   Gmode = 0,
#'   Gthres = 1
#' )
#'
#' fit <- badgers_mcmc(
#'   ytrain = sim$y, Wtrain = sim$W,
#'   ytune = sim$ytune, Wtune = sim$Wtune,
#'   ytest = sim$ytest, Wtest = sim$Wtest,
#'   adj_graph = sim$graph,
#'   mu = 3,
#'   nu = 3,
#'   n_iter = 100,
#'   burnin = 50
#' )
#'
#' cat("Tuning MSPE:", fit$MSPE_tune, "\n")
#' cat("Test MSPE:", fit$MSPE_test, "\n")
#' cat("True nonzero coefficients:", which(sim$beta != 0), "\n")
#' cat("Selected by 95% interval:", which(fit$CI95 != 0), "\n")
#' cat("Selected by 90% interval:", which(fit$CI90 != 0), "\n")
#' cat("Selected by 80% interval:", which(fit$CI80 != 0), "\n")
#'
#' @export
badgers_mcmc <- function(
    ytrain,
    Wtrain,
    ytune,
    Wtune,
    ytest,
    Wtest,
    adj_graph,
    mu,
    nu,
    n_iter = 10000,
    burnin = 5000,
    k = 50,
    a_sigma = 1,
    b_sigma = 1,
    a_iota = 1,
    b_iota = 1e-2,
    a_rho = 1,
    b_rho = 1,
    eta = 10,
    epsilon = 0.2,
    c_graph = 2,
    lite = TRUE,
    verbose = FALSE
) {

  check_badgers_inputs(ytrain, Wtrain, ytune, Wtune, ytest, Wtest, adj_graph)
  check_mcmc_controls(n_iter, burnin, k, a_sigma, b_sigma, a_iota, b_iota,
                      a_rho, b_rho, eta, epsilon, c_graph, lite, verbose)
  check_scalar(mu, "mu")
  check_scalar(nu, "nu", 0, open_lower = TRUE)

  start.time <- Sys.time()

  out <- badgers_mcmc_cpp(
    ytrain = ytrain, Wtrain_ = Wtrain, ytune = ytune, Wtune_ = Wtune,
    ytest = ytest, Wtest_ = Wtest, adj_graph = adj_graph,
    n_iter = n_iter, burnin = burnin, k = k,
    a_sigma = a_sigma, b_sigma = b_sigma, a_iota = a_iota, b_iota = b_iota,
    a_rho = a_rho, b_rho = b_rho, eta = eta, epsilon = epsilon,
    mu = mu, nu = nu, c_graph = c_graph, lite = lite,
    verbose = verbose
  )

  end.time <- Sys.time()
  out$time <- as.numeric(difftime(end.time, start.time, units = "mins"))

  out
}


#' Fit BADGERS and Tune the Graph Hyperparameters
#'
#' @description Runs \code{badgers_mcmc()} over a grid of \code{(mu, nu)} values and returns the model with the smallest validation mean squared prediction error.
#'
#' @param ytrain Numeric vector of training responses, of length \code{n}.
#'
#' @param Wtrain Three-dimensional array of replicate measurements of the training predictors, with dimension \code{num_rep} by \code{n} by \code{p}. Predictors may have different numbers of replicates \eqn{q_j}: the first \eqn{q_j} rows of \code{Wtrain[, , j]} hold the replicates and the remaining rows are \code{NA}, with the same pattern for every observation.
#'
#' @param ytune Numeric vector of validation (tuning) responses.
#'
#' @param Wtune Replicate array of the validation predictors, with dimension \code{num_rep} by \code{ntune} by \code{p} and the same \code{NA} pattern as \code{Wtrain}.
#'
#' @param ytest Numeric vector of test responses.
#'
#' @param Wtest Replicate array of the test predictors, with dimension \code{num_rep} by \code{ntest} by \code{p} and the same \code{NA} pattern as \code{Wtrain}.
#'
#' @param adj_graph Symmetric finite numeric or logical \code{p} by \code{p} 0/1 adjacency matrix. The diagonal is ignored.
#'
#' @param mu_grid Nonempty vector of finite candidate log-shrinkage means.
#'
#' @param nu_grid Nonempty vector of finite, strictly positive candidate scales.
#'
#' @param n_iter Positive integer number of MCMC iterations.
#'
#' @param burnin Nonnegative integer strictly smaller than \code{n_iter}. Summaries use iterations \code{burnin + 1} through \code{n_iter}.
#'
#' @param k Positive integer adaptation interval. Adaptation continues throughout the chain; \code{k > n_iter} disables adaptation.
#'
#' @param a_sigma,b_sigma Shape and scale of the inverse-Gamma prior on the residual variance \eqn{\sigma^2}.
#'
#' @param a_iota,b_iota Shape \code{a_iota} of the Gamma prior on the shrinkage parameters \eqn{\iota_{jk}} of the latent-covariate precision matrix, and its rate \code{b_iota} for the diagonal entries. The rates of the off-diagonal entries are determined by the working graph and \code{c_graph}.
#'
#' @param a_rho,b_rho Shape and scale of the inverse-Gamma prior on the measurement error variance \eqn{\rho^2}.
#'
#' @param eta,epsilon Shape control \eqn{\eta} and scale control \eqn{\zeta} (called \code{epsilon} here) of the constrained Wishart prior on the precision matrix of the log-shrinkage parameters.
#'
#' @param c_graph Strength of the graph information, a scalar strictly between 0 and 6. The off-diagonal rates of the \eqn{\iota_{jk}} prior are \eqn{\max(10^{-6}, d_{jk}^{-1} 10^{-6 + c})}, where \eqn{d_{jk}} is the shortest-path distance between predictors \eqn{j} and \eqn{k} in \code{adj_graph}, and \eqn{10^{-6}} if they are not connected.
#'
#' @param verbose Logical; print one line per \code{(mu, nu)} grid point. The sampler itself is always silent here.
#'
#' @param lite Logical; if \code{TRUE}, full parameter traces are not stored.
#'
#' @returns An object of class \code{badgers} with the selected \code{best_mu} and \code{best_nu}, the full result of the best run in \code{best_result}, and the tuning and test MSPE grids \code{MSPE_tune_grid} and \code{MSPE_test_grid}. The result also includes \code{c_graph} and a character matrix \code{grid_errors}, with \code{NA} for successful grid points.
#'
#' @details \code{badgers()} is the primary interface for fitting BADGERS (Bayesian
#' ADaptive Graph-guided Errors-in-variables Regression Solution). For each
#' combination of \code{mu_grid} and \code{nu_grid}, the function calls
#' \code{badgers_mcmc()}, records the validation (tuning) and test mean squared
#' prediction errors (MSPE), and keeps the run with the smallest tuning MSPE.
#' Invalid shared inputs are rejected before the grid runs. Failed or non-finite
#' grid points are skipped with a warning and recorded in \code{grid_errors}; an
#' error reports the first failure if no usable grid point remains. Ties select
#' the first point in traversal order (\code{mu} outer, \code{nu} inner).
#'
#' The outcome is modeled as a linear function of latent covariates
#' \eqn{x_{ij}}, which are observed through replicates
#' \eqn{w_{ijr} = x_{ij} + \kappa_{ijr}} with a common measurement error
#' variance \eqn{\rho^2}. The coefficients have Laplace priors with
#' predictor-specific rate multipliers \eqn{\lambda_j}, and the log-shrinkage
#' parameters \eqn{\alpha_j = \log \lambda_j} follow
#' \eqn{N(\mu 1, \nu \Phi^{-1})}. The precision matrix \eqn{\Phi} has a
#' constrained Wishart prior whose off-diagonal entries are zero for pairs of
#' predictors without a direct edge in \code{adj_graph}, so that connected
#' predictors receive related shrinkage. The hyperparameters \code{mu} and
#' \code{nu} are tuned on the grid.
#'
#' The defaults of the other arguments are the settings used in the paper:
#' 10000 iterations with 5000 burn-in iterations, and the coarse grid
#' \code{c(1, 3, 5, 7, 9)} for both \code{mu} and \code{nu}. The examples use 100 iterations
#' with 50 burn-in iterations to keep the run time short. A fine grid around the
#' coarse optimum can be searched in a second call.
#'
#' For large grids, \code{badgers_mcmc()} can be run for each \code{(mu, nu)} pair
#' in parallel. A complete template that works on Windows, macOS and Linux is
#' installed with the package and can be found with
#' \code{system.file("examples", "parallel-runs.R", package = "badgers")}; it is
#' also described in the parallel-runs vignette.
#'
#' The returned coefficient summaries \code{CI95}, \code{CI90} and \code{CI80} in
#' \code{best_result} are posterior means set to zero when the corresponding
#' 95\%, 90\% or 80\% credible interval covers zero. Other levels can be
#' computed from the coefficient trace \code{best_result$BETA}.
#'
#' Responses and observed replicate entries must be finite. Each array dimension
#' must be positive. Only trailing \code{NA} padding is allowed; \code{NaN},
#' \code{Inf}, and missing responses are rejected. Prior shape/rate/scale controls
#' \code{a_sigma}, \code{b_sigma}, \code{a_iota}, \code{b_iota}, \code{a_rho},
#' \code{b_rho}, \code{eta}, and \code{epsilon} must be finite and strictly positive.
#' An inverse-Gamma(shape, scale) variable is defined as the reciprocal of a
#' Gamma(shape, rate = scale) variable.
#'
#' There is no intercept or automatic centering/scaling. The latent-covariate
#' prior has mean zero. Preprocessing should be estimated on training data and
#' applied consistently to held-out data. Validation and test responses are used
#' only for scoring, not for latent-covariate sampling. Select tuning parameters
#' using validation MSPE and reserve test MSPE for final evaluation.
#'
#' The \code{CI95}, \code{CI90}, and \code{CI80} fields are thresholded posterior
#' means, not interval bounds or frequentist confidence intervals. They use
#' equal-tailed empirical quantiles (R's \code{quantile(type = 7)}). No standard
#' errors or convergence diagnostics are computed. Short examples illustrate the
#' interface; assess mixing, effective sample sizes, and sensitivity with longer,
#' independent chains before inference. Iteration counts alone do not establish
#' convergence.
#'
#' @references Hu, Z., Qu, K., Long, Q., Yan, J., and Li, W. A Graph-Guided Bayesian Regression Model for High-Dimensional Data with Replicated Covariates and Measurement Error. Manuscript.
#'
#' @examples
#' library(badgers)
#'
#' sim <- dg(
#'   seed = 123,
#'   p = 8,
#'   g = 2,
#'   gsm = 2,
#'   ntrain = 12,
#'   ntune = 12,
#'   ntest = 12,
#'   sigma2 = 2,
#'   rho2 = 2,
#'   num_rep = 3,
#'   qmin = 2,
#'   qmax = 3,
#'   edii = 0.6,
#'   ediu = 0,
#'   eduu = 0.3,
#'   Gmode = 0,
#'   Gthres = 1
#' )
#'
#' fit <- badgers(
#'   ytrain = sim$y, Wtrain = sim$W,
#'   ytune = sim$ytune, Wtune = sim$Wtune,
#'   ytest = sim$ytest, Wtest = sim$Wtest,
#'   adj_graph = sim$graph,
#'   mu_grid = c(1, 3),
#'   nu_grid = c(1, 3),
#'   n_iter = 100,
#'   burnin = 50,
#'   verbose = FALSE
#' )
#'
#' print(fit)
#'
#' cat("Test MSPE:", fit$best_result$MSPE_test, "\n")
#' cat("True nonzero coefficients:", which(sim$beta != 0), "\n")
#' cat("Selected by 95% interval:", which(fit$best_result$CI95 != 0), "\n")
#' cat("Selected by 90% interval:", which(fit$best_result$CI90 != 0), "\n")
#' cat("Selected by 80% interval:", which(fit$best_result$CI80 != 0), "\n")
#'
#' @export
badgers <- function(
    ytrain,
    Wtrain,
    ytune,
    Wtune,
    ytest,
    Wtest,
    adj_graph,
    mu_grid = c(1, 3, 5, 7, 9),
    nu_grid = c(1, 3, 5, 7, 9),
    n_iter = 10000,
    burnin = 5000,
    k = 50,
    a_sigma = 1,
    b_sigma = 1,
    a_iota = 1,
    b_iota = 1e-2,
    a_rho = 1,
    b_rho = 1,
    eta = 10,
    epsilon = 0.2,
    c_graph = 2,
    verbose = TRUE,
    lite = TRUE
) {

  check_badgers_inputs(ytrain, Wtrain, ytune, Wtune, ytest, Wtest, adj_graph)
  check_mcmc_controls(n_iter, burnin, k, a_sigma, b_sigma, a_iota, b_iota,
                      a_rho, b_rho, eta, epsilon, c_graph, lite, verbose)
  check_grid(mu_grid, "mu_grid")
  check_grid(nu_grid, "nu_grid", positive = TRUE)

  n_mu <- length(mu_grid)
  n_nu <- length(nu_grid)

  MSPE_tune_grid <- matrix(
    NA_real_, n_mu, n_nu,
    dimnames = list(mu = mu_grid, nu = nu_grid)
  )
  MSPE_test_grid <- MSPE_tune_grid

  best_mspe_tune <- Inf
  best_result <- NULL
  best_mu <- NA_real_
  best_nu <- NA_real_
  grid_errors <- matrix(NA_character_, n_mu, n_nu,
                        dimnames = dimnames(MSPE_tune_grid))

  for (a in seq_len(n_mu)) for (b in seq_len(n_nu)) {

    if (verbose) {
      cat(sprintf("mu=%.2f nu=%.2f ... ", mu_grid[a], nu_grid[b]))
    }

    out <- tryCatch(
      badgers_mcmc(
        ytrain = ytrain, Wtrain = Wtrain,
        ytune = ytune, Wtune = Wtune,
        ytest = ytest, Wtest = Wtest,
        adj_graph = adj_graph,
        n_iter = n_iter, burnin = burnin, k = k,
        a_sigma = a_sigma, b_sigma = b_sigma,
        a_iota = a_iota, b_iota = b_iota,
        a_rho = a_rho, b_rho = b_rho,
        eta = eta, epsilon = epsilon,
        mu = mu_grid[a], nu = nu_grid[b],
        c_graph = c_graph, lite = lite, verbose = FALSE
      ),
      error = function(e) {
        grid_errors[a, b] <<- conditionMessage(e)
        if (verbose) cat("ERROR:", conditionMessage(e), "\n")
        NULL
      }
    )

    if (is.null(out)) next
    if (!is.finite(out$MSPE_tune) || !is.finite(out$MSPE_test)) {
      grid_errors[a, b] <- "Sampler returned a non-finite prediction error."
      next
    }

    MSPE_tune_grid[a, b] <- out$MSPE_tune
    MSPE_test_grid[a, b] <- out$MSPE_test

    if (verbose) {
      cat(sprintf("MSPE_tune=%.3f MSPE_test=%.3f\n",
                  out$MSPE_tune, out$MSPE_test))
    }

    if (out$MSPE_tune < best_mspe_tune) {
      best_mspe_tune <- out$MSPE_tune
      best_result <- out
      best_mu <- mu_grid[a]
      best_nu <- nu_grid[b]
    }
  }

  if (is.null(best_result)) {
    stop("All grid points failed. First error: ",
         grid_errors[!is.na(grid_errors)][1], call. = FALSE)
  }
  if (any(!is.na(grid_errors))) {
    warning("Some grid points failed; inspect grid_errors in the returned fit.",
            call. = FALSE)
  }

  if (verbose) {
    cat(sprintf("\nBest: mu=%.2f nu=%.2f  MSPE_tune=%.3f MSPE_test=%.3f\n",
                best_mu, best_nu, best_mspe_tune, best_result$MSPE_test))
  }

  structure(
    list(
      best_mu = best_mu,
      best_nu = best_nu,
      c_graph = c_graph,
      best_result = best_result,
      MSPE_tune_grid = MSPE_tune_grid,
      MSPE_test_grid = MSPE_test_grid,
      grid_errors = grid_errors
    ),
    class = "badgers"
  )
}


#' Generate Batches of BADGERS Data Sets
#'
#' @description Repeatedly calls \code{dg()} and saves each data set to a file.
#'
#' @param R Number of data sets to generate.
#'
#' @param head Directory in which the data sets are saved; it must exist.
#'
#' @param seed Base random seed; data set \code{i} uses \code{seed + batch + i}.
#'
#' @param p,g,gsm,ntrain,ntune,ntest,sigma2,rho2,num_rep,qmin,qmax,edii,ediu,eduu,Gmode,Gthres Passed to \code{dg()}; see there for the defaults.
#'
#' @param batch Offset added to the data set index and the seed.
#'
#' @param verbose Logical; print progress messages.
#'
#' @returns Invisibly \code{NULL}. The data sets are saved as \code{head/data001}, \code{head/data002}, and so on, each holding one list named \code{data}.
#'
#' @details \code{R} and \code{batch} are nonnegative integers. The base seed and all
#' resulting seeds must lie in \code{[-.Machine$integer.max, .Machine$integer.max]}. Other simulation
#' constraints are described in \code{dg()}. Existing output files are overwritten;
#' use a separate directory or batch offset to preserve previous simulations.
#'
#' @examples
#' out_dir <- file.path(tempdir(), "badgers_data")
#' dir.create(out_dir)
#'
#' dg_batch(
#'   R = 2, head = out_dir, seed = 100,
#'   p = 30, g = 5, gsm = 8, ntrain = 20, ntune = 20, ntest = 20,
#'   verbose = FALSE
#' )
#'
#' load(file.path(out_dir, "data001"))
#' names(data)[1:6]
#'
#' @export
dg_batch <- function(
    R,
    head,
    seed,
    p = 200,
    g = 20,
    gsm = 10,
    ntrain = 50,
    ntune = 50,
    ntest = 50,
    sigma2 = 2,
    rho2 = 2,
    num_rep = 3,
    qmin = 4,
    qmax = 8,
    edii = 0.6,
    ediu = 0,
    eduu = 0.3,
    Gmode = 0,
    Gthres = 1,
    batch = 0,
    verbose = TRUE
) {

  check_scalar(R, "R", 0, .Machine$integer.max, integer = TRUE)
  check_scalar(batch, "batch", 0, .Machine$integer.max - R, integer = TRUE)
  check_scalar(seed, "seed", -.Machine$integer.max, .Machine$integer.max - batch - R, integer = TRUE)
  check_flag(verbose, "verbose")
  if (!is.character(head) || length(head) != 1L || is.na(head) ||
      !dir.exists(head)) stop("head must name an existing directory.", call. = FALSE)

  if (verbose) {
    cat("Generating mode:", Gmode, "\n")
    cat("Corresponding partial correlations threshold in %:", Gthres, "\n")
  }

  for (i in seq_len(R)) {
    if (verbose) cat(i, "\n")
    fname <- sprintf("%s/data%03d", head, batch + i)
    dg(seed + batch + i, p, g, gsm, ntrain, ntune, ntest, sigma2, rho2,
       num_rep, qmin, qmax, edii, ediu, eduu, Gmode, Gthres,
       savefile = fname)
  }

  invisible(NULL)
}


#' Generate BADGERS Data
#'
#' @description Generates a true graph among predictors, multivariate normal predictors with the corresponding precision structure, replicate measurements with measurement error, and a sparse linear response. A working graph, either equal to the true graph or contaminated with random edges, is returned for the sampler.
#'
#' @param seed Random seed.
#'
#' @param p Number of predictors.
#'
#' @param g Number of pathways.
#'
#' @param gsm Mean number of predictors per pathway (Poisson).
#'
#' @param ntrain Training sample size.
#'
#' @param ntune Validation (tuning) sample size.
#'
#' @param ntest Test sample size.
#'
#' @param sigma2 Residual variance.
#'
#' @param rho2 Variance of the measurement error in the replicates.
#'
#' @param num_rep Number of replicate measurements per predictor.
#'
#' @param qmin Minimum number of nonzero coefficients.
#'
#' @param qmax Maximum number of nonzero coefficients.
#'
#' @param edii Edge probability between two predictors with nonzero coefficients in the same pathway.
#'
#' @param ediu Edge probability between a predictor with a nonzero coefficient and one with a zero coefficient in the same pathway.
#'
#' @param eduu Edge probability between two predictors with zero coefficients in the same pathway.
#'
#' @param Gmode \code{0} uses the true graph as the working graph; any other value uses a noisy working graph.
#'
#' @param Gthres Fraction of true edges kept in the noisy working graph, in [0, 1]; the remaining edges are random so that the number of edges is unchanged.
#'
#' @param savefile Optional file name; if given, the generated list is saved there with \code{save()}.
#'
#' @returns A list containing the true and working graphs (\code{E0}, \code{E}, \code{graph}), the coefficients \code{beta}, and the response and replicate arrays for the training, tuning and test sets (\code{y}, \code{W}, \code{ytune}, \code{Wtune}, \code{ytest}, \code{Wtest}). The replicate arrays have dimension \code{num_rep} by sample size by \code{p}. Other fields include \code{q}, \code{pathway}, \code{ngene}, \code{pcor}, \code{nE0}, \code{nE}, latent predictors \code{X}, \code{Xtune}, \code{Xtest}, the unit-diagonal covariance \code{Cov_norm}, its inverse \code{Prec_norm}, and simulation settings.
#'
#' @details The defaults are the simulation settings of the paper: \code{p = 200}
#' predictors in \code{g = 20} pathways, three replicates per predictor, 4 to 8
#' nonzero coefficients with effect sizes from U(0.5, 2), residual variance 2, and
#' training, validation and test sets of size 50. The four scenarios of the paper
#' are obtained with \code{rho2} and the working graph as follows, using
#' \code{seed = 100 + i} (scenarios 1 and 3) or \code{seed = 200 + i} (scenarios 2
#' and 4) for data set \code{i}:
#'
#' \tabular{llll}{
#' Scenario \tab Working graph \tab \code{rho2} \tab Arguments \cr
#' 1 \tab true graph \tab 2 \tab \code{Gmode = 0} \cr
#' 2 \tab true graph \tab 5 \tab \code{Gmode = 0} \cr
#' 3 \tab 75\% of true edges plus random edges \tab 2 \tab \code{Gmode = 1, Gthres = 0.75} \cr
#' 4 \tab 75\% of true edges plus random edges \tab 5 \tab \code{Gmode = 1, Gthres = 0.75} \cr
#' }
#'
#' The generated predictors \code{X} are the latent covariates; the sampler uses
#' only the replicates \code{W} and the working graph.
#'
#' Counts \code{p}, \code{g}, sample sizes, and \code{num_rep} must be positive
#' integers. The seed is an integer in \code{[-.Machine$integer.max, .Machine$integer.max]}.
#' Sparsity bounds are integers satisfying \code{0 <= qmin <= qmax <= p};
#' \code{gsm} and both variances are finite and nonnegative. Edge probabilities
#' and \code{Gthres} lie in [0, 1]. Pathway sizes remain Poisson draws; a draw
#' exceeding \code{p} raises an informative error instead of truncating or
#' resampling it. Increase \code{p} or reduce \code{gsm} in that case.
#'
#' Noisy graphs retain \code{round(Gthres * nE0 / 2)} undirected edges with the
#' largest absolute partial correlations. The remaining edges are sampled
#' uniformly without replacement from pairs outside the retained set. They may
#' include discarded true edges, so \code{Gthres} is a retained fraction, not
#' necessarily the final graph overlap. Both graphs are symmetric with a zero
#' diagonal and equal edge counts. Fixing noisy-graph sampling changes its
#' seeded results relative to earlier versions.
#'
#' This function calls \code{set.seed(seed)} and advances the session RNG state.
#' With zero variances it returns noiseless responses or replicates.
#'
#' @examples
#' sim <- dg(
#'   seed = 1,
#'   p = 30,
#'   g = 5,
#'   gsm = 8,
#'   ntrain = 20,
#'   ntune = 20,
#'   ntest = 20
#' )
#'
#' dim(sim$W)
#' sum(sim$graph) / 2
#' which(sim$beta != 0)
#'
#' @export
dg <- function(
    seed,
    p = 200,
    g = 20,
    gsm = 10,
    ntrain = 50,
    ntune = 50,
    ntest = 50,
    sigma2 = 2,
    rho2 = 2,
    num_rep = 3,
    qmin = 4,
    qmax = 8,
    edii = 0.6,
    ediu = 0,
    eduu = 0.3,
    Gmode = 0,
    Gthres = 1,
    savefile = NULL
) {

  counts <- list(p = p, g = g, ntrain = ntrain, ntune = ntune,
                 ntest = ntest, num_rep = num_rep)
  for (nm in names(counts)) {
    check_scalar(counts[[nm]], nm, 1, .Machine$integer.max, integer = TRUE)
  }
  check_scalar(seed, "seed", -.Machine$integer.max, .Machine$integer.max, integer = TRUE)
  check_scalar(gsm, "gsm", 0)
  check_scalar(qmin, "qmin", 0, p, integer = TRUE)
  check_scalar(qmax, "qmax", qmin, p, integer = TRUE)
  check_scalar(sigma2, "sigma2", 0)
  check_scalar(rho2, "rho2", 0)
  for (nm in c("edii", "ediu", "eduu", "Gthres")) {
    check_scalar(get(nm), nm, 0, 1)
  }
  check_scalar(Gmode, "Gmode")
  if (!is.null(savefile) && (!is.character(savefile) || length(savefile) != 1L ||
                            is.na(savefile) || !nzchar(savefile))) {
    stop("savefile must be NULL or a nonempty file name.", call. = FALSE)
  }

  set.seed(seed)

  ## ---- Number of nonzero coefficients, uniform on qmin:qmax ----
  q <- qmin + sample.int(qmax - qmin + 1L, 1L) - 1L

  n <- ntrain
  pathway <- matrix(FALSE, p, g)
  ngene <- rpois(g, gsm)
  ngene[1] <- max(ngene[1], q)
  if (any(ngene > p)) {
    stop("A Poisson pathway size exceeds p; increase p or reduce gsm.", call. = FALSE)
  }
  pathway[seq_len(ngene[1]), 1] <- TRUE
  for (i in seq_len(g)[-1L]) {
    pathway[sample.int(p, ngene[i]), i] <- TRUE
  }

  ## ---- True graph G0: edges inside pathways ----
  E0 <- matrix(0, 0, 2)

  for (i in 1:g) {
    Eg <- matrix(0, 0, 2)
    mbr <- which(pathway[, i])
    for (j in mbr) {
      for (k in mbr) {
        if (j <= k) break
        if (j <= q) {
          if (runif(1) < edii | k + 1 == j) {
            Eg <- rbind(Eg, matrix(c(j, k, k, j), 2, 2))
          }
        } else if (k <= q) {
          if (runif(1) < ediu) {
            Eg <- rbind(Eg, matrix(c(j, k, k, j), 2, 2))
          }
        } else {
          if (runif(1) < eduu) {
            Eg <- rbind(Eg, matrix(c(j, k, k, j), 2, 2))
          }
        }
      }
    }
    E0 <- rbind(E0, Eg)
  }
  E0 <- unique(E0[order(E0[, 1], E0[, 2]), ])
  nE0 <- nrow(E0)

  idx <- E0[, 1] + p * (E0[, 2] - 1)

  ## ---- Precision matrix A, strictly diagonally dominant ----
  A <- diag(1, p, p)
  A[idx] <- sample(c(-1, 1), nE0, TRUE)
  A <- A * t(A)
  A[seq_len(q), seq_len(q)] <- abs(A[seq_len(q), seq_len(q), drop = FALSE])

  ## Off-diagonals are divided by 1.1 * max(sumA_j, sumA_k), where sumA is
  ## degree + 1, so row j has sum_k |A_jk| < 1 and A is positive definite.
  sumA <- apply(abs(A), 1, sum)
  maxsum <- pmax(rep(sumA, p), rep(sumA, each = p)) * 1.1
  A <- -A / maxsum
  diag(A) <- 1

  pcor <- -A[idx] / sqrt(A[E0[, 1] + p * (E0[, 1] - 1)] *
                           A[E0[, 2] + p * (E0[, 2] - 1)])

  ## ---- Normalized covariance (unit diagonal) from the precision matrix ----
  chol_A <- chol(A)
  chol_A_t <- backsolve(chol_A, diag(p))
  dchol_A_t <- apply(chol_A_t^2, 1, sum)
  chol_A_t_norm <- chol_A_t / sqrt(dchol_A_t)
  Prec_norm <- solve(tcrossprod(chol_A_t_norm))
  Cov_norm <- tcrossprod(chol_A_t_norm)

  beta <- c(runif(q, 0.5, 2), rep(0, p - q))

  ## ---- Training, tuning and test data ----
  X <- matrix(rnorm(n * p), n) %*% t(chol_A_t_norm)
  W <- apply(X, c(1, 2), function(x) {
    x + rnorm(num_rep, 0, sd = sqrt(rho2))
  })
  W <- array(W, dim = c(num_rep, n, p))
  y <- drop(X %*% beta + sqrt(sigma2) * rnorm(n))

  Xtune <- matrix(rnorm(ntune * p), ntune) %*% t(chol_A_t_norm)
  Wtune <- apply(Xtune, c(1, 2), function(x) {
    x + rnorm(num_rep, 0, sd = sqrt(rho2))
  })
  Wtune <- array(Wtune, dim = c(num_rep, ntune, p))
  ytune <- drop(Xtune %*% beta + sqrt(sigma2) * rnorm(ntune))

  Xtest <- matrix(rnorm(ntest * p), ntest) %*% t(chol_A_t_norm)
  Wtest <- apply(Xtest, c(1, 2), function(x) {
    x + rnorm(num_rep, 0, sd = sqrt(rho2))
  })
  Wtest <- array(Wtest, dim = c(num_rep, ntest, p))
  ytest <- drop(Xtest %*% beta + sqrt(sigma2) * rnorm(ntest))

  ## ---- Working graph G given to the sampler ----
  if (Gmode == 0) {
    E <- E0
  } else {
    ## Keep the top Gthres fraction of G0 edges ranked by |partial
    ## correlation|, then add random edges so that |G| = |G0|.
    E_tmp <- E0[E0[, 1] < E0[, 2], , drop = FALSE]
    pcor_tmp <- pcor[E0[, 1] < E0[, 2]]
    odr_tmp <- order(abs(pcor_tmp), decreasing = TRUE)
    num_remain <- round(Gthres * nrow(E_tmp))
    retained <- E_tmp[odr_tmp[seq_len(num_remain)], , drop = FALSE]
    candidates <- which(upper.tri(matrix(FALSE, p, p)), arr.ind = TRUE)
    key <- function(x) x[, 1] + p * (x[, 2] - 1)
    candidates <- candidates[!key(candidates) %in% key(retained), , drop = FALSE]
    added <- candidates[sample.int(nrow(candidates), nrow(E_tmp) - num_remain),
                        , drop = FALSE]
    E <- rbind(retained, added)
    E <- rbind(E, E[, 2:1, drop = FALSE])
    E <- unname(E[order(E[, 1], E[, 2]), , drop = FALSE])
  }
  nE <- nrow(E)

  graph <- matrix(0, nrow = p, ncol = p)
  for (e in seq_len(nrow(E))) {
    graph[E[e, 1], E[e, 2]] <- 1
  }

  data <- list(
    p = p, q = q, g = g, n = n, ntune = ntune, ntest = ntest,
    sigma2 = sigma2, edii = edii, ediu = ediu, eduu = eduu,
    Gmode = Gmode, Gthres = Gthres
  )
  data <- c(data, list(
    pathway = pathway, ngene = ngene, pcor = pcor, E0 = E0, nE0 = nE0,
    E = E, nE = nE, Prec_norm = Prec_norm, Cov_norm = Cov_norm,
    graph = graph
  ))
  data <- c(data, list(
    beta = beta, W = W, X = X, y = y,
    Wtune = Wtune, Xtune = Xtune, ytune = ytune,
    Wtest = Wtest, Xtest = Xtest, ytest = ytest
  ))

  if (!is.null(savefile)) {
    save(data, file = savefile)
  }

  data
}
