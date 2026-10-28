\name{badgers_mcmc}
\alias{badgers_mcmc}
\title{Run the BADGERS MCMC Sampler}
\usage{
badgers_mcmc(ytrain, Wtrain, ytune, Wtune, ytest, Wtest, adj_graph,
  mu, nu, n_iter = 10000, burnin = 5000, k = 50,
  a_sigma = 1, b_sigma = 1, a_iota = 1, b_iota = 1e-2,
  a_rho = 1, b_rho = 1, eta = 10, epsilon = 0.2,
  c_graph = 2, lite = TRUE, verbose = FALSE)
}
\arguments{
\item{ytrain}{Numeric vector of training responses, of length \code{n}.}
\item{Wtrain}{Three-dimensional array of replicate measurements of the training predictors, with dimension \code{num_rep} by \code{n} by \code{p}. Predictors may have different numbers of replicates \eqn{q_j}: the first \eqn{q_j} rows of \code{Wtrain[, , j]} hold the replicates and the remaining rows are \code{NA}, with the same pattern for every observation.}
\item{ytune}{Numeric vector of validation (tuning) responses.}
\item{Wtune}{Replicate array of the validation predictors, with dimension \code{num_rep} by \code{ntune} by \code{p} and the same \code{NA} pattern as \code{Wtrain}.}
\item{ytest}{Numeric vector of test responses.}
\item{Wtest}{Replicate array of the test predictors, with dimension \code{num_rep} by \code{ntest} by \code{p} and the same \code{NA} pattern as \code{Wtrain}.}
\item{adj_graph}{A \code{p} by \code{p} 0/1 adjacency matrix of the working graph among the predictors.}
\item{mu}{Common mean of the log-shrinkage parameters \eqn{\alpha_j = \log \lambda_j}, where \eqn{\beta_j} has a Laplace prior with scale \eqn{\sigma / \lambda_j}.}
\item{nu}{Scale of the log-shrinkage parameters, \eqn{\alpha \sim N(\mu 1, \nu \Phi^{-1})}.}
\item{n_iter}{Number of MCMC iterations.}
\item{burnin}{Number of burn-in iterations.}
\item{k}{Number of iterations between adaptations of the Metropolis-Hastings proposal standard deviations for the log-shrinkage parameters.}
\item{a_sigma, b_sigma}{Shape and scale of the inverse-Gamma prior on the residual variance \eqn{\sigma^2}.}
\item{a_iota, b_iota}{Shape \code{a_iota} of the Gamma prior on the shrinkage parameters \eqn{\iota_{jk}} of the latent-covariate precision matrix, and its rate \code{b_iota} for the diagonal entries. The rates of the off-diagonal entries are determined by the working graph and \code{c_graph}.}
\item{a_rho, b_rho}{Shape and scale of the inverse-Gamma prior on the measurement error variance \eqn{\rho^2}.}
\item{eta, epsilon}{Shape control \eqn{\eta} and scale control \eqn{\zeta} (called \code{epsilon} here) of the constrained Wishart prior on the precision matrix of the log-shrinkage parameters.}
\item{c_graph}{Strength of the graph information, a scalar strictly between 0 and 6. The off-diagonal rates of the \eqn{\iota_{jk}} prior are \eqn{d_{jk}^{-1} 10^{-6 + c}}, where \eqn{d_{jk}} is the shortest-path distance between predictors \eqn{j} and \eqn{k} in \code{adj_graph}, and \eqn{10^{-6}} if they are not connected.}
\item{lite}{Logical; if \code{TRUE}, full parameter traces are not stored.}
\item{verbose}{Logical; if \code{TRUE}, print the iteration number every \code{k} iterations. The default is silent.}
}
\value{
A list containing the validation and test MSPE (\code{MSPE_tune}, \code{MSPE_test}), the posterior mean \code{beta_mean}, the thresholded coefficient estimates \code{CI95}, \code{CI90} and \code{CI80}, the coefficient trace \code{BETA} (a \code{p} by \code{n_iter} matrix), the sampler adaptation summaries, and the run time in minutes (\code{time}). With \code{lite = FALSE} the traces of the latent predictors and the remaining parameters are also returned.
}
\description{
Runs the fixed-\code{(mu, nu)} BADGERS Markov chain Monte Carlo sampler.
}
\details{
The MSPE values are computed from the posterior mean of
\eqn{x_i^\top \beta} over the post-burn-in draws, using the latent covariates
sampled for the validation and test observations. The sampler is silent
unless \code{verbose = TRUE}.

The inputs are checked before sampling: the arrays must have matching
dimensions, and the \code{NA} padding of the replicates must be left-packed and
identical for all observations.

For large grids, \code{badgers_mcmc()} can be run for each \code{(mu, nu)} pair
in parallel. A complete template that works on Windows, macOS and Linux is
installed with the package and can be found with
\code{system.file("examples", "parallel-runs.R", package = "badgers")}; it is
also described in the parallel-runs vignette.
}
\examples{
library(badgers)

sim <- dg(
  seed = 123,
  p = 30,
  g = 5,
  gsm = 8,
  ntrain = 40,
  ntune = 40,
  ntest = 40,
  sigma2 = 2,
  rho2 = 2,
  num_rep = 3,
  qmin = 4,
  qmax = 6,
  edii = 0.6,
  ediu = 0,
  eduu = 0.3,
  Gmode = 0,
  Gthres = 1
)

fit <- badgers_mcmc(
  ytrain = sim$y, Wtrain = sim$W,
  ytune = sim$ytune, Wtune = sim$Wtune,
  ytest = sim$ytest, Wtest = sim$Wtest,
  adj_graph = sim$graph,
  mu = 3,
  nu = 3,
  n_iter = 1000,
  burnin = 500
)

cat("Tuning MSPE:", fit$MSPE_tune, "\n")
cat("Test MSPE:", fit$MSPE_test, "\n")
cat("True nonzero coefficients:", which(sim$beta != 0), "\n")
cat("Selected by 95\% interval:", which(fit$CI95 != 0), "\n")
cat("Selected by 90\% interval:", which(fit$CI90 != 0), "\n")
cat("Selected by 80\% interval:", which(fit$CI80 != 0), "\n")
}
\references{
Hu, Z., Qu, K., Long, Q., Yan, J., and Li, W. A Graph-Guided Bayesian Regression Model for High-Dimensional Data with Replicated Covariates and Measurement Error. Manuscript.
}
