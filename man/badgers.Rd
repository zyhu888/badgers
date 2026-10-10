\name{badgers}
\alias{badgers}
\title{Fit BADGERS and Tune the Graph Hyperparameters}
\usage{
badgers(ytrain, Wtrain, ytune, Wtune, ytest, Wtest, adj_graph,
  mu_grid = c(1, 3, 5, 7, 9), nu_grid = c(1, 3, 5, 7, 9),
  n_iter = 10000, burnin = 5000, k = 50,
  a_sigma = 1, b_sigma = 1, a_iota = 1, b_iota = 1e-2,
  a_rho = 1, b_rho = 1, eta = 10, epsilon = 0.2,
  c_graph = 2, verbose = TRUE, lite = TRUE)
}
\arguments{
\item{ytrain}{Numeric vector of training responses, of length \code{n}.}
\item{Wtrain}{Three-dimensional array of replicate measurements of the training predictors, with dimension \code{num_rep} by \code{n} by \code{p}. Predictors may have different numbers of replicates \eqn{q_j}: the first \eqn{q_j} rows of \code{Wtrain[, , j]} hold the replicates and the remaining rows are \code{NA}, with the same pattern for every observation.}
\item{ytune}{Numeric vector of validation (tuning) responses.}
\item{Wtune}{Replicate array of the validation predictors, with dimension \code{num_rep} by \code{ntune} by \code{p} and the same \code{NA} pattern as \code{Wtrain}.}
\item{ytest}{Numeric vector of test responses.}
\item{Wtest}{Replicate array of the test predictors, with dimension \code{num_rep} by \code{ntest} by \code{p} and the same \code{NA} pattern as \code{Wtrain}.}
\item{adj_graph}{A \code{p} by \code{p} 0/1 adjacency matrix of the working graph among the predictors.}
\item{mu_grid}{Candidate values of \code{mu}.}
\item{nu_grid}{Candidate values of \code{nu}.}
\item{n_iter}{Number of MCMC iterations.}
\item{burnin}{Number of burn-in iterations.}
\item{k}{Number of iterations between adaptations of the Metropolis-Hastings proposal standard deviations for the log-shrinkage parameters.}
\item{a_sigma, b_sigma}{Shape and scale of the inverse-Gamma prior on the residual variance \eqn{\sigma^2}.}
\item{a_iota, b_iota}{Shape \code{a_iota} of the Gamma prior on the shrinkage parameters \eqn{\iota_{jk}} of the latent-covariate precision matrix, and its rate \code{b_iota} for the diagonal entries. The rates of the off-diagonal entries are determined by the working graph and \code{c_graph}.}
\item{a_rho, b_rho}{Shape and scale of the inverse-Gamma prior on the measurement error variance \eqn{\rho^2}.}
\item{eta, epsilon}{Shape control \eqn{\eta} and scale control \eqn{\zeta} (called \code{epsilon} here) of the constrained Wishart prior on the precision matrix of the log-shrinkage parameters.}
\item{c_graph}{Strength of the graph information, a scalar strictly between 0 and 6. The off-diagonal rates of the \eqn{\iota_{jk}} prior are \eqn{d_{jk}^{-1} 10^{-6 + c}}, where \eqn{d_{jk}} is the shortest-path distance between predictors \eqn{j} and \eqn{k} in \code{adj_graph}, and \eqn{10^{-6}} if they are not connected.}
\item{verbose}{Logical; print progress messages.}
\item{lite}{Logical; if \code{TRUE}, full parameter traces are not stored.}
}
\value{
An object of class \code{badgers} with the selected \code{best_mu} and \code{best_nu}, the full result of the best run in \code{best_result}, and the tuning and test MSPE grids \code{MSPE_tune_grid} and \code{MSPE_test_grid}.
}
\description{
Runs \code{badgers_mcmc()} over a grid of \code{(mu, nu)} values and returns the model with the smallest validation mean squared prediction error.
}
\details{
\code{badgers()} is the primary interface for fitting BADGERS (Bayesian
ADaptive Graph-guided Errors-in-variables Regression Solution). For each
combination of \code{mu_grid} and \code{nu_grid}, the function calls
\code{badgers_mcmc()}, records the validation (tuning) and test mean squared
prediction errors (MSPE), and keeps the run with the smallest tuning MSPE.
A grid point that fails is skipped; an error is raised only if every grid
point fails.

The outcome is modeled as a linear function of latent covariates
\eqn{x_{ij}}, which are observed through replicates
\eqn{w_{ijr} = x_{ij} + \kappa_{ijr}} with a common measurement error
variance \eqn{\rho^2}. The coefficients have Laplace priors with
predictor-specific rate multipliers \eqn{\lambda_j}, and the log-shrinkage
parameters \eqn{\alpha_j = \log \lambda_j} follow
\eqn{N(\mu 1, \nu \Phi^{-1})}. The precision matrix \eqn{\Phi} has a
constrained Wishart prior whose off-diagonal entries are zero for pairs of
predictors that are not connected in \code{adj_graph}, so that connected
predictors receive related shrinkage. The hyperparameters \code{mu} and
\code{nu} are tuned on the grid.

The defaults of the other arguments are the settings used in the paper:
10000 iterations with 5000 burn-in iterations (also used in the example), and the coarse grid
\code{c(1, 3, 5, 7, 9)} for both \code{mu} and \code{nu}. A fine grid around the
coarse optimum can be searched in a second call.

For large grids, \code{badgers_mcmc()} can be run for each \code{(mu, nu)} pair
in parallel. A complete template that works on Windows, macOS and Linux is
installed with the package and can be found with
\code{system.file("examples", "parallel-runs.R", package = "badgers")}; it is
also described in the parallel-runs vignette.

The returned coefficient summaries \code{CI95}, \code{CI90} and \code{CI80} in
\code{best_result} are posterior means set to zero when the corresponding
95\%, 90\% or 80\% credible interval covers zero. Other levels can be
computed from the coefficient trace \code{best_result$BETA}.
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

fit <- badgers(
  ytrain = sim$y, Wtrain = sim$W,
  ytune = sim$ytune, Wtune = sim$Wtune,
  ytest = sim$ytest, Wtest = sim$Wtest,
  adj_graph = sim$graph,
  mu_grid = c(1, 3),
  nu_grid = c(1, 3),
  n_iter = 10000,
  burnin = 5000,
  verbose = FALSE
)

print(fit)

cat("Test MSPE:", fit$best_result$MSPE_test, "\n")
cat("True nonzero coefficients:", which(sim$beta != 0), "\n")
cat("Selected by 95\% interval:", which(fit$best_result$CI95 != 0), "\n")
cat("Selected by 90\% interval:", which(fit$best_result$CI90 != 0), "\n")
cat("Selected by 80\% interval:", which(fit$best_result$CI80 != 0), "\n")
}
\references{
Hu, Z., Qu, K., Long, Q., Yan, J., and Li, W. A Graph-Guided Bayesian Regression Model for High-Dimensional Data with Replicated Covariates and Measurement Error. Manuscript.
}
