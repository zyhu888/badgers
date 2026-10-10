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
\item{Wtrain}{Three-dimensional array of replicate measurements of the training predictors, with dimension \code{num_rep} by \code{n} by \code{p}.}
\item{ytune}{Numeric vector of validation (tuning) responses.}
\item{Wtune}{Three-dimensional array of replicate measurements of the validation predictors, with dimension \code{num_rep} by \code{ntune} by \code{p}.}
\item{ytest}{Numeric vector of test responses.}
\item{Wtest}{Three-dimensional array of replicate measurements of the test predictors, with dimension \code{num_rep} by \code{ntest} by \code{p}.}
\item{adj_graph}{A \code{p} by \code{p} 0/1 adjacency matrix of the working graph among the predictors.}
\item{mu_grid}{Candidate values of \code{mu}.}
\item{nu_grid}{Candidate values of \code{nu}.}
\item{n_iter}{Number of MCMC iterations.}
\item{burnin}{Number of burn-in iterations.}
\item{k}{Number of iterations between adaptations of the Metropolis-Hastings proposal standard deviations.}
\item{a_sigma, b_sigma}{Shape and rate hyperparameters of the prior on the residual variance.}
\item{a_iota, b_iota}{Shape and diagonal rate hyperparameters of the prior on the precision matrix penalties.}
\item{a_rho, b_rho}{Shape and rate hyperparameters of the prior on the measurement error variance.}
\item{eta, epsilon}{Tuning constants of the graph weight update.}
\item{c_graph}{Edge-level rate constant of the precision matrix prior; a scalar strictly between 0 and 6.}
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
\code{badgers()} is the primary interface for fitting BADGERS. For each
combination of \code{mu_grid} and \code{nu_grid}, the function calls
\code{badgers_mcmc()}, records the validation (tuning) and test mean squared
prediction errors (MSPE), and keeps the run with the smallest tuning MSPE.
A grid point that fails is skipped; an error is raised only if every grid
point fails.

The predictors are observed through replicates: each of \code{Wtrain},
\code{Wtune} and \code{Wtest} is a \code{num_rep} by \code{n} by \code{p}
array of noisy measurements of the latent predictors. The working graph
\code{adj_graph} encodes known or partially known relationships among the
predictors and enters the prior on the predictor-specific shrinkage
parameters. The hyperparameters \code{mu} and \code{nu} of this prior are
tuned on the grid.

The returned coefficient summaries \code{CI95} and \code{CI90} in
\code{best_result} are posterior means set to zero when the corresponding
credible interval covers zero.
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
  n_iter = 200,
  burnin = 100,
  verbose = FALSE
)

print(fit)

fit$best_result$MSPE_test
which(fit$best_result$CI95 != 0)
which(sim$beta != 0)
}
\references{
XXX. References will be added after the manuscript is available.
}
