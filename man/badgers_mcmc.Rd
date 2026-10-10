\name{badgers_mcmc}
\alias{badgers_mcmc}
\title{Run the BADGERS MCMC Sampler}
\usage{
badgers_mcmc(ytrain, Wtrain, ytune, Wtune, ytest, Wtest, adj_graph,
  mu, nu, n_iter = 10000, burnin = 5000, k = 50,
  a_sigma = 1, b_sigma = 1, a_iota = 1, b_iota = 1e-2,
  a_rho = 1, b_rho = 1, eta = 10, epsilon = 0.2,
  c_graph = 2, lite = TRUE)
}
\arguments{
\item{ytrain}{Numeric vector of training responses, of length \code{n}.}
\item{Wtrain}{Three-dimensional array of replicate measurements of the training predictors, with dimension \code{num_rep} by \code{n} by \code{p}.}
\item{ytune}{Numeric vector of validation (tuning) responses.}
\item{Wtune}{Three-dimensional array of replicate measurements of the validation predictors.}
\item{ytest}{Numeric vector of test responses.}
\item{Wtest}{Three-dimensional array of replicate measurements of the test predictors.}
\item{adj_graph}{A \code{p} by \code{p} 0/1 adjacency matrix of the working graph among the predictors.}
\item{mu}{Hyperparameter of the graph-informed prior on the shrinkage parameters (prior mean).}
\item{nu}{Hyperparameter of the graph-informed prior on the shrinkage parameters (prior scale).}
\item{n_iter}{Number of MCMC iterations.}
\item{burnin}{Number of burn-in iterations.}
\item{k}{Number of iterations between adaptations of the Metropolis-Hastings proposal standard deviations.}
\item{a_sigma, b_sigma}{Shape and rate hyperparameters of the prior on the residual variance.}
\item{a_iota, b_iota}{Shape and diagonal rate hyperparameters of the prior on the precision matrix penalties.}
\item{a_rho, b_rho}{Shape and rate hyperparameters of the prior on the measurement error variance.}
\item{eta, epsilon}{Tuning constants of the graph weight update.}
\item{c_graph}{Edge-level rate constant of the precision matrix prior; a scalar strictly between 0 and 6.}
\item{lite}{Logical; if \code{TRUE}, full parameter traces are not stored.}
}
\value{
A list containing the validation and test MSPE (\code{MSPE_tune}, \code{MSPE_test}), the posterior mean \code{beta_mean}, the thresholded coefficient estimates \code{CI95} and \code{CI90}, the coefficient trace \code{BETA}, the sampler adaptation summaries, and the run time in minutes (\code{time}). With \code{lite = FALSE} the traces of the latent predictors and the remaining parameters are also returned.
}
\description{
Runs the fixed-\code{(mu, nu)} BADGERS Markov chain Monte Carlo sampler.
}
