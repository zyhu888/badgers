\name{dg}
\alias{dg}
\title{Generate BADGERS Data}
\usage{
dg(seed, p, g, gsm, ntrain, ntune, ntest, sigma2, rho2, num_rep,
  qmin, qmax, edii, ediu, eduu, Gmode, Gthres, savefile = NULL)
}
\arguments{
\item{seed}{Random seed.}
\item{p}{Number of predictors.}
\item{g}{Number of pathways.}
\item{gsm}{Mean number of predictors per pathway.}
\item{ntrain}{Training sample size.}
\item{ntune}{Validation (tuning) sample size.}
\item{ntest}{Test sample size.}
\item{sigma2}{Residual variance.}
\item{rho2}{Variance of the measurement error in the replicates.}
\item{num_rep}{Number of replicate measurements per predictor.}
\item{qmin}{Minimum number of nonzero coefficients.}
\item{qmax}{Maximum number of nonzero coefficients.}
\item{edii}{Edge density between important predictors in the true graph.}
\item{ediu}{Edge density between important and unimportant predictors.}
\item{eduu}{Edge density between unimportant predictors.}
\item{Gmode}{\code{0} uses the true graph as the working graph; any other value uses a noisy working graph.}
\item{Gthres}{Fraction of true edges kept in the noisy working graph, in (0, 1]; the remaining edges are random so that the number of edges is unchanged.}
\item{savefile}{Optional file name; if given, the generated list is saved there with \code{save()}.}
}
\value{
A list containing the true and working graphs (\code{E0}, \code{E}, \code{graph}), the coefficients \code{beta}, and the response and replicate arrays for the training, tuning and test sets (\code{y}, \code{W}, \code{ytune}, \code{Wtune}, \code{ytest}, \code{Wtest}).
}
\description{
Generates a true graph among predictors, multivariate normal predictors with the corresponding precision structure, replicate measurements with measurement error, and a sparse linear response. A working graph, either equal to the true graph or contaminated with random edges, is returned for the sampler.
}
