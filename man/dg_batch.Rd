\name{dg_batch}
\alias{dg_batch}
\title{Generate Batches of BADGERS Data Sets}
\usage{
dg_batch(R, head, seed, p, g, gsm, ntrain, ntune, ntest, sigma2, rho2,
  num_rep, qmin, qmax, edii, ediu, eduu, Gmode, Gthres, batch = 0,
  verbose = TRUE)
}
\arguments{
\item{R}{Number of data sets to generate.}
\item{head}{Directory in which the data sets are saved; it must exist.}
\item{seed}{Base random seed; data set \code{i} uses \code{seed + batch + i}.}
\item{p, g, gsm, ntrain, ntune, ntest, sigma2, rho2, num_rep, qmin, qmax, edii, ediu, eduu, Gmode, Gthres}{Passed to \code{dg()}.}
\item{batch}{Offset added to the data set index and the seed.}
\item{verbose}{Logical; print progress messages.}
}
\value{
Invisibly \code{NULL}. The data sets are saved as \code{head/data001}, \code{head/data002}, and so on.
}
\description{
Repeatedly calls \code{dg()} and saves each data set to a file.
}
