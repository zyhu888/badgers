\name{dg_batch}
\alias{dg_batch}
\title{Generate Batches of BADGERS Data Sets}
\usage{
dg_batch(R, head, seed, p = 200, g = 20, gsm = 10, ntrain = 50,
  ntune = 50, ntest = 50, sigma2 = 2, rho2 = 2, num_rep = 3,
  qmin = 4, qmax = 8, edii = 0.6, ediu = 0, eduu = 0.3, Gmode = 0,
  Gthres = 1, batch = 0, verbose = TRUE)
}
\arguments{
\item{R}{Number of data sets to generate.}
\item{head}{Directory in which the data sets are saved; it must exist.}
\item{seed}{Base random seed; data set \code{i} uses \code{seed + batch + i}.}
\item{p, g, gsm, ntrain, ntune, ntest, sigma2, rho2, num_rep, qmin, qmax, edii, ediu, eduu, Gmode, Gthres}{Passed to \code{dg()}; see there for the defaults.}
\item{batch}{Offset added to the data set index and the seed.}
\item{verbose}{Logical; print progress messages.}
}
\value{
Invisibly \code{NULL}. The data sets are saved as \code{head/data001}, \code{head/data002}, and so on, each holding one list named \code{data}.
}
\description{
Repeatedly calls \code{dg()} and saves each data set to a file.
}
\examples{
out_dir <- file.path(tempdir(), "badgers_data")
dir.create(out_dir)

dg_batch(
  R = 2, head = out_dir, seed = 100,
  p = 30, g = 5, gsm = 8, ntrain = 20, ntune = 20, ntest = 20,
  verbose = FALSE
)

load(file.path(out_dir, "data001"))
names(data)[1:6]
}
