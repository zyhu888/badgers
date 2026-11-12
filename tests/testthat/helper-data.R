small_sim <- function(...) {
  dg(seed = 17, p = 4, g = 1, gsm = 0, qmin = 2, qmax = 2,
     ntrain = 5, ntune = 3, ntest = 2, ...)
}

fit_args <- function(sim = small_sim()) {
  list(ytrain = sim$y, Wtrain = sim$W, ytune = sim$ytune,
       Wtune = sim$Wtune, ytest = sim$ytest, Wtest = sim$Wtest,
       adj_graph = sim$graph, mu = 1, nu = 1, n_iter = 12,
       burnin = 6, k = 4)
}
