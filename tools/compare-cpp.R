# Run each package in a separate R process to avoid sharing native symbols.
# Rscript tools/compare-cpp.R record /path/to/baseline-lib /tmp/baseline.rds
# Rscript tools/compare-cpp.R record /path/to/candidate-lib /tmp/candidate.rds
# Rscript tools/compare-cpp.R compare /tmp/baseline.rds /tmp/candidate.rds /tmp/comparison.csv
# Timing is excluded; all other returned fields and final RNG states are checked.

record_cpp_results <- function(lib, output) {
  Sys.setenv(TZ = "UTC")
  library("badgers", lib.loc = lib, character.only = TRUE)
  internal <- function(name) getFromNamespace(name, "badgers")
  results <- list()
  capture <- function(name, seed, expression, scope = "block") {
    set.seed(seed)
    value <- force(expression)
    results[[name]] <<- list(scope = scope, value = value, rng = .Random.seed)
    cat("Recorded", name, "\n")
    flush.console()
  }
  for (p in c(1L, 3L, 8L, 30L)) {
    set.seed(100 + p)
    n <- 6L
    X <- matrix(rnorm(n * p), n, p)
    y <- rnorm(n)
    beta <- seq_len(p) / p
    tau <- seq_len(p) / p + 0.5
    alpha <- seq_len(p) / (2 * p)
    S <- crossprod(X)
    gamma <- diag(p) + matrix(0.05, p, p)
    inv_gamma <- solve(gamma)
    upper <- which(upper.tri(gamma))
    ij <- arrayInd(upper, dim(gamma))
    lower <- ij[, 2] + p * (ij[, 1] - 1)
    not_j <- vapply(seq_len(p), function(j) setdiff(seq_len(p), j), integer(p - 1L))
    dim(not_j) <- c(p - 1L, p)
    W <- array(rnorm(3 * n * p), c(3, n, p))
    W[3, , 1] <- NA
    counts <- c(2, rep(3, p - 1L))
    prefix <- paste0("block/p", p, "/")
    capture(paste0(prefix, "sigma2"), 201,
            internal("update_sigma2_cpp")(y, X, beta, tau, 1, 1 + (n + p) / 2))
    capture(paste0(prefix, "beta"), 202,
            internal("update_beta_cpp")(X, y, tau, 1.7))
    capture(paste0(prefix, "tau"), 203,
            internal("update_tau_cpp")(beta, alpha, 1.7))
    capture(paste0(prefix, "tau_zero"), 203,
            internal("update_tau_cpp")(rep(0, p), alpha, 1.7))
    capture(paste0(prefix, "rho2"), 204,
            internal("update_rho2_cpp")(W, X, counts, 1, 1 + n * sum(counts) / 2))
    capture(paste0(prefix, "latent_X"), 205,
            internal("update_X_cpp")(W, W, W, y, beta, counts, gamma, 1.7, 2))
    for (graph in c("empty", "path", "complete")) {
      adj <- diag(p)
      if (graph == "complete") adj[] <- 1
      if (graph == "path" && p > 1L) {
        for (j in seq_len(p - 1L)) adj[j, j + 1L] <- adj[j + 1L, j] <- 1
      }
      phi <- diag(p) + 0.05 * (adj - diag(p))
      capture(paste0(prefix, "phi_", graph), 206,
              internal("update_phi_cpp")(phi, solve(phi), alpha, 0, 10, 0.2, 1))
      capture(paste0(prefix, "alpha_", graph), 207,
              internal("update_alpha_cpp")(alpha, tau, phi, rep(1, p), 1, 0))
      capture(paste0(prefix, "density_", graph), 207,
              internal("p_density_mmh_cpp")(alpha[1] + 0.2, alpha[1], tau[1],
                                             alpha, phi, 1L, 1, 0))
    }
    capture(paste0(prefix, "gamma"), 208,
            internal("update_gamma_cpp")(X, n, p, gamma, inv_gamma, 1,
                                           rep(1e-4, length(upper)), 0.01,
                                           upper, lower, not_j))
    acs <- matrix(rep(c(0, 1, 1, 0, 1), length.out = p * 15L), p, 15L)
    capture(paste0(prefix, "adaptation"), 209,
            internal("adapt_alpha_sd_cpp")(acs, 15L, 5L, 0.8, 0.44, rep(1, p)))
  }
  for (p in c(1L, 4L, 8L, 30L)) for (seed in c(17L, 41L, 83L)) {
    n <- if (p == 1L) 1L else 12L
    sim <- dg(seed, p = p, g = 1, gsm = 0, qmin = min(2L, p),
              qmax = min(2L, p), ntrain = n, ntune = n, ntest = n)
    for (graph in c("empty", "path", "complete")) {
      adj <- matrix(0, p, p)
      if (graph == "complete") adj[] <- 1
      if (graph == "path" && p > 1L) {
        for (j in seq_len(p - 1L)) adj[j, j + 1L] <- adj[j + 1L, j] <- 1
      }
      sim$W[3, , 1] <- sim$Wtune[3, , 1] <- sim$Wtest[3, , 1] <- NA
      n_iter <- if (p == 30L) 500L else 1000L
      args <- list(ytrain = sim$y, Wtrain = sim$W, ytune = sim$ytune,
                   Wtune = sim$Wtune, ytest = sim$ytest, Wtest = sim$Wtest,
                   adj_graph = adj, mu = 1, nu = 1, n_iter = n_iter,
                   burnin = n_iter %/% 2, k = 25L, lite = p <= 4L)
      name <- sprintf("chain/p%d/seed%d/%s", p, seed, graph)
      capture(name, seed + 1000L, do.call(badgers_mcmc, args), "chain")
    }
  }
  # Match the default 10000/5000 sampling schedule at manageable dimensions.
  for (p in c(8L, 30L)) {
    sim <- dg(41, p = p, g = 1, gsm = 0, qmin = 2, qmax = 2,
              ntrain = 12, ntune = 8, ntest = 6)
    adj <- matrix(1, p, p)
    capture(paste0("chain/long_p", p), 1041,
            badgers_mcmc(sim$y, sim$W, sim$ytune, sim$Wtune, sim$ytest, sim$Wtest,
                         adj, mu = 3, nu = 3, n_iter = 10000, burnin = 5000,
                         lite = p > 8L), "chain")
  }
  sim <- dg(321, p = 42, g = 1, gsm = 0, qmin = 2, qmax = 2,
            ntrain = 12, ntune = 8, ntest = 6)
  adj <- matrix(0, 42, 42)
  adj[1:21, 22:42] <- adj[22:42, 1:21] <- 1
  capture("chain/bipartite42", 1321,
          badgers_mcmc(sim$y, sim$W, sim$ytune, sim$Wtune, sim$ytest, sim$Wtest,
                       adj, mu = 3, nu = 3, n_iter = 50, burnin = 25, lite = FALSE), "chain")
  sim <- dg(123)
  capture("chain/default_p200", 1123,
          badgers_mcmc(sim$y, sim$W, sim$ytune, sim$Wtune, sim$ytest, sim$Wtest,
                       sim$graph, mu = 3, nu = 3, n_iter = 5, burnin = 2), "chain")
  sim <- dg(123, p = 8, g = 1, gsm = 0, qmin = 2, qmax = 2,
            ntrain = 12, ntune = 8, ntest = 6)
  capture("grid/p8", 999,
          badgers(sim$y, sim$W, sim$ytune, sim$Wtune, sim$ytest, sim$Wtest,
                  sim$graph, mu_grid = c(0, 1), nu_grid = c(1, 2),
                  n_iter = 200, burnin = 100, verbose = FALSE), "chain")
  saveRDS(list(results = results, session = sessionInfo()), output)
  cat("Recorded", length(results), "cases in", output, "\n")
}

compare_cpp_results <- function(baseline, candidate, output) {
  old <- readRDS(baseline)$results
  new <- readRDS(candidate)$results
  if (!identical(names(old), names(new))) stop("Comparison cases differ.")
  rows <- list()
  walk <- function(x, y, path, tolerance, exact = FALSE) {
    same <- identical(x, y)
    if (!identical(attributes(x), attributes(y)) || !identical(typeof(x), typeof(y))) {
      rows[[length(rows) + 1L]] <<- data.frame(
        field = path, identical = FALSE, max_abs = NA_real_, max_scaled = NA_real_,
        tolerance = tolerance, passed = FALSE)
    } else if (is.list(x)) {
      for (i in seq_along(x)) {
        label <- names(x)[i]
        if (is.null(label)) label <- as.character(i)
        if (label == "time") next
        walk(x[[i]], y[[i]], paste(path, label, sep = "/"), tolerance, exact)
      }
    } else {
      abs_error <- scaled_error <- 0
      numerical <- is.numeric(x) && !exact
      if (numerical) {
        if (length(x) != length(y) || any(!is.finite(x)) || any(!is.finite(y))) {
          abs_error <- scaled_error <- Inf
        } else if (!same) {
          abs_error <- max(abs(x - y), 0)
          scaled_error <- max(abs(x - y) / pmax(1, abs(x), abs(y)), 0)
        }
      }
      passed <- if (numerical) scaled_error <= tolerance else same
      rows[[length(rows) + 1L]] <<- data.frame(
        field = path, identical = same, max_abs = abs_error,
        max_scaled = scaled_error, tolerance = if (numerical) tolerance else 0,
        passed = passed)
    }
  }
  for (name in names(old)) {
    tolerance <- if (old[[name]]$scope == "block") 1e-12 else 1e-8
    walk(old[[name]]$value, new[[name]]$value, name, tolerance)
    walk(old[[name]]$rng, new[[name]]$rng, paste0(name, "/RNG"), 0, TRUE)
    x <- old[[name]]$value
    y <- new[[name]]$value
    if (old[[name]]$scope == "chain") {
      if (inherits(x, "badgers")) {
        walk(c(x$best_mu, x$best_nu), c(y$best_mu, y$best_nu),
             paste0(name, "/selected_grid"), 0, TRUE)
        x <- x$best_result
        y <- y$best_result
      }
      walk(x$acs_alpha, y$acs_alpha, paste0(name, "/acceptance_exact"), 0, TRUE)
      for (level in c("CI95", "CI90", "CI80")) {
        walk(x[[level]] != 0, y[[level]] != 0,
             paste0(name, "/selection_", level), 0, TRUE)
      }
    }
  }
  report <- do.call(rbind, rows)
  write.csv(report, output, row.names = FALSE)
  cat(sum(report$passed), "of", nrow(report), "field comparisons passed;",
      sum(report$identical), "were bit-for-bit identical.\n")
  failures <- report[!report$passed, ]
  if (nrow(failures)) print(head(failures, 30), row.names = FALSE)
  if (nrow(failures)) quit(status = 1L)
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 3L && args[1] == "record") {
  record_cpp_results(args[2], args[3])
} else if (length(args) == 4L && args[1] == "compare") {
  compare_cpp_results(args[2], args[3], args[4])
} else {
  stop("Use: record LIB OUTPUT.rds, or compare BASELINE.rds CANDIDATE.rds REPORT.csv")
}
