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
