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
    lite = TRUE
) {

  start.time <- Sys.time()

  out <- badgers_mcmc_cpp(
    ytrain = ytrain, Wtrain_ = Wtrain, ytune = ytune, Wtune_ = Wtune,
    ytest = ytest, Wtest_ = Wtest, adj_graph = adj_graph,
    n_iter = n_iter, burnin = burnin, k = k,
    a_sigma = a_sigma, b_sigma = b_sigma, a_iota = a_iota, b_iota = b_iota,
    a_rho = a_rho, b_rho = b_rho, eta = eta, epsilon = epsilon,
    mu = mu, nu = nu, c_graph = c_graph, lite = lite
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
        c_graph = c_graph, lite = lite
      ),
      error = function(e) {
        if (verbose) cat("ERROR:", conditionMessage(e), "\n")
        NULL
      }
    )

    if (is.null(out)) next

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

  if (is.null(best_result)) stop("All grid points failed.")

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
      MSPE_test_grid = MSPE_test_grid
    ),
    class = "badgers"
  )
}


dg_batch <- function(
    R,
    head,
    seed,
    p,
    g,
    gsm,
    ntrain,
    ntune,
    ntest,
    sigma2,
    rho2,
    num_rep,
    qmin,
    qmax,
    edii,
    ediu,
    eduu,
    Gmode,
    Gthres,
    batch = 0,
    verbose = TRUE
) {

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
    p,
    g,
    gsm,
    ntrain,
    ntune,
    ntest,
    sigma2,
    rho2,
    num_rep,
    qmin,
    qmax,
    edii,
    ediu,
    eduu,
    Gmode,
    Gthres,
    savefile = NULL
) {

  set.seed(seed)

  ## ---- Number of nonzero coefficients, uniform on qmin:qmax ----
  q <- sample(qmin:qmax, 1)

  n <- ntrain
  pathway <- matrix(FALSE, p, g)
  ngene <- rpois(g, gsm)
  ngene[1] <- max(ngene[1], q)
  pathway[1:ngene[1], 1] <- TRUE
  for (i in 2:g) {
    pathway[sample(1:p, ngene[i]), i] <- TRUE
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
  A[1:q, 1:q] <- abs(A[1:q, 1:q])

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
  y <- drop(X %*% beta + sqrt(sigma2) * rnorm(n))

  Xtune <- matrix(rnorm(ntune * p), ntune) %*% t(chol_A_t_norm)
  Wtune <- apply(Xtune, c(1, 2), function(x) {
    x + rnorm(num_rep, 0, sd = sqrt(rho2))
  })
  ytune <- drop(Xtune %*% beta + sqrt(sigma2) * rnorm(ntune))

  Xtest <- matrix(rnorm(ntest * p), ntest) %*% t(chol_A_t_norm)
  Wtest <- apply(Xtest, c(1, 2), function(x) {
    x + rnorm(num_rep, 0, sd = sqrt(rho2))
  })
  ytest <- drop(Xtest %*% beta + sqrt(sigma2) * rnorm(ntest))

  ## ---- Working graph G given to the sampler ----
  if (Gmode == 0) {
    E <- E0
  } else {
    ## Keep the top Gthres fraction of G0 edges ranked by |partial
    ## correlation|, then add random edges so that |G| = |G0|.
    uni_index <- which(E0[, 1] < E0[, 2])
    E_tmp <- E0[uni_index, ]
    pcor_tmp <- pcor[uni_index]
    odr_tmp <- order(abs(pcor_tmp), decreasing = TRUE)
    num_remain <- round(Gthres * nrow(E_tmp))
    incl <- odr_tmp[1:num_remain]
    E <- E_tmp[incl, ]
    E <- rbind(E, cbind(E[, 2], E[, 1]))
    E1 <- E[order(E[, 1], E[, 2]), ]
    e <- nE0 / 2

    A_tilde <- Matrix::rsparsematrix(
      p, p, nnz = e + p, symmetric = TRUE, rand.x = function(n) return(1)
    )
    A_tilde[lower.tri(A_tilde, diag = TRUE)] <- 0
    E2 <- which(A_tilde == 1, arr.ind = TRUE)

    ## Drop random edges that are already in E1 (and any duplicates).
    key2 <- E2[, 1] + p * (E2[, 2] - 1)
    key1 <- E1[, 1] + p * (E1[, 2] - 1)
    E3 <- E2[!(key2 %in% key1) & !duplicated(key2), , drop = FALSE]
    E3 <- E3[sample(1:nrow(E3), e - nrow(E1) / 2), ]

    E <- rbind(E1, E3, cbind(E3[, 2], E3[, 1]))
    E <- E[order(E[, 1], E[, 2]), ]
  }
  nE <- nrow(E)

  graph <- matrix(0, nrow = p, ncol = p)
  for (e in 1:nrow(E)) {
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
