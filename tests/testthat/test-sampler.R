test_that("summaries and prediction errors agree with independent R calculations", {
  args <- fit_args()
  args$lite <- FALSE
  set.seed(41)
  fit <- do.call(badgers_mcmc, args)
  keep <- seq.int(args$burnin + 1, args$n_iter)
  draws <- fit$BETA[, keep, drop = FALSE]
  expect_equal(fit$beta_mean, rowMeans(draws), tolerance = 1e-12)
  for (level in c(95, 90, 80)) {
    tail <- (1 - level / 100) / 2
    ci <- t(apply(draws, 1, quantile, probs = c(tail, 1 - tail), type = 7))
    expected <- rowMeans(draws)
    expected[ci[, 1] <= 0 & ci[, 2] >= 0] <- 0
    expect_equal(fit[[paste0("CI", level)]], expected, tolerance = 1e-12)
  }
  for (set in c("tune", "test")) {
    pred <- vapply(keep, function(i) {
      drop(fit[[paste0("X", set)]][, , i] %*% fit$BETA[, i])
    }, numeric(length(args[[paste0("y", set)]])))
    expect_equal(fit[[paste0("MSPE_", set)]],
                 mean((args[[paste0("y", set)]] - rowMeans(pred))^2), tolerance = 1e-12)
  }
  args$lite <- TRUE
  set.seed(41)
  light <- do.call(badgers_mcmc, args)
  expect_identical(light$BETA, fit$BETA)
  expect_null(light$Xtrain)
  expect_equal(dim(fit$r_alpha), c(4, 3))
  for (i in seq_len(args$n_iter)) {
    expect_true(all(eigen(fit$GAMMA[, , i], symmetric = TRUE)$values > 0))
  }
})

test_that("held-out outcomes do not leak into fitting", {
  args <- fit_args()
  set.seed(45)
  fit <- do.call(badgers_mcmc, args)
  args$ytune <- args$ytune + 100
  args$ytest <- args$ytest - 100
  set.seed(45)
  other <- do.call(badgers_mcmc, args)
  expect_identical(fit$BETA, other$BETA)
  expect_false(identical(fit$MSPE_tune, other$MSPE_tune))
})

test_that("single predictor and minimal chains keep dimensions", {
  sim <- dg(4, p = 1, g = 1, gsm = 0, qmin = 1, qmax = 1,
            ntrain = 1, ntune = 1, ntest = 1, num_rep = 1)
  args <- fit_args(sim)
  args$n_iter <- 1
  args$burnin <- 0
  args$k <- 2
  fit <- do.call(badgers_mcmc, args)
  expect_equal(dim(fit$BETA), c(1, 1))
  expect_equal(dim(fit$r_alpha), c(1, 0))
  expect_true(is.finite(fit$MSPE_test))
})

test_that("large bipartite graphs initialize with positive definite precision", {
  p <- 42
  adj <- matrix(0, p, p)
  adj[1:21, 22:42] <- adj[22:42, 1:21] <- 1
  # I + .05 * adjacency has minimum eigenvalue -0.05.
  W <- array(0, c(1, 2, p))
  fit <- badgers_mcmc(c(0, 1), W, c(0, 1), W, c(0, 1), W, adj,
                       mu = 1, nu = 1, n_iter = 1, burnin = 0)
  expect_true(all(is.finite(fit$BETA)))
})

test_that("grid selection matches fixed runs and preserves failures", {
  args <- fit_args()
  args[c("mu", "nu")] <- NULL
  args$mu_grid <- c(0, 1)
  args$nu_grid <- c(1, 2)
  args$verbose <- FALSE
  set.seed(71)
  fit <- do.call(badgers, args)
  set.seed(71)
  expected <- matrix(NA_real_, 2, 2)
  fixed <- args
  fixed[c("mu_grid", "nu_grid")] <- NULL
  for (a in 1:2) for (b in 1:2) {
    expected[a, b] <- do.call(badgers_mcmc, c(fixed, list(mu = args$mu_grid[a],
                                                       nu = args$nu_grid[b])))$MSPE_tune
  }
  expect_equal(dim(fit$MSPE_tune_grid), c(2, 2, 1))
  expect_equal(unname(fit$MSPE_tune_grid[, , 1]), expected)
  expect_equal(fit$best_result$MSPE_tune, min(expected))
  expect_equal(fit$best_psi, 6)
  expect_true(all(is.na(fit$grid_errors)))
  original <- badgers_mcmc
  local_mocked_bindings(badgers_mcmc = function(..., mu) {
    if (mu == 0) stop("deliberate numerical failure")
    original(..., mu = mu)
  }, .package = "badgers")
  expect_warning(partial <- do.call(badgers, args), "Some grid points failed")
  expect_true(all(grepl("deliberate", partial$grid_errors[1, , ])))
  args$mu_grid <- 0
  expect_error(do.call(badgers, args), "First error: deliberate")
})
