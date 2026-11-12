test_that("Gaussian coefficient conditional has the analytic mean and covariance", {
  X <- matrix(c(1, 0, 1, 2, -1, 1), 3, 2)
  y <- c(2, -1, 0.5)
  tau <- c(0.8, 1.3)
  sigma2 <- 1.7
  V <- sigma2 * solve(crossprod(X) + diag(1 / tau^2))
  m <- drop(V %*% crossprod(X, y) / sigma2)
  N <- 2500L
  set.seed(201)
  draws <- replicate(N, badgers:::update_beta_cpp(X, y, tau, sigma2))
  # Independent normal draws: six analytic Monte Carlo standard errors.
  expect_true(all(abs(rowMeans(draws) - m) < 6 * sqrt(diag(V) / N)))
  covariance_se <- sqrt((V^2 + outer(diag(V), diag(V))) / (N - 1))
  expect_true(all(abs(cov(t(draws)) - V) < 6 * covariance_se))
})

test_that("variance conditionals use inverse-Gamma scale rather than reciprocal rate", {
  X <- matrix(c(1, 2, -1, 0, 1, 2), 3, 2)
  y <- c(1, -1, 2)
  beta <- c(0.5, -0.2)
  tau <- c(1, 2)
  shape <- 2 + (3 + 2) / 2
  scale <- 3 + (sum((y - X %*% beta)^2) + sum(beta^2 / tau^2)) / 2
  set.seed(202)
  expected <- 1 / rgamma(1, shape, rate = scale)
  set.seed(202)
  expect_equal(badgers:::update_sigma2_cpp(y, X, beta, tau, 3, shape), expected)
  W <- array(c(1, 2, 3, 4, 5, 6, 7, 8, NA, NA, NA, NA), c(2, 3, 2))
  # One replicate for predictor 2, with padding consistent across subjects.
  W[2, , 2] <- NA
  W[1, , 2] <- c(7, 8, 9)
  rho_shape <- 1 + 3 * 3 / 2
  rho_scale <- 2
  for (j in 1:2) for (r in seq_len(c(2, 1)[j])) {
    rho_scale <- rho_scale + sum((W[r, , j] - X[, j])^2) / 2
  }
  set.seed(203)
  expected <- 1 / rgamma(1, rho_shape, rate = rho_scale)
  set.seed(203)
  expect_equal(badgers:::update_rho2_cpp(W, X, c(2, 1), 2, rho_shape), expected)
})

test_that("zero coefficients use the exact Gamma local-scale limit", {
  N <- 3000
  set.seed(204)
  tau <- badgers:::update_tau_cpp(rep(0, N), rep(log(2), N), 1)
  expect_true(all(is.finite(tau) & tau > 0))
  # tau^2 ~ Gamma(1/2, rate=2), with mean 1/4 and variance 1/8.
  expect_lt(abs(mean(tau^2) - 1/4), 6 * sqrt(1 / (8 * N)))
})

test_that("precision updates preserve inverse and structural zero constraints", {
  for (adj in list(diag(3), matrix(1, 3, 3),
                   matrix(c(1, 1, 0, 1, 1, 1, 0, 1, 1), 3))) {
    phi <- diag(3) + 0.05 * (adj - diag(3))
    inv <- solve(phi)
    set.seed(205)
    for (i in 1:5) {
      result <- badgers:::update_phi_cpp(phi, inv, c(0.1, -0.3, 0.2), 0, 10, 0.2, 1)
      phi <- result$phi
      inv <- result$inv_phi
      expect_equal(phi, t(phi), tolerance = 1e-12)
      expect_true(all(eigen(phi, symmetric = TRUE)$values > 0))
      expect_equal(phi[adj == 0], rep(0, sum(adj == 0)))
      expect_equal(inv, solve(phi), tolerance = 1e-10)
    }
  }
})

test_that("Metropolis decisions remain valid when density ratios overflow", {
  for (baseline in c(350, 400)) for (seed in 1:20) {
    set.seed(seed)
    proposed <- rnorm(1, baseline, 20)
    u <- runif(1)
    delta <- proposed - baseline
    log_ratio <- 2 * delta - (exp(2 * delta) - 1) / 2 - delta^2 / 2
    expected <- if (log(u) < log_ratio) proposed else baseline
    set.seed(seed)
    actual <- badgers:::update_alpha_cpp(baseline, exp(-baseline), matrix(1), 20, 1, baseline)
    expect_equal(actual$alpha, expected)
  }
})

test_that("latent Gaussian conditionals have the analytic moments", {
  W <- array(c(1, 2, NA, -0.4, 0.3, 0.9), c(3, 1, 2))
  counts <- c(2, 3)
  beta <- c(1, 0.4)
  gamma <- matrix(c(2, 0.4, 0.4, 1), 2)
  rho2 <- 0.8
  sigma2 <- 1.2
  y <- 0.6
  heldout_precision <- gamma + diag(counts / rho2)
  train_precision <- heldout_precision + tcrossprod(beta) / sigma2
  heldout_covariance <- solve(heldout_precision)
  train_covariance <- solve(train_precision)
  measurement_rhs <- colSums(W[, 1, ], na.rm = TRUE) / rho2
  heldout_mean <- drop(heldout_covariance %*% measurement_rhs)
  train_mean <- drop(train_covariance %*% (measurement_rhs + beta * y / sigma2))
  N <- 2500L
  set.seed(206)
  draws <- replicate(N, {
    draw <- badgers:::update_X_cpp(W, W, W, y, beta, counts, gamma, sigma2, rho2)
    c(draw$xtrain[1, ], draw$xtune[1, ], draw$xtest[1, ])
  })
  # Independent Gaussian conditional draws: six analytic Monte Carlo SEs.
  for (block in 0:2) {
    V <- if (block == 0) train_covariance else heldout_covariance
    m <- if (block == 0) train_mean else heldout_mean
    sample <- draws[2 * block + 1:2, , drop = FALSE]
    expect_true(all(abs(rowMeans(sample) - m) < 6 * sqrt(diag(V) / N)))
    covariance_se <- sqrt((V^2 + outer(diag(V), diag(V))) / (N - 1))
    expect_true(all(abs(cov(t(sample)) - V) < 6 * covariance_se))
  }
})

test_that("compiled vector outputs have no matrix dimension attribute", {
  draw <- badgers:::update_beta_cpp(matrix(c(1, 2), 2, 1), c(1, 0), 1, 1)
  expect_null(dim(draw))
  adapted <- badgers:::adapt_alpha_sd_cpp(matrix(c(0, 1), 1), 2L, 2L, 0.8, 0.44, 1)
  expect_null(dim(adapted$alpha_sd_prop))
  expect_null(dim(adapted$r_alpha_col))
})
