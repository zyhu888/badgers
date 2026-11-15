test_that("psi = 6 is the default and other values change the prior rates", {
  args <- fit_args()
  set.seed(81)
  default <- do.call(badgers_mcmc, args)
  set.seed(81)
  explicit <- do.call(badgers_mcmc, c(args, list(psi = 6)))
  expect_identical(default$BETA, explicit$BETA)
  expect_equal(default$psi, 6)
  expect_equal(default$c, 2)
  set.seed(81)
  other <- do.call(badgers_mcmc, c(args, list(psi = 3)))
  expect_equal(other$psi, 3)
  expect_false(identical(default$BETA, other$BETA))
})

test_that("psi_grid is searched with psi outermost and recorded in the fit", {
  args <- fit_args()
  args[c("mu", "nu")] <- NULL
  args$mu_grid <- c(0, 1)
  args$nu_grid <- 1
  args$psi_grid <- c(4, 6)
  args$verbose <- FALSE
  set.seed(91)
  fit <- do.call(badgers, args)
  expect_equal(dim(fit$MSPE_tune_grid), c(2, 1, 2))
  expect_equal(dimnames(fit$MSPE_tune_grid)$psi, c("4", "6"))
  set.seed(91)
  fixed <- args
  fixed[c("mu_grid", "nu_grid", "psi_grid")] <- NULL
  expected <- array(NA_real_, c(2, 1, 2))
  for (h in 1:2) for (a in 1:2) {
    expected[a, 1, h] <- do.call(badgers_mcmc, c(fixed, list(
      mu = args$mu_grid[a], nu = 1, psi = args$psi_grid[h])))$MSPE_tune
  }
  expect_equal(unname(fit$MSPE_tune_grid), expected)
  best <- which(expected == min(expected), arr.ind = TRUE)[1, ]
  expect_equal(fit$best_mu, args$mu_grid[best[1]])
  expect_equal(fit$best_psi, args$psi_grid[best[3]])
  expect_equal(fit$best_result$psi, fit$best_psi)
  expect_output(print(fit), "Selected psi")
})
