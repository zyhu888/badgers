test_that("fixed sparsity, singleton dimensions and empty graphs work", {
  for (p in c(1, 4)) for (rep in c(1, 3)) {
    sim <- dg(1, p = p, g = 1, gsm = 0, qmin = 0, qmax = 0,
              num_rep = rep, ntrain = 1, ntune = 1, ntest = 1,
              Gmode = 1, Gthres = 0)
    expect_identical(dim(sim$W), c(as.integer(rep), 1L, as.integer(p)))
    expect_equal(sim$q, 0)
    expect_equal(sum(sim$graph), 0)
    expect_equal(dim(sim$E), c(0L, 2L))
  }
  for (seed in 1:10) {
    sim <- dg(seed, p = 8, g = 1, gsm = 0, qmin = 4, qmax = 4)
    expect_equal(sim$q, 4)
  }
})

test_that("noisy graphs preserve edge counts across boundary cases", {
  for (fraction in c(0, 0.1, 0.75, 1)) for (seed in 1:5) {
    sim <- dg(seed, p = 6, g = 1, gsm = 0, qmin = 6, qmax = 6,
              edii = 1, Gmode = 1, Gthres = fraction)
    expect_equal(sim$nE, sim$nE0)
    expect_equal(sum(sim$graph), sim$nE)
    expect_equal(sim$graph, t(sim$graph))
    expect_equal(diag(sim$graph), rep(0, 6))
    expect_equal(nrow(unique(sim$E)), sim$nE)
    if (fraction == 1) expect_equal(sim$E, sim$E0)
  }
  sim <- small_sim(Gmode = 1, Gthres = 0.1)
  expect_equal(sim$nE, sim$nE0) # one undirected edge, zero retained
})

test_that("simulation covariance and noiseless responses agree analytically", {
  sim <- small_sim(sigma2 = 0, rho2 = 0, num_rep = 1)
  expect_equal(diag(sim$Cov_norm), rep(1, sim$p), tolerance = 1e-12)
  expect_equal(sim$Cov_norm %*% sim$Prec_norm, diag(sim$p), tolerance = 1e-12)
  expect_true(all(eigen(sim$Cov_norm, symmetric = TRUE)$values > 0))
  expect_equal(sim$y, drop(sim$X %*% sim$beta))
  expect_equal(sim$W[1, , ], sim$X)
  expect_identical(small_sim(), small_sim())
})

test_that("invalid simulations fail clearly and batch files round-trip", {
  expect_error(dg(1, qmin = 9, qmax = 8), "qmax")
  expect_error(dg(1, p = 2, qmin = 0, qmax = 0, gsm = 100), "pathway size exceeds")
  expect_error(dg(1, num_rep = 0), "num_rep")
  expect_error(dg(1, edii = -1), "edii")
  expect_error(dg_batch(-1, tempdir(), 1), "R")
  expect_error(dg_batch(1, "nonexistent-badgers-directory", 1), "head")
  expect_error(dg_batch(2, tempdir(), .Machine$integer.max), "seed")
  expect_error(dg_batch(1, tempdir(), -1, batch = .Machine$integer.max), "batch")
  path <- tempfile("badgers-batch-")
  dir.create(path)
  on.exit(unlink(path, recursive = TRUE))
  dg_batch(1, path, 10, p = 4, g = 1, gsm = 0, qmin = 2, qmax = 2,
           verbose = FALSE)
  env <- new.env()
  load(file.path(path, "data001"), envir = env)
  expect_equal(env$data, dg(11, p = 4, g = 1, gsm = 0, qmin = 2, qmax = 2))
})

test_that("simulated latent covariance and measurement noise have the stated moments", {
  N <- 2000L
  sim <- dg(301, p = 4, g = 1, gsm = 0, qmin = 2, qmax = 2,
            ntrain = N, ntune = 2, ntest = 2, num_rep = 2, rho2 = 2)
  C <- sim$Cov_norm
  # Gaussian sample covariance: analytic standard error for each entry.
  covariance_se <- sqrt((C^2 + outer(diag(C), diag(C))) / (N - 1))
  expect_true(all(abs(cov(sim$X) - C) < 6 * covariance_se))
  expect_true(all(abs(colMeans(sim$X)) < 6 / sqrt(N)))
  errors <- sim$W - array(rep(as.vector(sim$X), each = 2), dim(sim$W))
  expect_lt(abs(mean(errors)), 6 * sqrt(2 / length(errors)))
  expect_lt(abs(var(as.vector(errors)) - 2), 6 * 2 * sqrt(2 / (length(errors) - 1)))
})


test_that("negative integer seeds retain R set.seed semantics", {
  sim <- dg(-1, p = 4, g = 1, gsm = 0, qmin = 2, qmax = 2)
  expect_identical(sim, dg(-1, p = 4, g = 1, gsm = 0, qmin = 2, qmax = 2))
})
