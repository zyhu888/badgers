test_that("invalid controls fail before compiled sampling", {
  invalid <- list(n_iter = c(0, -1, 1.5, NA, Inf), burnin = c(-1, 12, 20, 0.5),
                  k = c(0, -1, 1.5), nu = c(0, -1, NA, Inf),
                  mu = c(NA, Inf), c = c(0, -1, NA, Inf), psi = c(0, -1, NA, Inf),
                  a_sigma = c(0, -1), b_sigma = 0, a_iota = 0, b_iota_diag = 0,
                  a_rho = 0, b_rho = 0, eta = 0, zeta = 0,
                  lite = NA, verbose = 1)
  for (nm in names(invalid)) for (value in invalid[[nm]]) {
    args <- fit_args()
    args[[nm]] <- value
    expect_error(do.call(badgers_mcmc, args), nm, info = nm)
  }
  args <- fit_args()
  args$mu <- c(1, 2)
  expect_error(do.call(badgers_mcmc, args), "mu")
})

test_that("responses, dimensions, replicate layout and graph are checked", {
  for (value in list(c(NA, 1:4), c(Inf, 1:4), letters[1:5], matrix(1:5))) {
    args <- fit_args()
    args$ytrain <- value
    expect_error(do.call(badgers_mcmc, args), "ytrain")
  }
  for (value in c(Inf, NaN)) {
    args <- fit_args()
    args$Wtrain[1] <- value
    expect_error(do.call(badgers_mcmc, args), "Wtrain")
  }
  args <- fit_args()
  args$Wtrain <- array(numeric(), c(3, 0, 4))
  expect_error(do.call(badgers_mcmc, args), "Wtrain")
  for (mode in c("first", "hole", "inconsistent")) {
    args <- fit_args()
    if (mode == "first") args$Wtrain[1, , 1] <- NA
    if (mode == "hole") args$Wtrain[2, , 1] <- NA
    if (mode == "inconsistent") args$Wtrain[3, 1, 1] <- NA
    expect_error(do.call(badgers_mcmc, args), "replicate|padding|NA", ignore.case = TRUE)
  }
  for (graph in list(matrix(NA, 4, 4), matrix(2, 4, 4), diag(3),
                     matrix(c(0, 1, rep(0, 14)), 4, 4))) {
    args <- fit_args()
    args$adj_graph <- graph
    expect_error(do.call(badgers_mcmc, args), "adj_graph")
  }
})

test_that("named arrays and unequal replicate counts are supported", {
  args <- fit_args()
  for (nm in c("Wtrain", "Wtune", "Wtest")) {
    args[[nm]][3, , 1] <- NA
    dimnames(args[[nm]]) <- lapply(dim(args[[nm]]), function(n) as.character(seq_len(n)))
  }
  args$adj_graph <- args$adj_graph != 0
  expect_true(all(is.finite(do.call(badgers_mcmc, args)$BETA)))
})

test_that("bad grids and shared inputs retain actionable errors", {
  args <- fit_args()
  args[c("mu", "nu")] <- NULL
  args$verbose <- FALSE
  for (grid in list(numeric(), c(1, NA), c(1, Inf))) {
    expect_error(do.call(badgers, c(args, list(mu_grid = grid))), "mu_grid")
  }
  expect_error(do.call(badgers, c(args, list(nu_grid = 0))), "nu_grid")
  for (grid in list(numeric(), c(1, NA), c(0, 1), -2)) {
    expect_error(do.call(badgers, c(args, list(psi_grid = grid))), "psi_grid")
  }
  args$ytrain[1] <- NA
  expect_error(do.call(badgers, args), "ytrain")
})
