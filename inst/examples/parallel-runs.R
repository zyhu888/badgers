## Parallel BADGERS runs
##
## This example mirrors the small runnable example in ?badgers, but evaluates
## the (psi, mu, nu) grid in parallel by calling badgers_mcmc() directly. It
## then refines (mu, nu) around the coarse optimum at the selected psi, as in
## the paper workflow.
##
## It uses a PSOCK cluster from the parallel package, which works on Windows,
## macOS and Linux. It uses 1000 iterations with 500 burn-in to keep the run
## short; the runs are spread over the available cores. For real analyses
## use the defaults of badgers_mcmc(), which are the paper settings (10000
## iterations, 5000 burn-in).

## Optional: uncomment if you use a custom R library on a cluster.
## .libPaths("~/rlibs")

library(badgers)
library(parallel)

seed <- 123

n <- 40
p <- 30
burnin <- 500
n_iter <- 1000

sim <- dg(
  seed = seed,
  p = p,
  g = 5,
  gsm = 8,
  ntrain = n,
  ntune = n,
  ntest = n,
  sigma2 = 2,
  rho2 = 2,
  num_rep = 3,
  qmin = 4,
  qmax = 6,
  edii = 0.6,
  ediu = 0,
  eduu = 0.3,
  Gmode = 0,
  Gthres = 1
)

out_dir <- file.path(tempdir(), "badgers_parallel_example")
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

## Workers do not share memory with this session, so run_one() takes
## everything it needs as arguments.
run_one <- function(i, grid, tag, sim, out_dir, n_iter, burnin) {
  psi_i <- grid$psi[i]
  mu_i <- grid$mu[i]
  nu_i <- grid$nu[i]

  tryCatch({
    out_file <- file.path(
      out_dir, sprintf("%s_psi%g_mu%g_nu%g.rds", tag, psi_i, mu_i, nu_i)
    )

    if (file.exists(out_file)) {
      message(sprintf("psi %g mu %g nu %g already exists, skipping: %s",
                      psi_i, mu_i, nu_i, out_file))
      return(NULL)
    }

    fit <- badgers::badgers_mcmc(
      ytrain = sim$y, Wtrain = sim$W,
      ytune = sim$ytune, Wtune = sim$Wtune,
      ytest = sim$ytest, Wtest = sim$Wtest,
      adj_graph = sim$graph,
      mu = mu_i,
      nu = nu_i,
      psi = psi_i,
      n_iter = n_iter,
      burnin = burnin
    )

    tmp_file <- paste0(out_file, ".tmp")
    saveRDS(
      list(
        psi = psi_i, mu = mu_i, nu = nu_i,
        MSPE_tune = fit$MSPE_tune, MSPE_test = fit$MSPE_test,
        CI95 = fit$CI95, CI90 = fit$CI90, CI80 = fit$CI80,
        beta_mean = fit$beta_mean
      ),
      file = tmp_file
    )
    file.rename(tmp_file, out_file)

    invisible(NULL)
  }, error = function(e) {
    message(sprintf("ERROR: psi %g mu %g nu %g failed: %s",
                    psi_i, mu_i, nu_i, e$message))
    invisible(NULL)
  })
}

## ---- Start the cluster ----
n_cores <- detectCores()
if (is.na(n_cores)) n_cores <- 2
n_workers <- max(1, min(n_cores - 1, 25))

cl <- makeCluster(n_workers)

## Use the same library paths as this session, load badgers on the workers,
## and give each worker its own reproducible random number stream.
invisible(clusterCall(cl, function(paths) .libPaths(paths), .libPaths()))
clusterEvalQ(cl, library(badgers))
clusterSetRNGStream(cl, seed)

run_grid <- function(grid, tag) {
  parLapplyLB(
    cl,
    seq_len(nrow(grid)),
    run_one,
    grid = grid, tag = tag, sim = sim, out_dir = out_dir,
    n_iter = n_iter, burnin = burnin
  )

  res <- lapply(seq_len(nrow(grid)), function(i) {
    f <- file.path(out_dir, sprintf("%s_psi%g_mu%g_nu%g.rds",
                                    tag, grid$psi[i], grid$mu[i], grid$nu[i]))
    if (file.exists(f)) readRDS(f) else NULL
  })
  res[!vapply(res, is.null, logical(1))]
}

summarize <- function(res) {
  data.frame(
    psi = vapply(res, `[[`, numeric(1), "psi"),
    mu = vapply(res, `[[`, numeric(1), "mu"),
    nu = vapply(res, `[[`, numeric(1), "nu"),
    MSPE_tune = vapply(res, `[[`, numeric(1), "MSPE_tune"),
    MSPE_test = vapply(res, `[[`, numeric(1), "MSPE_test")
  )
}

## ---- Step 1: coarse grid over (psi, mu, nu) ----
coarse_grid <- expand.grid(psi = c(4, 6), mu = c(1, 3, 5), nu = c(1, 3, 5))
coarse_res <- run_grid(coarse_grid, "coarse")
coarse_table <- summarize(coarse_res)
coarse_best <- coarse_table[which.min(coarse_table$MSPE_tune), ]

## ---- Step 2: fine (mu, nu) grid centered at the coarse optimum ----
make_fine_values <- function(center) {
  seq(from = max(0.5, center - 1), to = center + 1, by = 0.5)
}

fine_grid <- expand.grid(
  psi = coarse_best$psi,
  mu = make_fine_values(coarse_best$mu),
  nu = make_fine_values(coarse_best$nu)
)
fine_res <- run_grid(fine_grid, "fine")
fine_table <- summarize(fine_res)

stopCluster(cl)

best_idx <- which.min(fine_table$MSPE_tune)
best <- fine_res[[best_idx]]

cat("Selected mu:", best$mu, "\n")
cat("Selected nu:", best$nu, "\n")
cat("Selected psi:", best$psi, "\n")
cat("Test MSPE:", best$MSPE_test, "\n")
cat("True nonzero coefficients:", which(sim$beta != 0), "\n")
cat("Selected by 95% interval:", which(best$CI95 != 0), "\n")
cat("Selected by 90% interval:", which(best$CI90 != 0), "\n")
cat("Selected by 80% interval:", which(best$CI80 != 0), "\n")
