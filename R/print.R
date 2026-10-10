print.badgers <- function(x, ...) {
  cat("BADGERS fit\n")
  cat("Selected mu:", x$best_mu, "\n")
  cat("Selected nu:", x$best_nu, "\n")
  cat("Tuning MSPE:", x$best_result$MSPE_tune, "\n")
  cat("Test MSPE:", x$best_result$MSPE_test, "\n")
  if (!is.null(x$MSPE_tune_grid)) {
    cat("\nTuning MSPE over the (mu, nu) grid:\n")
    print(x$MSPE_tune_grid)
  }
  invisible(x)
}
