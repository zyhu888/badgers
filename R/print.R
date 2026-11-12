#' @export
print.badgers <- function(x, ...) {
  cat("BADGERS fit\n")
  cat("Selected mu:", x$best_mu, "\n")
  cat("Selected nu:", x$best_nu, "\n")
  cat("Tuning MSPE:", x$best_result$MSPE_tune, "\n")
  cat("Test MSPE:", x$best_result$MSPE_test, "\n")

  ci <- x$best_result[c("CI95", "CI90", "CI80")]
  ci <- ci[!vapply(ci, is.null, logical(1))]
  if (length(ci) > 0) {
    cat("\nNonzero coefficients by credible interval:\n")
    for (nm in names(ci)) {
      cat(sprintf("  %s%%: %d\n", substring(nm, 3), sum(ci[[nm]] != 0)))
    }
  }

  if (!is.null(x$MSPE_tune_grid)) {
    cat("\nTuning MSPE over the (mu, nu) grid:\n")
    print(x$MSPE_tune_grid)
  }
  invisible(x)
}
