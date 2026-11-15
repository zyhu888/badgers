#' @export
print.badgers <- function(x, ...) {
  cat("BADGERS fit\n")
  cat("Selected mu:", x$best_mu, "\n")
  cat("Selected nu:", x$best_nu, "\n")
  cat("Selected psi:", x$best_psi, "\n")
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

  g <- x$MSPE_tune_grid
  if (!is.null(g)) {
    cat("\nTuning MSPE over the (mu, nu, psi) grid:\n")
    psi_vals <- dimnames(g)$psi
    for (h in seq_along(psi_vals)) {
      cat(sprintf("\npsi = %s\n", psi_vals[h]))
      print(matrix(g[, , h], dim(g)[1], dim(g)[2], dimnames = dimnames(g)[1:2]))
    }
  }
  invisible(x)
}
