# MINRES extraction by Zegers & ten Berge (1983), row-wise updates.
#
# C    -> covariance or correlation matrix
# q    -> number of factors
# conv -> convergence value
# N    -> number of subjects
#
# `heywood = TRUE` applies the Ten Berge & Nevels (1977) constrained update
# whenever a row communality would exceed 1 (Heywood case).
#' @noRd
.minres_engine <- function(C, q, conv, N, heywood) {
  A <- pcaR(C, q)$A

  n <- nrow(A)
  Co <- C - diag(diag(C))

  loss <- function(A) {
    AA <- A %*% t(A)
    ssq(Co - AA - diag(diag(AA)))
  }

  g <- loss(A)
  gold <- g * 2
  iter <- 0

  while (abs(gold - g) > conv) {
    # Ai is A with row i zeroed. Instead of copying A for every row, keep a
    # single working copy in step with A (zero row i before the update,
    # write the new row back after). The arithmetic is unchanged, so the
    # results are bit-for-bit identical to the row-by-row copy; calling
    # solve.default() directly also skips S3 dispatch.
    Ai <- A
    for (i in 1:n) {
      Ai[i, ] <- 0
      ci <- Co[, i, drop = FALSE]
      ai <- solve.default(crossprod(Ai), crossprod(Ai, ci))

      if (heywood && sum(ai^2) > 1) {
        ai <- BERGE(Ai, ci)
      }

      A[i, ] <- ai
      Ai[i, ] <- ai
    }
    gold <- g
    g <- loss(A)
    iter <- iter + 1
    if (iter > 1000) {
      gold <- g
    }
  }

  Re <- (A %*% t(A)) + diag(diag(C)) - diag(diag(A %*% t(A)))
  Um <- (N - 1) * (log(det(Re)) - log(det(C)))
  if (Um < 0) {
    Um <- 0
  }
  list("Loadings" = A, "Index of Fit" = Um)
}

#' @noRd
my_minres <- function(C, q, conv, N) {
  .minres_engine(C, q, conv, N, heywood = FALSE)
}
