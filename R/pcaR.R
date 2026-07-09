# Principal component step used to initialise/iterate the MINRES fits.
#
# R -> correlation (or reduced correlation) matrix
# q -> number of components
# B -> component weights (F = Z %*% B)
# A -> loading matrix
#' @noRd
pcaR <- function(R, q) {
  tmp0 <- svd(R)
  Q <- tmp0$v

  Kq <- Q[, 1:q, drop = FALSE]
  Dq <- diag(tmp0$d[1:q], nrow = q)

  B <- Kq %*% sqrt(solve(Dq))
  A <- t(solve(t(B) %*% R %*% B) %*% t(B) %*% R)

  list(B = B, A = A)
}
