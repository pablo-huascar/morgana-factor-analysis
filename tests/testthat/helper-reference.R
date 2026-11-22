# Frozen copy of the MINRES engine as of morganafa 0.1.1, before the
# performance changes. Used to check that optimisations leave the results
# bit-for-bit identical.
reference_minres_engine <- function(C, q, conv, N, heywood) {
  A <- morganafa:::pcaR(C, q)$A
  n <- nrow(A)
  Co <- C - diag(diag(C))
  loss <- function(A) {
    AA <- A %*% t(A)
    sum((Co - AA - diag(diag(AA)))^2)
  }
  g <- loss(A)
  gold <- g * 2
  iter <- 0
  while (abs(gold - g) > conv) {
    for (i in 1:n) {
      Ai <- A
      Ai[i, ] <- 0
      ci <- Co[, i, drop = FALSE]
      ai <- solve(crossprod(Ai), crossprod(Ai, ci))
      if (heywood && sum(ai^2) > 1) {
        ai <- morganafa:::BERGE(Ai, ci)
      }
      A[i, ] <- ai
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
