# Moore-Penrose pseudoinverse via SVD.
pinv <- function(A, tol = sqrt(.Machine$double.eps)) {
  A <- as.matrix(A)
  s <- svd(A)
  d <- s$d
  if (length(d) == 0) {
    return(matrix(0, ncol(A), nrow(A)))
  }

  threshold <- max(tol * d[1], 0)
  keep <- d > threshold
  if (!any(keep)) {
    return(matrix(0, ncol(A), nrow(A)))
  }

  s$v[, keep, drop = FALSE] %*% (t(s$u[, keep, drop = FALSE]) / d[keep])
}
