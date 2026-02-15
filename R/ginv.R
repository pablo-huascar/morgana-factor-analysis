# Backward-compatible alias used by original MORGANA routines.
ginv <- function(A, tol = sqrt(.Machine$double.eps)) {
  pinv(A, tol = tol)
}
