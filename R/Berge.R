# Constrained least squares used for Heywood-case rows in MINRES:
# minimise ||f %*% t - phi|| subject to t't = 1.
#
# This is the one-column case of Ten Berge & Nevels (1977), "A general
# solution to Mosier's oblique Procrustes problem", Psychometrika.
#
# Note: the R translation shipped with the original MORGANA code contained
# several MATLAB-to-R indexing faults (e.g. `1:m-q`, scalar `%*%` matrix)
# that made it fail whenever more than one factor was extracted. The
# version below solves the same stationary equations
#   (f'f - b I) t = f'phi,  ||t|| = 1,  b < lambda_min(f'f)
# by direct root finding, which is equivalent to the iterative scheme in
# the paper.
#' @noRd
BERGE <- function(f, PHI) {
  phi <- as.matrix(PHI)[, 1]
  m <- ncol(f)

  ee <- eigen(crossprod(f), symmetric = TRUE)
  lambda <- ee$values
  U <- ee$vectors
  x <- drop(crossprod(U, crossprod(f, phi)))

  lam_min <- lambda[m]
  tol <- max(abs(lambda)) * 1e-10 + 1e-12
  minimal <- (lambda - lam_min) <= tol

  if (sum(x[minimal]^2) > tol) {
    w <- .berge_root(lambda, x, lam_min)
  } else {
    x[minimal] <- 0
    keep <- !minimal
    g_at_min <- sum(x[keep]^2 / (lambda[keep] - lam_min)^2) - 1
    if (g_at_min <= 0) {
      # The unconstrained part does not exhaust the unit norm: place the
      # remainder on the minimal eigenspace (case 3 in the paper).
      w <- numeric(m)
      w[keep] <- x[keep] / (lambda[keep] - lam_min)
      w[which(minimal)[1]] <- sqrt(-g_at_min)
    } else {
      w <- .berge_root(lambda, x, lam_min)
    }
  }

  matrix(U %*% w, ncol = 1)
}

# Root of g(b) = sum(x^2 / (lambda - b)^2) - 1 on (-Inf, lam_min).
#' @noRd
.berge_root <- function(lambda, x, lam_min) {
  g <- function(b) sum(x^2 / (lambda - b)^2) - 1

  hi <- lam_min - 1e-9 * max(1, abs(lam_min))
  while (g(hi) < 0) {
    hi <- lam_min - (lam_min - hi) / 10
  }

  lo <- lam_min - max(1, sqrt(sum(x^2)))
  while (g(lo) > 0) {
    lo <- lam_min - 2 * (lam_min - lo)
  }

  b <- uniroot(g, c(lo, hi), tol = 1e-12)$root
  x / (lambda - b)
}
