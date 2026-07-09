#' MORGANA estimation with fixed doublets
#'
#' Core MORGANA routine: fits an unrestricted factor model in which the
#' residual correlations of the supplied doublet pairs are fixed at their
#' EREC estimates (signed by the corresponding partial correlations), using
#' iterated eigendecomposition of the doubly-reduced matrix.
#'
#' This is a low-level function kept close to the original code released by
#' the authors of the method; most users should call [solutionMEFA()] or
#' [solutionMCFA()] instead.
#'
#' @param X Numeric data matrix (N observations x m variables).
#' @param r Number of factors.
#' @param conver Convergence criterion for the iterative fit.
#' @param maxiter Maximum number of iterations.
#' @param corr Dispersion matrix to analyse: 1 or `"covariance"`, 2 or
#'   `"pearson"`, 3 or `"polychoric"`.
#' @param doblets Matrix (`nd` x 3) of doublets: columns are the two
#'   variable indices and the (absolute) EREC value, as returned by
#'   [LargestDoubletsEREC()] or [DoubletsEREC()].
#' @param R Optional pre-computed dispersion matrix (see
#'   [LargestDoubletsEREC()]).
#'
#' @return A list with:
#'   \item{L}{Loading matrix (m x r), unrotated.}
#'   \item{RESCOV}{Residual variance/covariance matrix.}
#'   \item{RESCOR}{Residual correlation matrix (identity except for the
#'     freed doublets).}
#'
#' @inherit LargestDoubletsEREC references
#' @seealso [solutionMEFA()], [solutionMCFA()]
#' @export
Morgana <- function(X, r, conver, maxiter, corr, doblets, R = NULL) {
  X <- as.matrix(X)
  R <- .morgana_dispersion(X, corr, R)
  m <- ncol(R)

  doblets <- as.matrix(doblets)
  if (ncol(doblets) < 3 || nrow(doblets) < 1) {
    stop("`doblets` must be a matrix with at least one row and columns P1, P2, EREC.",
      call. = FALSE
    )
  }
  d <- nrow(doblets)

  # Matrix of partial correlations: AI is Guttman's anti-image matrix and
  # P the matrix of partial correlations.
  Rinv <- solve(R)
  D <- diag(diag(Rinv))
  tmp0 <- diag(D^(-1 / 2))
  AI <- diag(tmp0) %*% Rinv %*% diag(tmp0)
  P <- (2 * eye(m)) - AI

  # Target residual correlations: EREC magnitude, signed by the partial
  # correlation of the pair.
  Pr <- eye(m)
  for (i in 1:d) {
    Pr[doblets[i, 1], doblets[i, 2]] <-
      doblets[i, 3] * sign(P[doblets[i, 1], doblets[i, 2]])
    Pr[doblets[i, 2], doblets[i, 1]] <- Pr[doblets[i, 1], doblets[i, 2]]
  }

  # Initial estimated solution (no iterations). RESI includes PSI^2 in the
  # diagonal and the correlated residuals off-diagonal.
  U <- diag(diag(D^(-1 / 2)))
  RED <- R - U %*% Pr %*% U
  L <- pcaR(RED, r)$A

  f <- sum(diag(t(L) %*% L))
  fold <- f * 10
  iter <- 1

  while (abs(f - fold) > conver) {
    W <- R - (L %*% t(L))
    U <- sqrt(diag(diag(W)))
    RED <- R - (U %*% Pr %*% U)
    L <- pcaR(RED, r)$A
    iter <- iter + 1
    fold <- f
    if (iter < maxiter) {
      f <- sum(diag(t(L) %*% L))
    }
  }

  RESCOV <- R - (L %*% t(L))
  RESI <- solve(U) %*% RESCOV %*% solve(U)
  RESCOR <- eye(m)
  for (i in 1:d) {
    RESCOR[doblets[i, 1], doblets[i, 2]] <- RESI[doblets[i, 1], doblets[i, 2]]
    RESCOR[doblets[i, 2], doblets[i, 1]] <- RESI[doblets[i, 2], doblets[i, 1]]
  }

  list(L = L, RESCOV = RESCOV, RESCOR = RESCOR)
}
