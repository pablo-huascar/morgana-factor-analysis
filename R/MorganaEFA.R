#' MORGANA extraction in exploratory mode
#'
#' Detects the largest doublets with [LargestDoubletsEREC()] and then fits
#' the MORGANA model with those doublets freed ([Morgana()]).
#'
#' This is a low-level function kept close to the original code released by
#' the authors of the method; most users should call [solutionMEFA()]
#' instead.
#'
#' @inheritParams Morgana
#' @param cor Dispersion matrix to analyse: 1 or `"covariance"`, 2 or
#'   `"pearson"`, 3 or `"polychoric"`.
#'
#' @return A list with:
#'   \item{doublets}{Matrix (`nd` x 3) of detected doublets (`P1`, `P2`,
#'     `EREC`).}
#'   \item{A}{Loading matrix (m x r), unrotated.}
#'   \item{RESCOV}{Residual variance/covariance matrix.}
#'   \item{RESCOR}{Residual correlation matrix.}
#'
#' @inherit LargestDoubletsEREC references
#' @seealso [solutionMEFA()], [MorganaCFA()]
#' @export
MorganaEFA <- function(X, r, conver, maxiter, cor, R = NULL) {
  X <- as.matrix(X)
  cor_code <- .morgana_cor_code(cor)
  R <- .morgana_dispersion(X, cor_code, R)

  doublets <- LargestDoubletsEREC(X, r, cor_code, R = R)$doublets
  res <- Morgana(X, r, conver, maxiter, cor_code, doublets, R = R)

  list(doublets = doublets, A = res$L, RESCOV = res$RESCOV, RESCOR = res$RESCOR)
}
