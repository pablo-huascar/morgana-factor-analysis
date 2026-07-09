#' MORGANA extraction in confirmatory (mild-restricted) mode
#'
#' Computes EREC values for the doublet pairs supplied by the user
#' ([DoubletsEREC()]) and then fits the MORGANA model with those doublets
#' freed ([Morgana()]).
#'
#' This is a low-level function kept close to the original code released by
#' the authors of the method; most users should call [solutionMCFA()]
#' instead.
#'
#' @inheritParams MorganaEFA
#' @param doblets Two-column matrix of variable index pairs whose residual
#'   correlations should be freed.
#'
#' @return A list with:
#'   \item{doublets}{Matrix (`nd` x 3) of doublets (`P1`, `P2`, `EREC`).}
#'   \item{A}{Loading matrix (m x r), unrotated.}
#'   \item{RESCOV}{Residual variance/covariance matrix.}
#'   \item{RESCOR}{Residual correlation matrix.}
#'
#' @references
#' Lorenzo-Seva, U. (2025). Mild-restricted confirmatory factor analysis: a
#' new approach to assessing population factor models. *Communications in
#' Statistics - Simulation and Computation*.
#' \doi{10.1080/03610918.2025.2453473}
#'
#' @seealso [solutionMCFA()], [MorganaEFA()]
#' @export
MorganaCFA <- function(X, r, conver, maxiter, cor, doblets, R = NULL) {
  X <- as.matrix(X)
  cor_code <- .morgana_cor_code(cor)
  R <- .morgana_dispersion(X, cor_code, R)

  doublets <- DoubletsEREC(X, r, cor_code, doblets, R = R)$doublets
  res <- Morgana(X, r, conver, maxiter, cor_code, doublets, R = R)

  list(doublets = doublets, A = res$L, RESCOV = res$RESCOV, RESCOR = res$RESCOR)
}
