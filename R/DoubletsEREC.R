#' EREC values for user-specified doublets
#'
#' Computes the EREC index (see [LargestDoubletsEREC()]) for pairs of
#' variables specified by the user, as required by the confirmatory
#' (mild-restricted) MORGANA workflow.
#'
#' This is a low-level function kept close to the original code released by
#' the authors of the method; most users should call [solutionMCFA()]
#' instead.
#'
#' @inheritParams LargestDoubletsEREC
#' @param PairsDoublets Two-column matrix (or data frame) of variable index
#'   pairs whose residual correlations should be freed.
#'
#' @return A list with:
#'   \item{doublets}{Matrix (`nd` x 3) with columns `P1`, `P2`, and `EREC`.}
#'   \item{nd}{Number of pairs retained (pairs beyond the admissible
#'     maximum are dropped with a warning).}
#'   \item{EREC}{The full m x m matrix of EREC values.}
#'
#' @inherit LargestDoubletsEREC references
#' @seealso [solutionMCFA()], [LargestDoubletsEREC()]
#' @export
#' @examples
#' \donttest{
#' database <- as.matrix(read.table(
#'   system.file("extdata", "database.dat", package = "morganafa")
#' ))
#' pairs <- rbind(c(1, 2), c(6, 7))
#' DoubletsEREC(database, r = 2, cor = 2, PairsDoublets = pairs)$doublets
#' }
DoubletsEREC <- function(X, r, cor, PairsDoublets, R = NULL) {
  X <- as.matrix(X)
  N <- nrow(X)
  m <- ncol(X)
  cor_code <- .morgana_cor_code(cor)

  PairsDoublets <- as.matrix(PairsDoublets)
  if (ncol(PairsDoublets) < 2 || nrow(PairsDoublets) < 1) {
    stop("`PairsDoublets` must be a matrix with two columns of variable indices.",
      call. = FALSE
    )
  }
  if (any(PairsDoublets[, 1:2] < 1) || any(PairsDoublets[, 1:2] > m) ||
    any(PairsDoublets[, 1] == PairsDoublets[, 2])) {
    stop("`PairsDoublets` must contain pairs of distinct variable indices between 1 and ", m, ".",
      call. = FALSE
    )
  }

  nd_max <- .erec_max_doublets(m, r)
  nd <- nrow(PairsDoublets)
  if (nd > nd_max) {
    warning(
      "Only ", nd_max, " doublet(s) can be freed with m = ", m,
      " variables and r = ", r, " factors; the first ", nd_max,
      " pair(s) were kept.",
      call. = FALSE
    )
    nd <- nd_max
    PairsDoublets <- PairsDoublets[1:nd, , drop = FALSE]
  }

  R <- .morgana_dispersion(X, cor_code, R)
  EREC <- .erec_matrix(R, r, N, cor_code)

  doublets <- cbind(
    PairsDoublets[, 1],
    PairsDoublets[, 2],
    abs(EREC[PairsDoublets[, 1:2, drop = FALSE]])
  )
  colnames(doublets) <- c("P1", "P2", "EREC")

  list(doublets = doublets, nd = nd, EREC = EREC)
}
