#' Detect the largest doublets via the EREC index
#'
#' Computes the EREC index (Expected REsidual Correlation direct change;
#' EPC1 in Ferrando, Hernandez-Dorado & Lorenzo-Seva, 2022) for every pair
#' of variables and returns the pairs with the largest absolute values,
#' limited by the Ledermann bound and the available degrees of freedom.
#'
#' This is a low-level function kept close to the original code released by
#' the authors of the method; most users should call [solutionMEFA()]
#' instead.
#'
#' @param X Numeric data matrix (N observations x m variables).
#' @param r Number of factors.
#' @param cor Dispersion matrix to analyse: 1 or `"covariance"`, 2 or
#'   `"pearson"`, 3 or `"polychoric"`.
#' @param R Optional pre-computed dispersion matrix. When supplied, `cor` is
#'   only used to select the extraction variant and `R` is not recomputed
#'   from `X` (useful to avoid repeating an expensive polychoric estimation).
#'
#' @return A list with:
#'   \item{doublets}{Matrix (`nd` x 3) with columns `P1`, `P2` (variable
#'     indices) and `EREC` (absolute EREC value), sorted by decreasing EREC.}
#'   \item{nd}{Number of doublet slots evaluated.}
#'   \item{EREC}{The full m x m matrix of EREC values.}
#'
#' @references
#' Ferrando, P. J., Hernandez-Dorado, A., & Lorenzo-Seva, U. (2022).
#' Detecting correlated residuals in exploratory factor analysis: New
#' proposals and a comparison of procedures. *Structural Equation Modeling*,
#' 29(4), 630-638. \doi{10.1080/10705511.2021.2004543}
#'
#' @seealso [solutionMEFA()], [DoubletsEREC()], [ParallelAnalysisEREC()]
#' @export
#' @examples
#' \donttest{
#' database <- as.matrix(read.table(
#'   system.file("extdata", "database.dat", package = "morganafa")
#' ))
#' LargestDoubletsEREC(database, r = 2, cor = 2)$doublets
#' }
LargestDoubletsEREC <- function(X, r, cor, R = NULL) {
  X <- as.matrix(X)
  N <- nrow(X)
  m <- ncol(X)
  cor_code <- .morgana_cor_code(cor)

  nd <- .erec_max_doublets(m, r)
  if (nd < 1) {
    stop(
      "No doublet can be freed with m = ", m, " variables and r = ", r,
      " factors: the Ledermann bound / degrees of freedom leave no room.",
      call. = FALSE
    )
  }

  R <- .morgana_dispersion(X, cor_code, R)
  EREC <- .erec_matrix(R, r, N, cor_code)

  ut <- which(upper.tri(EREC), arr.ind = TRUE)
  erec_vals <- abs(EREC[ut])
  J <- order(erec_vals, decreasing = TRUE)[1:nd]

  doublets <- cbind(ut[J, 1], ut[J, 2], erec_vals[J])
  colnames(doublets) <- c("P1", "P2", "EREC")

  list(doublets = doublets, nd = nd, EREC = EREC)
}
