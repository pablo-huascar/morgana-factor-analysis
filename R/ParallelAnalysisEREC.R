#' Parallel Analysis reference thresholds for the EREC index
#'
#' Generates a null EREC distribution by Monte Carlo simulation from a
#' factor model with no correlated residuals, and returns Parallel Analysis
#' cutoffs for each rank position so observed doublets can be compared
#' against chance (see the "Reference values and upper bounds" section of
#' Ferrando, Hernandez-Dorado & Lorenzo-Seva, 2022).
#'
#' A doublet is flagged as *significant* when its observed EREC exceeds the
#' PA cutoff at the same rank position.
#'
#' Each replication requires m*(m-1)/2 model fits, so this can be slow for
#' large m; consider `pa_reps = 100` for initial explorations.
#'
#' @inheritParams LargestDoubletsEREC
#' @param pa_reps Number of Monte Carlo replications (default 500).
#' @param percentile Percentile used as the cutoff (default 95).
#'
#' @return An object of class `morgana_pa`: a list with
#'   \item{doublets}{Matrix (`nd` x 3) of observed doublets (`P1`, `P2`,
#'     `EREC`).}
#'   \item{pa_thresholds}{PA cutoff at the requested percentile for each
#'     rank position.}
#'   \item{significant}{Logical: observed EREC > cutoff, per rank.}
#'   \item{nd, pa_reps, percentile}{Settings used.}
#'   \item{item_names}{Column names of `X`.}
#'
#' @inherit LargestDoubletsEREC references
#' @seealso [solutionMEFA()], [LargestDoubletsEREC()]
#' @export
#' @examples
#' \donttest{
#' database <- as.matrix(read.table(
#'   system.file("extdata", "database.dat", package = "morganafa")
#' ))
#' set.seed(1)
#' pa <- ParallelAnalysisEREC(database, r = 2, cor = 2, pa_reps = 50)
#' pa
#' }
ParallelAnalysisEREC <- function(X, r, cor = 2, pa_reps = 500,
                                 percentile = 95, R = NULL) {
  X <- .morgana_check_data(X, arg = "X")
  N <- nrow(X)
  m <- ncol(X)
  conver <- 0.00001
  cor_code <- .morgana_cor_code(cor)

  R <- .morgana_dispersion(X, cor_code, R)

  # Observed EREC and the number of doublet slots
  obs_result <- LargestDoubletsEREC(X, r, cor_code, R = R)
  nd <- obs_result$nd
  obs_doublets <- obs_result$doublets

  # Fit the baseline model (no doublets) to define the null population
  if (cor_code > 1) {
    A_null <- as.matrix(my_minres_hc(R, r, conver, N)$Loadings)
  } else {
    A_null <- as.matrix(my_minres(R, r, conver, N)$Loadings)
  }

  # Null correlation matrix: AA' + diag(uniqueness), unit diagonal
  AA <- A_null %*% t(A_null)
  Rr <- AA - diag(diag(AA)) + diag(m)

  # Ensure positive definiteness before the Cholesky factorisation
  min_ev <- min(eigen(Rr, symmetric = TRUE, only.values = TRUE)$values)
  if (min_ev < 0.001) {
    Rr <- Rr + diag(rep(0.001 - min_ev, m))
  }
  L_chol <- chol(Rr)

  # Large pool of null observations, sub-sampled at each replication
  NP <- 100000
  XP <- matrix(rnorm(NP * m), ncol = m) %*% L_chol

  # Monte Carlo: collect the nd ranked EREC values from each replication.
  # The null data are continuous, so Pearson correlations are used.
  pa_erec <- matrix(0, nrow = pa_reps, ncol = nd)
  for (h in seq_len(pa_reps)) {
    Xi <- XP[sample.int(NP, N), ]
    pa_erec[h, ] <- LargestDoubletsEREC(Xi, r, 2)$doublets[, 3]
  }

  pa_thresholds <- apply(pa_erec, 2, quantile, probs = percentile / 100)

  structure(
    list(
      doublets = obs_doublets,
      pa_thresholds = unname(pa_thresholds),
      significant = unname(obs_doublets[, 3] > pa_thresholds),
      nd = nd,
      pa_reps = pa_reps,
      percentile = percentile,
      item_names = colnames(X)
    ),
    class = "morgana_pa"
  )
}
