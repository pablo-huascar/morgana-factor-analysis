
###############################################
###                                         ###
###      ParallelAnalysisEREC               ###
###                                         ###
###############################################

# Parallel Analysis reference distribution for EREC doublet detection.
#
# Generates a null EREC distribution by Monte Carlo simulation from a
# factor model with no correlated residuals, then returns PA cutoffs for
# each rank position so the caller can assess which observed doublets
# exceed chance.
#
# IN:
#   X           Matrix of observed data (N x m)
#   r           Number of factors
#   cor         1 -> covariance  2 -> Pearson  3 -> polychoric
#   pa_reps     Number of Monte Carlo replications (default 500)
#   percentile  Percentile used as the cutoff (default 95)
#
# OUT:
#   doublets      Matrix (nd x 3): pairs and observed EREC values
#   pa_thresholds Vector (nd): PA cutoff at the requested percentile per rank
#   significant   Logical vector (nd): TRUE when observed EREC > PA cutoff
#   nd            Number of doublet positions evaluated
#   pa_reps       Replications used
#   percentile    Percentile used

ParallelAnalysisEREC <- function(X, r, cor = 2, pa_reps = 500, percentile = 95) {

  N <- nrow(X)
  m <- ncol(X)
  conver <- 0.00001

  # Observed EREC and the number of doublet slots
  obs_result <- LargestDoubletsEREC(X, r, cor)
  nd          <- obs_result$nd
  obs_doublets <- obs_result$doublets

  # Fit baseline model (no doublets) to define the null population
  R <- switch(as.character(cor),
    "1" = cov(X),
    "2" = cor(X),
    "3" = psych::polychoric(X)$rho,
    stop("`cor` must be 1, 2, or 3.")
  )

  if (cor > 1) {
    A_null <- my_minres_hc(R, r, conver, N)$Loadings
  } else {
    A_null <- as.matrix(my_minres(R, r, conver, N)$Loadings)
  }

  # Null correlation matrix: AA' + diag(uniqueness), unit diagonal
  AA     <- A_null %*% t(A_null)
  psi_d  <- pmax(diag(diag(m) - diag(diag(AA))), 0.001)
  Rr     <- AA + diag(psi_d)
  Rr     <- Rr - diag(diag(Rr)) + diag(m)

  # Ensure positive definiteness before Cholesky
  min_ev <- min(eigen(Rr, only.values = TRUE)$values)
  if (min_ev < 0.001) {
    Rr <- Rr + diag(rep(0.001 - min_ev, m))
  }
  L_chol <- chol(Rr)

  # Large pool of null observations for efficient sub-sampling
  NP <- 100000
  Zi  <- matrix(rnorm(NP * m), ncol = m)
  XP  <- Zi %*% L_chol

  # Monte Carlo: collect the nd ranked EREC values from each replication
  pa_erec <- matrix(0, nrow = pa_reps, ncol = nd)
  for (h in seq_len(pa_reps)) {
    idx      <- order(runif(NP))[seq_len(N)]
    Xi       <- XP[idx, ]
    rep_res  <- LargestDoubletsEREC(Xi, r, 2)   # Pearson: null data are continuous
    pa_erec[h, ] <- rep_res$doublets[, 3]
  }

  # Cutoff at the requested percentile, one per rank position
  pa_thresholds <- apply(pa_erec, 2, function(x) quantile(x, percentile / 100))

  list(
    doublets      = obs_doublets,
    pa_thresholds = pa_thresholds,
    significant   = obs_doublets[, 3] > pa_thresholds,
    nd            = nd,
    pa_reps       = pa_reps,
    percentile    = percentile
  )
}
