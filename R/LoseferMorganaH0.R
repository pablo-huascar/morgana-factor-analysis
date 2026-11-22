#' Empirical H0 distribution for the LOSEFER correction
#'
#' Simulates the null distribution of the chi-square statistic under the
#' fitted MORGANA model (population without sampling error) via Monte Carlo,
#' and returns the parameters used by [LoseferMorgana()] to correct the
#' observed statistic.
#'
#' This is a low-level function kept close to the original code released by
#' the authors of the method; most users should call [solutionMEFA()] or
#' [solutionMCFA()] instead.
#'
#' @param A Loading matrix (m x r), unrotated.
#' @param RESCOR Residual correlation matrix returned by the MORGANA fit.
#' @param doblets Two-column matrix with the freed doublet pairs.
#' @param N Number of observations.
#' @param montecarlo Number of Monte Carlo replications (e.g. 5000). The
#'   running time grows linearly with it; a progress bar is shown in
#'   interactive sessions.
#' @param cores Number of parallel worker processes for the Monte Carlo
#'   replications. `1` (default) runs sequentially. With more cores the
#'   results are identical to the sequential run with the same seed (see
#'   [ParallelAnalysisEREC()] for details). The default can be set once per
#'   session with `options(morganafa.cores = 8)`.
#'
#' @return A list with the location (`v`), scale (`vv`), normalising
#'   coefficients (`a_n`, `b1_n`, `b2_n`) and Fleishman coefficients
#'   (`a_f`, `b_f`, `c_f`, `d_f`) of the empirical H0 distribution.
#'
#' @inherit LargestDoubletsEREC references
#' @seealso [LoseferMorgana()], [solutionMEFA()]
#' @export
LoseferMorganaH0 <- function(A, RESCOR, doblets, N, montecarlo,
                             cores = getOption("morganafa.cores", 1L)) {
  cl <- .morgana_cluster(cores)
  if (!is.null(cl)) {
    on.exit(parallel::stopCluster(cl), add = TRUE)
  }
  .losefer_h0(A, RESCOR, doblets, N, montecarlo, cl)
}

# Workhorse of LoseferMorganaH0(); `cl` is a PSOCK cluster or NULL
# (sequential), so solutionMEFA() can share one cluster across steps.
#' @noRd
.losefer_h0 <- function(A, RESCOR, doblets, N, montecarlo, cl = NULL) {
  conver <- 0.0001

  m <- nrow(A)
  r <- ncol(A)
  doblets <- as.matrix(doblets)
  nd <- nrow(doblets)

  # Model-implied population correlation matrix
  AA <- A %*% t(A)
  PSI2 <- diag(m) - diag(diag(AA))
  PSI <- sqrt(PSI2)
  RESCORZ <- RESCOR - diag(diag(RESCOR))
  Rr <- AA + PSI %*% RESCORZ %*% PSI
  Rr <- Rr - diag(diag(Rr)) + diag(m)
  df <- ((m - r) * (m - r + 1)) / 2 - m - nd

  # Ensure positive definiteness before the Cholesky factorisation
  min_ev <- min(eigen(Rr, symmetric = TRUE, only.values = TRUE)$values)
  if (min_ev < 0.001) {
    Rr <- Rr + diag(rep(0.001 - min_ev, m))
  }
  L <- chol(Rr)

  # Large pool of null observations, sub-sampled at each replication
  NP <- 100000
  XP <- matrix(rnorm(NP * m), ncol = m) %*% L

  cli::cli_progress_bar("LOSEFER H0 distribution", total = montecarlo)
  if (is.null(cl)) {
    DCHI <- numeric(montecarlo)
    h <- 1
    while (h <= montecarlo) {
      Xi <- XP[sample.int(NP, N), ]
      c2i <- .h0_c2(Xi, r, doblets, conver)
      if (c2i > 0.01) {
        DCHI[h] <- c2i
        h <- h + 1
        cli::cli_progress_update()
      }
    }
  } else {
    # Draw exactly as many samples as are still needed, in the sequential
    # order. Every one of them would also be drawn by the sequential loop,
    # so the accepted statistics and the final RNG state are the same.
    parallel::clusterCall(cl, .set_pool, XP)
    args <- list(r = r, doblets = doblets, conver = conver)
    DCHI <- numeric(0)
    while (length(DCHI) < montecarlo) {
      need <- montecarlo - length(DCHI)
      idx <- lapply(seq_len(need), function(k) sample.int(NP, N))
      c2 <- unlist(.pool_map(cl, idx, .h0_worker, args))
      DCHI <- c(DCHI, c2[c2 > 0.01])
      cli::cli_progress_update(set = length(DCHI))
    }
  }
  cli::cli_progress_done()

  # Cube-root transform, trimming extreme replicates
  md <- median(DCHI)
  sx <- sd(DCHI)
  x <- DCHI[DCHI < (md + 5.5 * sx)]^(1 / 3)

  res <- Trans2Normal(x)

  skew <- sqrt(8 / df)
  kurt <- 12 / df
  coeffs <- fit_fleishman_from_sk(skew, kurt)

  list(
    v = res$md, vv = res$sx,
    a_n = res$a, b1_n = res$b1, b2_n = res$b2,
    a_f = -coeffs$B, b_f = coeffs$A, c_f = coeffs$B, d_f = coeffs$C
  )
}
