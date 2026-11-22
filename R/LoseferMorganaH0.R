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
#'
#' @return A list with the location (`v`), scale (`vv`), normalising
#'   coefficients (`a_n`, `b1_n`, `b2_n`) and Fleishman coefficients
#'   (`a_f`, `b_f`, `c_f`, `d_f`) of the empirical H0 distribution.
#'
#' @inherit LargestDoubletsEREC references
#' @seealso [LoseferMorgana()], [solutionMEFA()]
#' @export
LoseferMorganaH0 <- function(A, RESCOR, doblets, N, montecarlo) {
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

  DCHI <- numeric(montecarlo)
  h <- 1
  cli::cli_progress_bar("LOSEFER H0 distribution", total = montecarlo)
  while (h <= montecarlo) {
    Xi <- XP[sample.int(NP, N), ]
    Ri <- cor(Xi)

    result <- fufa1(Xi, r, doblets, conver, 2, R = Ri)
    Ai <- result$L
    RESCORi <- result$RESCOR

    AAi <- Ai %*% t(Ai)
    PSIi <- sqrt(abs(eye(m) - diag(diag(AAi))))
    RESi <- Ri - (AAi + PSIi %*% (RESCORi - eye(m)) %*% PSIi)

    c2i <- sum(RESi[upper.tri(RESi)]^2) * (N - 1)
    if (c2i > 0.01) {
      DCHI[h] <- c2i
      h <- h + 1
      cli::cli_progress_update()
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
