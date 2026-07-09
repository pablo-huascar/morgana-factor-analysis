#' LOSEFER-corrected fit indices for a MORGANA solution
#'
#' Computes the LOSEFER-adjusted chi-square and the derived RMSEA, CFI and
#' NNFI/TLI, using the empirical H0 distribution parameters estimated by
#' [LoseferMorganaH0()].
#'
#' This is a low-level function kept close to the original code released by
#' the authors of the method; most users should call [solutionMEFA()] or
#' [solutionMCFA()] instead.
#'
#' @param S Dispersion (correlation/covariance) matrix that was analysed.
#' @param A Loading matrix (m x r), unrotated.
#' @param N Number of observations.
#' @param RESCOR Residual correlation matrix returned by the MORGANA fit.
#' @param nd Number of freed doublets.
#' @param mchi,sxchi Location/scale of the cube-root H0 chi-square values
#'   (`v` and `vv` from [LoseferMorganaH0()]).
#' @param a_n,b1_n,b2_n Normalising polynomial coefficients from
#'   [LoseferMorganaH0()].
#' @param a_f,b_f,c_f,d_f Fleishman coefficients from [LoseferMorganaH0()].
#' @param display If `1`, prints a short summary of the fit indices.
#'
#' @return A list with `Um` (adjusted chi-square), `d` (degrees of freedom),
#'   `d0` and `C0` (null model), `RMSEA`, `CFI`, and `TL` (NNFI/TLI).
#'
#' @inherit LargestDoubletsEREC references
#' @seealso [LoseferMorganaH0()], [solutionMEFA()]
#' @export
LoseferMorgana <- function(S, A, N, RESCOR, nd, mchi, sxchi,
                           a_n, b1_n, b2_n, a_f, b_f, c_f, d_f,
                           display = 0) {
  m <- nrow(A)
  r <- ncol(A)

  # Residuals of the fitted model
  AA <- A %*% t(A)
  PSI2 <- eye(m) - diag(diag(AA))
  PSI <- sqrt(PSI2)
  RESCORZ <- RESCOR - eye(m)
  RES <- S - (AA + PSI %*% RESCORZ %*% PSI)

  FF <- sum(RES[upper.tri(RES)]^2)
  c2 <- FF * (N - 1)
  d <- ((m - r) * (m - r + 1)) / 2 - m - nd

  # Map the observed statistic through the empirical H0 distribution
  z <- (c2^(1 / 3) - mchi) / sxchi
  zn <- a_n + (b1_n * z) + (b2_n * z * z)
  Um <- (a_f + zn * (b_f + zn * (c_f + zn * d_f))) * sqrt(2 * d) + d

  if (Um > d + 300 * sqrt(d * 2)) {
    Um <- d + 300 * sqrt(d * 2) + runif(1) * mchi * 4 - (mchi * 2)
  }

  # Null (independence) model
  C0 <- sum(S[upper.tri(S)]^2) * (N - 1)
  d0 <- m * (m - 1) / 2

  CFI <- abs(((C0 - d0) - (Um - d)) / (C0 - d0))
  if (CFI > .999) {
    CFI <- .999
  } else if (CFI < 0) {
    CFI <- 0.01
  }

  Mo <- C0 / d0
  Mk <- Um / d
  TL <- (Mo - Mk) / (Mo - 1)
  if (TL > .999) {
    TL <- .999
  } else if (TL < 0) {
    TL <- CFI
  }

  FFo <- Um / (N - 1) - (d / (N - 1))
  if (FFo < 0) {
    FFo <- 0
  }
  RMSEA <- sqrt(FFo / d)

  if (display == 1) {
    cat(sprintf(
      "Chi Square adjusted by LOSEFER with %i degrees of freedom = %.3f\n",
      as.integer(d), Um
    ))
    cat(sprintf("Goodness of Fit (Null Model) = %.2f\n", C0))
    cat(sprintf("Degrees of freedom = %.0f\n", d0))
    cat(sprintf("RMSEA (Root Mean Square Error of Approximation) = %.3f\n", RMSEA))
    cat(sprintf("CFI (Comparative Fit Index) = %.3f\n", CFI))
    cat(sprintf("NNFI (Tucker-Lewis Index) = %.3f\n", TL))
  }

  list(Um = Um, d = d, d0 = d0, C0 = C0, RMSEA = RMSEA, CFI = CFI, TL = TL)
}
