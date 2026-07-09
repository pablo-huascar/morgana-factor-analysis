# Shared workhorse for the EREC index.
#
# For every pair (i, j) the model is refitted on the m - 2 remaining
# variables (the "core set"), the loadings of i and j are recovered by
# extension analysis, and the residual covariance of the pair is rescaled
# by the smallest admissible uniquenesses across pairs. The result is the
# full m x m matrix of EREC values (Expected REsidual Correlation direct
# change; EPC1 in Ferrando, Hernandez-Dorado & Lorenzo-Seva, 2022).
#
# R        -> dispersion matrix (correlation/covariance), computed once
# r        -> number of factors
# N        -> number of observations
# cor_code -> 1 covariance, 2 Pearson, 3 polychoric (selects the Heywood-
#             protected MINRES for correlation input, as in the original)
#' @noRd
.erec_matrix <- function(R, r, N, cor_code) {
  m <- ncol(R)
  conver <- 0.00001
  fit_fun <- if (cor_code > 1) my_minres_hc else my_minres

  total_pairs <- m * (m - 1) / 2
  RESFREE <- eye(m)
  IPSIs <- matrix(0, nrow = total_pairs, ncol = m)

  h <- 0
  for (i in 1:(m - 1)) {
    for (j in (i + 1):m) {
      h <- h + 1
      Core <- setdiff(1:m, c(i, j))

      Ac <- as.matrix(fit_fun(R[Core, Core], r, conver, N)$Loadings)

      # Extension analysis: loadings of the left-out pair
      AtA <- crossprod(Ac)
      lam <- matrix(0, nrow = 2, ncol = r)
      lam[1, ] <- solve(AtA, crossprod(Ac, R[Core, i]))
      lam[2, ] <- solve(AtA, crossprod(Ac, R[Core, j]))

      # Uniquenesses implied by this pair's solution (for correcting EPC1)
      Li <- rbind(Ac, lam)[order(c(Core, i, j)), , drop = FALSE]
      IPSIs[h, ] <- (1 - rowSums(Li^2))^(-0.5)

      # Residual covariance of the pair under the core-set model
      RESFREE[i, j] <- R[i, j] - crossprod(lam[1, ], lam[2, ])
      RESFREE[j, i] <- RESFREE[i, j]
    }
  }

  # Guard against Heywood rows: non-finite or near-zero uniqueness weights
  # are replaced by the largest finite weight (the original code replaced
  # values below 0.001 by the maximum).
  bad <- !is.finite(IPSIs) | IPSIs < 0.001
  if (any(bad)) {
    IPSIs[bad] <- max(IPSIs[!bad], 1)
  }

  IPSI <- diag(colmin(IPSIs))
  IPSI %*% RESFREE %*% IPSI
}

# Maximum number of doublets that can be freed, given the Ledermann bound
# and the model degrees of freedom. Returns the number of available slots
# (possibly 0 or negative when the model is too large for m variables).
#' @noRd
.erec_max_doublets <- function(m, r) {
  Ledermann <- floor((2 * m + 1 - sqrt(8 * m + 1)) / 2)
  df <- ((m - r) * (m - r + 1)) / 2 - m
  nd <- Ledermann - r
  if ((df - nd) < 1) {
    nd <- df - 1
  }
  nd
}
