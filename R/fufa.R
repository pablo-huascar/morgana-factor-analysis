# Free Unrestricted Factor Analysis (FUFA), case 1.
#
# Used inside the LOSEFER Monte Carlo simulation: fits the factor model on
# a core set of variables free of doublets, extends the solution to the
# remaining variables, and estimates the residual correlations of the
# doublet pairs.
#
# X       -> data matrix (N x m)
# r       -> number of factors
# doblets -> matrix (p x 2) with pairs of variables with correlated errors
# conver  -> convergence value for MINRES
# cor     -> 1 covariance, 2 Pearson, 3 polychoric
#
# Note: the original translation used `switch(cor, 2, {...}, ...)`, which
# by positional matching only worked for cor = 2; it now handles all codes.
#' @noRd
fufa1 <- function(X, r, doblets, conver, cor, R = NULL) {
  doblets <- as.matrix(doblets)
  N <- nrow(X)
  m <- ncol(X)

  R <- .morgana_dispersion(X, cor, R)

  p <- nrow(doblets)
  Pairs <- c(t(doblets))

  # Core variables: those not involved in any doublet
  Core <- setdiff(1:m, Pairs)
  m_core <- length(Core)
  if (m_core < r + 1) {
    stop("Too many doublet variables: the core set is smaller than r + 1.",
      call. = FALSE
    )
  }

  # Sort the pairs by the average squared partial correlation of the item
  # that will be kept in the core solution.
  PairsSort <- zeros(p, 2)
  rpar <- numeric(p)
  for (h in 1:p) {
    J <- c(
      doblets[h, 1], doblets[h, 2],
      Pairs[Pairs != doblets[h, 1] & Pairs != doblets[h, 2]],
      Core
    )
    Ra <- PartCor(R[J, J])

    r1 <- mean(Ra[3:m, 1]^2)
    r2 <- mean(Ra[3:m, 2]^2)
    if (r1 < r2) {
      PairsSort[h, ] <- J[1:2]
      rpar[h] <- r1
    } else {
      PairsSort[h, ] <- J[2:1]
      rpar[h] <- r2
    }
  }
  PairsSort <- PairsSort[order(rpar), , drop = FALSE]

  # One item of each pair enters the fitted solution (PairIn); its partner
  # is prohibited so that no doublet has both members in the core fit.
  PairIn <- PairsSort[1, 1]
  PairProhibit <- PairsSort[1, 2]
  if (p > 1) {
    for (i in 2:p) {
      s1 <- PairsSort[i, 1]
      s2 <- PairsSort[i, 2]
      if (!(s1 %in% PairIn) && !(s2 %in% PairIn) && !(s1 %in% PairProhibit)) {
        PairIn <- c(PairIn, s1)
        PairProhibit <- c(PairProhibit, s2)
      } else if (!(s2 %in% c(PairIn, PairProhibit))) {
        PairIn <- c(PairIn, s2)
        PairProhibit <- c(PairProhibit, s1)
      }
    }
  }
  PairOut <- setdiff(1:m, c(Core, PairIn))
  m_OUT <- length(PairOut)

  # Step 1: fit the nuclear matrix without doublet partners
  SelCore <- c(Core, PairIn)
  Rcore <- R[SelCore, SelCore]
  CANCORE <- my_minres(Rcore, r, conver, N)$Loadings

  # Step 2: estimate the remaining loadings by extension from the core
  P <- CANCORE[1:m_core, , drop = FALSE]
  PtP <- crossprod(P)
  CANEXT <- zeros(m_OUT, r)
  for (i in seq_len(m_OUT)) {
    vecr <- R[Core, PairOut[i]]
    CANEXT[i, ] <- solve(PtP, crossprod(P, vecr))
  }

  # Step 3: assemble the complete pattern
  L <- rbind(CANCORE, CANEXT)
  L <- L[order(c(SelCore, PairOut)), , drop = FALSE]

  # Residual standard deviations. The original translation divided by
  # sqrt(communality) instead of sqrt(uniqueness); for typical loadings
  # (h2 close to .5) the two nearly coincide, but the uniqueness is the
  # correct scaling for a residual correlation.
  uniq <- pmax(diag(R) - rowSums(L^2), 0.001)
  IPSI <- diag(1 / sqrt(uniq))

  # Step 4: residuals
  RESCOV <- R - (L %*% t(L))
  RESI <- IPSI %*% RESCOV %*% IPSI
  RESCOR <- diag(m)
  for (i in 1:p) {
    RESCOR[doblets[i, 1], doblets[i, 2]] <- RESI[doblets[i, 1], doblets[i, 2]]
    RESCOR[doblets[i, 2], doblets[i, 1]] <- RESI[doblets[i, 2], doblets[i, 1]]
  }

  list(L = L, RESCOV = RESCOV, RESCOR = RESCOR, RESI = RESI)
}
