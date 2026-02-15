# High-level exploratory MORGANA solution.
solutionMEFA <- function(database,
                         r,
                         cor = 2,
                         conver = 0.0001,
                         maxiter = 500,
                         display = 1,
                         montecarlo = 5000) {
  Z <- as.matrix(database)
  if (!is.numeric(Z)) {
    stop("`database` must be a numeric matrix or data.frame.")
  }

  N <- nrow(Z)
  R <- .morgana_dispersion(Z, cor)

  extracted <- MorganaEFA(Z, r = r, conver = conver, maxiter = maxiter, cor = cor)
  A <- posit(extracted$A)
  RESCOR <- extracted$RESCOR
  doubletsEFA <- extracted$doublets
  colnames(doubletsEFA) <- c("P1", "P2", "EREC")

  h0 <- LoseferMorganaH0(A, RESCOR, doubletsEFA[, c(1, 2), drop = FALSE], N, montecarlo)
  fit <- LoseferMorgana(
    R, A, N, RESCOR, nrow(doubletsEFA),
    h0$v, h0$vv, h0$a_n, h0$b1_n, h0$b2_n,
    h0$a_f, h0$b_f, h0$c_f, h0$d_f,
    display = display
  )

  list(
    loadings = A,
    residual_covariance = extracted$RESCOV,
    residual_correlation = RESCOR,
    doublets = doubletsEFA,
    fit = fit,
    h0 = h0
  )
}

# High-level confirmatory MORGANA solution with user-provided doublets.
solutionMCFA <- function(database,
                         r,
                         doublets,
                         cor = 2,
                         conver = 0.0001,
                         maxiter = 500,
                         display = 1,
                         montecarlo = 5000) {
  Z <- as.matrix(database)
  if (!is.numeric(Z)) {
    stop("`database` must be a numeric matrix or data.frame.")
  }

  if (is.null(dim(doublets)) || ncol(as.matrix(doublets)) != 2) {
    stop("`doublets` must be a two-column matrix/data.frame with item index pairs.")
  }

  provided_doublets <- as.matrix(doublets)
  N <- nrow(Z)
  R <- .morgana_dispersion(Z, cor)

  extracted <- MorganaCFA(
    Z,
    r = r,
    conver = conver,
    maxiter = maxiter,
    cor = cor,
    doblets = provided_doublets
  )

  A <- posit(extracted$A)
  RESCOR <- extracted$RESCOR
  doubletsCFA <- extracted$doublets
  colnames(doubletsCFA) <- c("P1", "P2", "EREC")

  h0 <- LoseferMorganaH0(A, RESCOR, doubletsCFA[, c(1, 2), drop = FALSE], N, montecarlo)
  fit <- LoseferMorgana(
    R, A, N, RESCOR, nrow(doubletsCFA),
    h0$v, h0$vv, h0$a_n, h0$b1_n, h0$b2_n,
    h0$a_f, h0$b_f, h0$c_f, h0$d_f,
    display = display
  )

  list(
    loadings = A,
    residual_covariance = extracted$RESCOV,
    residual_correlation = RESCOR,
    doublets = doubletsCFA,
    fit = fit,
    h0 = h0
  )
}

.morgana_dispersion <- function(X, cor) {
  switch(as.character(cor),
    "1" = cov(X),
    "2" = cor(X),
    "3" = psych::polychoric(X)$rho,
    stop("`cor` must be 1 (covariance), 2 (Pearson), or 3 (polychoric).")
  )
}
