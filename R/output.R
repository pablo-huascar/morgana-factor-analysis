# Construction of the user-facing result objects and their print methods.

# Assemble a morgana_efa / morgana_cfa object from the pieces produced by
# the low-level routines.
#' @noRd
.morgana_solution <- function(class, A, extracted, doublets, fit, h0, pa,
                              rotate, item_names, labels, settings) {
  r <- ncol(A)
  factor_names <- paste0("F", seq_len(r))

  rot <- .morgana_rotate(A, rotate)

  dimnames(A) <- list(item_names, factor_names)
  dimnames(rot$loadings) <- list(item_names, factor_names)
  if (!is.null(rot$Phi)) {
    dimnames(rot$Phi) <- list(factor_names, factor_names)
  }

  RESCOR <- extracted$RESCOR
  RESCOV <- extracted$RESCOV
  dimnames(RESCOR) <- list(item_names, item_names)
  dimnames(RESCOV) <- list(item_names, item_names)

  h2 <- rowSums(A^2)
  u2 <- diag(RESCOV)

  loadings_tbl <- .morgana_loadings_tbl(rot$loadings, item_names, labels, h2, u2)
  loadings_unrot_tbl <- .morgana_loadings_tbl(A, item_names, labels, h2, u2)
  doublets_tbl <- .morgana_doublets_tbl(doublets, item_names, labels, RESCOR, pa)
  fit_tbl <- .morgana_fit_tbl(fit)

  structure(
    list(
      loadings = loadings_tbl,
      Phi = rot$Phi,
      loadings_unrotated = loadings_unrot_tbl,
      doublets = doublets_tbl,
      fit = fit_tbl,
      residual_correlation = RESCOR,
      residual_covariance = RESCOV,
      loading_matrix = rot$loadings,
      pa = pa,
      h0 = h0,
      settings = settings
    ),
    class = c(class, "morgana_solution")
  )
}

# Rotate an (already sign-reflected) unrotated loading matrix.
#' @noRd
.morgana_rotate <- function(A, rotate) {
  rotate <- match.arg(tolower(rotate[1]), c("oblimin", "varimax", "none"))
  r <- ncol(A)
  if (rotate == "none" || r == 1) {
    return(list(loadings = A, Phi = NULL, rotate = "none"))
  }

  if (rotate == "oblimin") {
    res <- GPArotation::oblimin(A)
    P <- unclass(res$loadings)
    Phi <- res$Phi
  } else {
    res <- GPArotation::Varimax(A)
    P <- unclass(res$loadings)
    Phi <- diag(r)
  }

  # The sign of each rotated factor is arbitrary: reflect so that the sum
  # of loadings is positive, adjusting Phi accordingly.
  signs <- ifelse(colSums(P) < 0, -1, 1)
  P <- sweep(P, 2, signs, `*`)
  Phi <- diag(signs, r) %*% Phi %*% diag(signs, r)

  list(loadings = P, Phi = Phi, rotate = rotate)
}

#' @noRd
.morgana_loadings_tbl <- function(L, item_names, labels, h2, u2) {
  out <- tibble::tibble(item = item_names)
  if (!is.null(labels)) {
    out$label <- labels
  }
  for (k in seq_len(ncol(L))) {
    out[[colnames(L)[k]]] <- unname(L[, k])
  }
  out$h2 <- unname(h2)
  out$u2 <- unname(u2)
  out
}

#' @noRd
.morgana_doublets_tbl <- function(doublets, item_names, labels, RESCOR, pa) {
  i1 <- doublets[, 1]
  i2 <- doublets[, 2]
  out <- tibble::tibble(
    rank = seq_len(nrow(doublets)),
    item_a = item_names[i1],
    item_b = item_names[i2]
  )
  if (!is.null(labels)) {
    out$label_a <- labels[i1]
    out$label_b <- labels[i2]
  }
  out$erec <- unname(doublets[, 3])
  out$residual_r <- RESCOR[cbind(i1, i2)]
  if (!is.null(pa)) {
    out$pa_threshold <- pa$pa_thresholds[seq_len(nrow(out))]
    out$significant <- pa$significant[seq_len(nrow(out))]
  }
  out
}

#' @noRd
.morgana_fit_tbl <- function(fit) {
  tibble::tibble(
    measure = c(
      "chi2_losefer", "df", "chi2_null", "df_null",
      "rmsea", "cfi", "nnfi_tli"
    ),
    value = c(fit$Um, fit$d, fit$C0, fit$d0, fit$RMSEA, fit$CFI, fit$TL)
  )
}

#' @noRd
.morgana_rule <- function(title = NULL, width = 72) {
  if (is.null(title)) {
    cat(strrep("=", width), "\n", sep = "")
  } else {
    pad <- max(width - nchar(title) - 4, 0)
    cat("-- ", title, " ", strrep("-", pad), "\n", sep = "")
  }
}

#' @export
print.morgana_solution <- function(x, digits = 3, n_doublets = NULL, ...) {
  s <- x$settings
  title <- if (inherits(x, "morgana_efa")) {
    "MORGANA Exploratory Factor Analysis"
  } else {
    "MORGANA Confirmatory (Mild-Restricted) Factor Analysis"
  }

  .morgana_rule()
  cat(title, "\n")
  cat("Method: Ferrando, Hernandez-Dorado & Lorenzo-Seva (unofficial R implementation)\n")
  .morgana_rule()
  cat(sprintf(
    "N = %d observations | %d items | %d factor(s) | %s %s\n",
    s$n, s$m, s$r, s$cor,
    if (s$cor == "covariance") "matrix" else "correlations"
  ))
  cat("Rotation:", s$rotate, "\n\n")

  .morgana_rule("Doublets (correlated residuals)")
  print(x$doublets, n = n_doublets)
  cat("\n")

  .morgana_rule("Goodness of fit (LOSEFER-corrected)")
  fit <- stats::setNames(x$fit$value, x$fit$measure)
  cat(sprintf(
    "Chi-square = %.3f on %d df (null model: %.2f on %d df)\n",
    fit["chi2_losefer"], as.integer(fit["df"]),
    fit["chi2_null"], as.integer(fit["df_null"])
  ))
  cat(sprintf(
    "RMSEA = %.3f | CFI = %.3f | NNFI/TLI = %.3f\n\n",
    fit["rmsea"], fit["cfi"], fit["nnfi_tli"]
  ))

  loading_title <- if (s$rotate == "none") {
    "Loadings (unrotated)"
  } else {
    sprintf("Loadings (%s pattern)", s$rotate)
  }
  .morgana_rule(loading_title)
  print(x$loadings, n = Inf)
  cat("\n")

  if (!is.null(x$Phi) && s$rotate == "oblimin") {
    .morgana_rule("Factor correlations (Phi)")
    print(round(x$Phi, digits))
    cat("\n")
  }

  cat("Residual correlations: $residual_correlation | unrotated loadings: $loadings_unrotated\n")
  if (!is.null(x$pa)) {
    cat("Parallel Analysis thresholds: $pa\n")
  }
  invisible(x)
}

#' @export
print.morgana_pa <- function(x, digits = 3, ...) {
  .morgana_rule("Parallel Analysis for EREC doublet detection")
  cat(sprintf(
    "%d replication(s) | cutoff at the %gth percentile per rank position\n\n",
    x$pa_reps, x$percentile
  ))

  tbl <- tibble::tibble(
    rank = seq_len(x$nd),
    item_a = x$item_names[x$doublets[, 1]],
    item_b = x$item_names[x$doublets[, 2]],
    erec = unname(x$doublets[, 3]),
    pa_threshold = x$pa_thresholds,
    significant = x$significant
  )
  print(tbl, n = Inf)

  n_sig <- sum(x$significant)
  cat(sprintf(
    "\n%d of %d doublet(s) exceed the chance threshold.\n", n_sig, x$nd
  ))
  invisible(x)
}
