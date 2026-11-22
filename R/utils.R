# Internal helpers shared across the MORGANA routines.
#
# Most of these mimic MATLAB built-ins used in the original code by
# Hernandez-Dorado, Lorenzo-Seva & Ferrando so that the translated
# algorithms stay close to their source.

# Identity matrix (MATLAB eye)
#' @noRd
eye <- function(a) {
  diag(a)
}

# Matrix of zeros (MATLAB zeros)
#' @noRd
zeros <- function(rows, cols) {
  matrix(0, nrow = rows, ncol = cols)
}

# Matrix of ones (MATLAB ones)
#' @noRd
ones <- function(nrow, ncol = nrow) {
  matrix(1, nrow = nrow, ncol = ncol)
}

# Column minima (MATLAB min over rows)
#' @noRd
colmin <- function(x) {
  apply(x, 2, min)
}

# Sum of squared elements (MATLAB ssq)
#' @noRd
ssq <- function(a) {
  sum(a^2)
}

# Replace NaN entries by 0 and force a unit diagonal
#' @noRd
Nanas <- function(matriz) {
  matriz[is.nan(matriz)] <- 0
  diag(matriz) <- 1
  matriz
}

# Drop the -1 placeholders used by the original code to mark unused slots
#' @noRd
ClearBuits <- function(B) {
  B <- as.matrix(B)
  B[B != -1]
}

# Moore-Penrose pseudoinverse via SVD
#' @noRd
pinv <- function(A, tol = sqrt(.Machine$double.eps)) {
  A <- as.matrix(A)
  s <- svd(A)
  d <- s$d
  if (length(d) == 0) {
    return(matrix(0, ncol(A), nrow(A)))
  }
  keep <- d > max(tol * d[1], 0)
  if (!any(keep)) {
    return(matrix(0, ncol(A), nrow(A)))
  }
  s$v[, keep, drop = FALSE] %*% (t(s$u[, keep, drop = FALSE]) / d[keep])
}

# Backward-compatible alias used by original MORGANA routines
#' @noRd
ginv <- function(A, tol = sqrt(.Machine$double.eps)) {
  pinv(A, tol = tol)
}

# Inverse with pseudoinverse fallback for singular matrices
#' @noRd
inverse <- function(mat) {
  inv_mat <- tryCatch(solve(mat), error = function(e) NULL)
  if (is.null(inv_mat)) {
    list(result = pinv(mat), method = "pseudo-inverse")
  } else {
    list(result = inv_mat, method = "inverse")
  }
}

# Guttman's partial correlation (anti-image) matrix
#' @noRd
PartCor <- function(R) {
  m <- nrow(R)
  ir <- as.matrix(inverse(R)$result)
  Dd <- diag(diag(ir)^(-1 / 2))
  AI <- Dd %*% ir %*% Dd
  Ra <- 2 * diag(m) - AI
  Nanas(Ra)
}

# Normalise the `cor` argument to the integer codes used internally:
# 1 = covariance, 2 = Pearson, 3 = polychoric.
#' @noRd
.morgana_cor_code <- function(cor) {
  if (is.character(cor)) {
    code <- switch(tolower(cor[1]),
      "covariance" = 1L, "cov" = 1L,
      "pearson" = 2L, "cor" = 2L,
      "polychoric" = 3L, "poly" = 3L,
      stop("`cor` must be \"covariance\", \"pearson\", \"polychoric\" (or 1, 2, 3).",
        call. = FALSE
      )
    )
    return(code)
  }
  cor <- as.integer(cor[1])
  if (!cor %in% 1:3) {
    stop("`cor` must be 1 (covariance), 2 (Pearson), or 3 (polychoric).",
      call. = FALSE
    )
  }
  cor
}

#' @noRd
.morgana_cor_label <- function(code) {
  c("covariance", "Pearson", "polychoric")[code]
}

# Dispersion matrix used everywhere. Computing it once per analysis and
# passing it down avoids repeated (and, for polychoric, very expensive)
# recomputation inside the lower-level routines.
#' @noRd
.morgana_dispersion <- function(X, cor, R = NULL) {
  if (!is.null(R)) {
    return(as.matrix(R))
  }
  code <- .morgana_cor_code(cor)
  switch(code,
    cov(X),
    cor(X),
    .morgana_polychoric(X)
  )
}

# psych::polychoric() parallelises via parallel::mcmapply(), which forks the
# R session. Some front-ends (e.g. Positron) forbid forking, so there we
# force a single core; mclapply() then falls back to a plain lapply().
#' @noRd
.morgana_polychoric <- function(X) {
  if (nzchar(Sys.getenv("POSITRON"))) {
    old <- options(mc.cores = 1L)
    on.exit(options(old), add = TRUE)
  }
  psych::polychoric(X)$rho
}

# Validate/standardise the input data matrix
#' @noRd
.morgana_check_data <- function(data, arg = "data") {
  X <- as.matrix(data)
  if (!is.numeric(X)) {
    stop(sprintf("`%s` must be a numeric matrix or data frame.", arg),
      call. = FALSE
    )
  }
  if (anyNA(X)) {
    stop(sprintf("`%s` contains missing values; remove or impute them first.", arg),
      call. = FALSE
    )
  }
  if (is.null(colnames(X))) {
    colnames(X) <- paste0("V", seq_len(ncol(X)))
  }
  X
}
