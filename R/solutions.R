#' Complete MORGANA exploratory factor analysis
#'
#' Runs the full exploratory MORGANA pipeline: EREC-based doublet detection
#' ([LargestDoubletsEREC()]), MORGANA extraction with the detected doublets
#' freed ([Morgana()]), LOSEFER-corrected fit indices ([LoseferMorgana()]),
#' optional factor rotation, and optional Parallel Analysis reference
#' thresholds for the EREC values ([ParallelAnalysisEREC()]).
#'
#' @section Disclaimer:
#' This package is an independent, unofficial implementation of methods
#' proposed by Ferrando, Hernandez-Dorado and Lorenzo-Seva. The original
#' authors are not involved in this package, and its results have not yet
#' been formally compared with the official implementation in the FACTOR
#' program. See `help("morganafa-package")`.
#'
#' @param data Numeric data matrix or data frame (observations x items).
#'   Column names are used as item names in all outputs.
#' @param r Number of factors to extract.
#' @param cor Dispersion matrix to analyse: `"pearson"` (default),
#'   `"polychoric"`, or `"covariance"`. The numeric codes used by the
#'   original code (1 = covariance, 2 = Pearson, 3 = polychoric) are also
#'   accepted.
#' @param rotate Factor rotation applied for interpretation: `"oblimin"`
#'   (default when more than one factor is extracted), `"varimax"`, or
#'   `"none"`. Rotation does not
#'   affect the fit indices, which are computed from the unrotated solution.
#' @param item_labels Optional character vector with the text of each item
#'   (e.g. `"I have felt nervous lately"`). Either unnamed and in column
#'   order, or named with the column names of `data`. Labels are carried
#'   into the loading and doublet tables.
#' @param montecarlo Number of Monte Carlo replications for the LOSEFER H0
#'   distribution (default 5000).
#' @param pa_reps Number of Parallel Analysis replications for the EREC
#'   thresholds (default 500). Use `0` to skip the Parallel Analysis.
#' @param pa_percentile Percentile for the PA cutoffs (default 95).
#' @param conver Convergence criterion for the MORGANA iterations.
#' @param maxiter Maximum number of MORGANA iterations.
#'
#' @return An object of class `morgana_efa`. Its main components are
#'   [tibble][tibble::tibble]s, ready for further processing:
#'   \item{loadings}{Tibble with one row per item: `item`, optional `label`,
#'     the rotated loadings (`F1`, ..., `Fr`), communality `h2` and
#'     uniqueness `u2`.}
#'   \item{Phi}{Factor correlation matrix (oblique rotations only).}
#'   \item{loadings_unrotated}{Tibble with the unrotated MORGANA loadings.}
#'   \item{doublets}{Tibble with one row per freed doublet: the items
#'     involved, their optional labels, the `erec` estimate, the fitted
#'     residual correlation `residual_r` and, when the Parallel Analysis is
#'     run, `pa_threshold` and `significant`.}
#'   \item{fit}{Tibble of LOSEFER-corrected fit measures (chi-square, df,
#'     RMSEA, CFI, NNFI/TLI).}
#'   \item{residual_correlation, residual_covariance}{Residual matrices with
#'     item names.}
#'   \item{pa}{The [ParallelAnalysisEREC()] result (or `NULL`).}
#'   \item{settings}{List with the options used.}
#'
#' @references
#' Ferrando, P. J., Hernandez-Dorado, A., & Lorenzo-Seva, U. (2022).
#' Detecting correlated residuals in exploratory factor analysis: New
#' proposals and a comparison of procedures. *Structural Equation Modeling*,
#' 29(4), 630-638. \doi{10.1080/10705511.2021.2004543}
#'
#' Lorenzo-Seva, U. (2025). Mild-restricted confirmatory factor analysis: a
#' new approach to assessing population factor models. *Communications in
#' Statistics - Simulation and Computation*.
#' \doi{10.1080/03610918.2025.2453473}
#'
#' @seealso [solutionMCFA()] for the confirmatory workflow.
#' @export
#' @examples
#' \donttest{
#' database <- as.matrix(read.table(
#'   system.file("extdata", "database.dat", package = "morganafa")
#' ))
#' set.seed(123)
#' res <- solutionMEFA(database, r = 2, montecarlo = 500, pa_reps = 0)
#' res
#' res$loadings
#' res$doublets
#' }
solutionMEFA <- function(data,
                         r,
                         cor = "pearson",
                         rotate = if (r > 1) "oblimin" else "none",
                         item_labels = NULL,
                         montecarlo = 5000,
                         pa_reps = 500,
                         pa_percentile = 95,
                         conver = 0.0001,
                         maxiter = 500) {
  X <- .morgana_check_data(data)
  cor_code <- .morgana_cor_code(cor)
  item_names <- colnames(X)
  labels <- .morgana_item_labels(item_labels, item_names)

  N <- nrow(X)
  R <- .morgana_dispersion(X, cor_code)

  extracted <- MorganaEFA(X, r = r, conver = conver, maxiter = maxiter,
                          cor = cor_code, R = R)
  A <- posit(extracted$A)
  doublets <- extracted$doublets

  h0 <- LoseferMorganaH0(A, extracted$RESCOR,
                         doublets[, 1:2, drop = FALSE], N, montecarlo)
  fit <- LoseferMorgana(
    R, A, N, extracted$RESCOR, nrow(doublets),
    h0$v, h0$vv, h0$a_n, h0$b1_n, h0$b2_n,
    h0$a_f, h0$b_f, h0$c_f, h0$d_f,
    display = 0
  )

  pa <- NULL
  if (pa_reps > 0) {
    pa <- ParallelAnalysisEREC(X, r, cor = cor_code, pa_reps = pa_reps,
                               percentile = pa_percentile, R = R)
  }

  .morgana_solution(
    class = "morgana_efa",
    A = A, extracted = extracted, doublets = doublets,
    fit = fit, h0 = h0, pa = pa,
    rotate = rotate, item_names = item_names, labels = labels,
    settings = list(
      n = N, m = ncol(X), r = r,
      cor = .morgana_cor_label(cor_code), rotate = rotate,
      montecarlo = montecarlo, pa_reps = pa_reps,
      pa_percentile = pa_percentile,
      conver = conver, maxiter = maxiter
    )
  )
}

#' Complete MORGANA confirmatory (mild-restricted) factor analysis
#'
#' Runs the confirmatory MORGANA pipeline with user-specified doublets:
#' EREC estimation for the given pairs ([DoubletsEREC()]), MORGANA
#' extraction ([Morgana()]), LOSEFER-corrected fit indices
#' ([LoseferMorgana()]), and optional rotation. This follows the
#' mild-restricted CFA approach of Lorenzo-Seva (2025).
#'
#' @inheritSection solutionMEFA Disclaimer
#'
#' @inheritParams solutionMEFA
#' @param doublets Pairs of items whose residual correlations should be
#'   freed: a two-column matrix/data frame of item indices, or of item
#'   names matching the column names of `data`.
#'
#' @return An object of class `morgana_cfa`, with the same components as
#'   the result of [solutionMEFA()] (except `pa`).
#'
#' @inherit solutionMEFA references
#' @seealso [solutionMEFA()] for the exploratory workflow.
#' @export
#' @examples
#' \donttest{
#' database <- as.matrix(read.table(
#'   system.file("extdata", "database.dat", package = "morganafa")
#' ))
#' pairs <- rbind(c(1, 2), c(6, 7))
#' set.seed(123)
#' res <- solutionMCFA(database, r = 2, doublets = pairs, montecarlo = 500)
#' res
#' }
solutionMCFA <- function(data,
                         r,
                         doublets,
                         cor = "pearson",
                         rotate = if (r > 1) "oblimin" else "none",
                         item_labels = NULL,
                         montecarlo = 5000,
                         conver = 0.0001,
                         maxiter = 500) {
  X <- .morgana_check_data(data)
  cor_code <- .morgana_cor_code(cor)
  item_names <- colnames(X)
  labels <- .morgana_item_labels(item_labels, item_names)

  provided <- .morgana_check_pairs(doublets, item_names)

  N <- nrow(X)
  R <- .morgana_dispersion(X, cor_code)

  extracted <- MorganaCFA(X, r = r, conver = conver, maxiter = maxiter,
                          cor = cor_code, doblets = provided, R = R)
  A <- posit(extracted$A)
  doublets_out <- extracted$doublets

  h0 <- LoseferMorganaH0(A, extracted$RESCOR,
                         doublets_out[, 1:2, drop = FALSE], N, montecarlo)
  fit <- LoseferMorgana(
    R, A, N, extracted$RESCOR, nrow(doublets_out),
    h0$v, h0$vv, h0$a_n, h0$b1_n, h0$b2_n,
    h0$a_f, h0$b_f, h0$c_f, h0$d_f,
    display = 0
  )

  .morgana_solution(
    class = "morgana_cfa",
    A = A, extracted = extracted, doublets = doublets_out,
    fit = fit, h0 = h0, pa = NULL,
    rotate = rotate, item_names = item_names, labels = labels,
    settings = list(
      n = N, m = ncol(X), r = r,
      cor = .morgana_cor_label(cor_code), rotate = rotate,
      montecarlo = montecarlo,
      conver = conver, maxiter = maxiter
    )
  )
}

# Resolve `item_labels` into a character vector aligned with item_names
# (or NULL when not supplied).
#' @noRd
.morgana_item_labels <- function(item_labels, item_names) {
  if (is.null(item_labels)) {
    return(NULL)
  }
  item_labels <- unlist(item_labels)
  if (!is.null(names(item_labels)) && all(nzchar(names(item_labels)))) {
    unknown <- setdiff(names(item_labels), item_names)
    if (length(unknown) > 0) {
      stop("`item_labels` has names that do not match the data columns: ",
        paste(unknown, collapse = ", "),
        call. = FALSE
      )
    }
    out <- rep(NA_character_, length(item_names))
    names(out) <- item_names
    out[names(item_labels)] <- as.character(item_labels)
    return(unname(out))
  }
  if (length(item_labels) != length(item_names)) {
    stop("`item_labels` must have one entry per item (", length(item_names),
      ") or be named with the data column names.",
      call. = FALSE
    )
  }
  as.character(item_labels)
}

# Validate user-provided doublet pairs; accepts indices or item names.
#' @noRd
.morgana_check_pairs <- function(doublets, item_names) {
  m <- length(item_names)
  pairs <- as.matrix(doublets)
  if (is.null(dim(pairs)) || ncol(pairs) != 2) {
    stop("`doublets` must be a two-column matrix/data.frame of item pairs.",
      call. = FALSE
    )
  }
  if (is.character(pairs)) {
    idx <- match(pairs, item_names)
    if (anyNA(idx)) {
      stop("Unknown item name(s) in `doublets`: ",
        paste(unique(pairs[is.na(idx)]), collapse = ", "),
        call. = FALSE
      )
    }
    pairs <- matrix(idx, ncol = 2)
  }
  storage.mode(pairs) <- "double"
  if (any(pairs < 1) || any(pairs > m) || any(pairs[, 1] == pairs[, 2])) {
    stop("`doublets` must contain pairs of distinct item indices between 1 and ",
      m, ".",
      call. = FALSE
    )
  }
  key <- paste(pmin(pairs[, 1], pairs[, 2]), pmax(pairs[, 1], pairs[, 2]))
  if (anyDuplicated(key)) {
    stop("`doublets` contains duplicated pairs.", call. = FALSE)
  }
  pairs
}
