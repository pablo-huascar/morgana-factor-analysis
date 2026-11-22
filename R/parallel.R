# Optional parallel execution of the Monte Carlo replications.
#
# Design: every random draw (the row indices of each replication) is made
# in the main R process, in exactly the same order as in the sequential
# loop. Only the deterministic model fits are sent to the workers, and the
# results are collected in order. The output is therefore identical to the
# sequential run with the same seed, whatever the number of workers.
#
# Workers are fresh R processes (PSOCK clusters), never forks, so this also
# works in front-ends that forbid forking (e.g. Positron) and on Windows.
# They load the *installed* morganafa, which must be the same version.

# Validate the `cores` argument and return it as an integer.
#' @noRd
.morgana_check_cores <- function(cores) {
  if (!is.numeric(cores) || length(cores) != 1 || is.na(cores) ||
    cores < 1 || cores != round(cores)) {
    stop("`cores` must be a single whole number >= 1.", call. = FALSE)
  }
  as.integer(cores)
}

# Start a PSOCK cluster with `cores` workers, or return NULL for
# sequential execution. The caller is responsible for stopping it.
#' @noRd
.morgana_cluster <- function(cores) {
  cores <- .morgana_check_cores(cores)
  if (cores < 2L) {
    return(NULL)
  }
  cl <- parallel::makePSOCKcluster(cores)
  parallel::clusterCall(cl, loadNamespace, "morganafa")
  cl
}

# Pool of simulated observations held by each worker, set once per
# simulation so that only row indices travel at each replication.
.pool_env <- new.env(parent = emptyenv())

#' @noRd
.set_pool <- function(XP) {
  assign("XP", XP, envir = .pool_env)
  invisible(NULL)
}

# Apply `worker(ix, args)` to every index vector in `idx` on cluster `cl`,
# in chunks so that `on_chunk(n_done)` can update a progress bar. Results
# are returned in the order of `idx`.
#' @noRd
.pool_map <- function(cl, idx, worker, args, on_chunk = function(n) NULL) {
  chunk_size <- 4L * length(cl)
  chunks <- split(seq_along(idx), ceiling(seq_along(idx) / chunk_size))
  out <- vector("list", length(idx))
  for (k in chunks) {
    out[k] <- parallel::parLapply(cl, idx[k], worker, args = args)
    on_chunk(length(k))
  }
  out
}

# Worker for ParallelAnalysisEREC(): ranked EREC values of one replication.
# The null data are continuous, so Pearson correlations are used.
#' @noRd
.pa_worker <- function(ix, args) {
  LargestDoubletsEREC(.pool_env$XP[ix, ], args$r, 2)$doublets[, 3]
}

# Worker for LoseferMorganaH0(): chi-square statistic of one replication.
#' @noRd
.h0_worker <- function(ix, args) {
  .h0_c2(.pool_env$XP[ix, ], args$r, args$doblets, args$conver)
}

# Chi-square statistic of the MORGANA model fitted to one simulated sample
# (used by both the sequential and the parallel LOSEFER simulation).
#' @noRd
.h0_c2 <- function(Xi, r, doblets, conver) {
  N <- nrow(Xi)
  m <- ncol(Xi)
  Ri <- cor(Xi)

  result <- fufa1(Xi, r, doblets, conver, 2, R = Ri)
  Ai <- result$L
  RESCORi <- result$RESCOR

  AAi <- Ai %*% t(Ai)
  PSIi <- sqrt(abs(eye(m) - diag(diag(AAi))))
  RESi <- Ri - (AAi + PSIi %*% (RESCORi - eye(m)) %*% PSIi)

  sum(RESi[upper.tri(RESi)]^2) * (N - 1)
}
