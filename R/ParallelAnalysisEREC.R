#' Parallel Analysis reference thresholds for the EREC index
#'
#' Generates a null EREC distribution by Monte Carlo simulation from a
#' factor model with no correlated residuals, and returns Parallel Analysis
#' cutoffs for each rank position so observed doublets can be compared
#' against chance (see the "Reference values and upper bounds" section of
#' Ferrando, Hernandez-Dorado & Lorenzo-Seva, 2022).
#'
#' A doublet is flagged as *significant* when its observed EREC exceeds the
#' PA cutoff at the same rank position.
#'
#' @section Computational cost:
#' Each replication refits the model once per item pair, i.e.
#' m(m-1)/2 MINRES fits, so the running time grows with the square of the
#' number of items and linearly with `pa_reps`. As a rough guide (Apple
#' silicon laptop, r = 3): one replication takes about 0.2 s with 20 items
#' and about 5 s with 58 items, so the default 500 replications take
#' minutes to most of an hour. Use e.g. `pa_reps = 100` while exploring
#' and keep the default for the final analysis. A progress bar is shown in
#' interactive sessions.
#'
#' @section Parallel execution:
#' With `cores > 1` the replications are distributed over that many fresh R
#' worker processes (a PSOCK cluster; no forking, so it works in Positron,
#' RStudio and on Windows). All random draws are still made in the main
#' process, in the same order as the sequential loop, and only the
#' deterministic model fits run on the workers. Results are therefore
#' **identical** to `cores = 1` with the same [set.seed()], whatever the
#' number of cores. Starting the workers takes a few seconds, so it only
#' pays off for long runs. A good choice is the number of physical cores
#' minus one (see [parallel::detectCores()]). The workers load the
#' installed version of morganafa.
#'
#' The dispersion matrix of the observed data is computed only once; pass
#' it through `R` when calling several functions on the same data (this
#' matters most for polychoric correlations).
#'
#' @inheritParams LargestDoubletsEREC
#' @param pa_reps Number of Monte Carlo replications (default 500).
#' @param percentile Percentile used as the cutoff (default 95).
#' @param cores Number of parallel worker processes. `1` (default) runs
#'   sequentially; see the "Parallel execution" section. The default can be
#'   set once per session with `options(morganafa.cores = 8)`.
#'
#' @return An object of class `morgana_pa`: a list with
#'   \item{doublets}{Matrix (`nd` x 3) of observed doublets (`P1`, `P2`,
#'     `EREC`).}
#'   \item{pa_thresholds}{PA cutoff at the requested percentile for each
#'     rank position.}
#'   \item{significant}{Logical: observed EREC > cutoff, per rank.}
#'   \item{nd, pa_reps, percentile}{Settings used.}
#'   \item{item_names}{Column names of `X`.}
#'
#' @inherit LargestDoubletsEREC references
#' @seealso [solutionMEFA()], [LargestDoubletsEREC()]
#' @export
#' @examples
#' \donttest{
#' database <- as.matrix(read.table(
#'   system.file("extdata", "database.dat", package = "morganafa")
#' ))
#' set.seed(1)
#' pa <- ParallelAnalysisEREC(database, r = 2, cor = 2, pa_reps = 50)
#' pa
#' }
ParallelAnalysisEREC <- function(X, r, cor = 2, pa_reps = 500,
                                 percentile = 95, R = NULL,
                                 cores = getOption("morganafa.cores", 1L)) {
  X <- .morgana_check_data(X, arg = "X")
  cor_code <- .morgana_cor_code(cor)
  R <- .morgana_dispersion(X, cor_code, R)
  cl <- .morgana_cluster(cores)
  if (!is.null(cl)) {
    on.exit(parallel::stopCluster(cl), add = TRUE)
  }
  .parallel_analysis_erec(X, r, cor_code, R, pa_reps, percentile, cl = cl)
}

# Workhorse of ParallelAnalysisEREC(). `obs_result` lets callers that have
# already run LargestDoubletsEREC() on the same X, r, cor and R (e.g.
# solutionMEFA() via MorganaEFA()) pass it in instead of recomputing it.
# That computation is deterministic, so reusing it leaves the results and
# the random-number stream unchanged.
#' @noRd
.parallel_analysis_erec <- function(X, r, cor_code, R, pa_reps, percentile,
                                    obs_result = NULL, cl = NULL) {
  N <- nrow(X)
  m <- ncol(X)
  conver <- 0.00001

  # Observed EREC and the number of doublet slots
  if (is.null(obs_result)) {
    obs_result <- LargestDoubletsEREC(X, r, cor_code, R = R)
  }
  nd <- obs_result$nd
  obs_doublets <- obs_result$doublets

  # Fit the baseline model (no doublets) to define the null population
  if (cor_code > 1) {
    A_null <- as.matrix(my_minres_hc(R, r, conver, N)$Loadings)
  } else {
    A_null <- as.matrix(my_minres(R, r, conver, N)$Loadings)
  }

  # Null correlation matrix: AA' + diag(uniqueness), unit diagonal
  AA <- A_null %*% t(A_null)
  Rr <- AA - diag(diag(AA)) + diag(m)

  # Ensure positive definiteness before the Cholesky factorisation
  min_ev <- min(eigen(Rr, symmetric = TRUE, only.values = TRUE)$values)
  if (min_ev < 0.001) {
    Rr <- Rr + diag(rep(0.001 - min_ev, m))
  }
  L_chol <- chol(Rr)

  # Large pool of null observations, sub-sampled at each replication
  NP <- 100000
  XP <- matrix(rnorm(NP * m), ncol = m) %*% L_chol

  # Monte Carlo: collect the nd ranked EREC values from each replication.
  # The null data are continuous, so Pearson correlations are used.
  pa_erec <- matrix(0, nrow = pa_reps, ncol = nd)
  cli::cli_progress_bar("Parallel Analysis (EREC)", total = pa_reps,
    auto_terminate = FALSE
  )
  if (is.null(cl)) {
    for (h in seq_len(pa_reps)) {
      Xi <- XP[sample.int(NP, N), ]
      pa_erec[h, ] <- LargestDoubletsEREC(Xi, r, 2)$doublets[, 3]
      cli::cli_progress_update()
    }
  } else {
    # All draws happen here, in the sequential order; workers only fit.
    idx <- lapply(seq_len(pa_reps), function(h) sample.int(NP, N))
    parallel::clusterCall(cl, .set_pool, XP)
    pb_env <- environment()
    res <- .pool_map(cl, idx, .pa_worker, list(r = r),
      on_chunk = function(n) cli::cli_progress_update(inc = n, .envir = pb_env)
    )
    pa_erec[] <- do.call(rbind, res)
  }
  cli::cli_progress_done()

  pa_thresholds <- apply(pa_erec, 2, quantile, probs = percentile / 100)

  structure(
    list(
      doublets = obs_doublets,
      pa_thresholds = unname(pa_thresholds),
      significant = unname(obs_doublets[, 3] > pa_thresholds),
      nd = nd,
      pa_reps = pa_reps,
      percentile = percentile,
      item_names = colnames(X)
    ),
    class = "morgana_pa"
  )
}
