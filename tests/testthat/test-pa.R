test_that("ParallelAnalysisEREC returns thresholds per rank position", {
  X <- morgana_example_data()

  set.seed(42)
  pa <- ParallelAnalysisEREC(X, r = 2, cor = 2, pa_reps = 20)

  expect_s3_class(pa, "morgana_pa")
  expect_equal(pa$nd, 4)
  expect_length(pa$pa_thresholds, 4)
  expect_length(pa$significant, 4)
  expect_true(all(pa$pa_thresholds > 0))

  # The two true doublets have EREC around .5; 20 null replications are
  # more than enough to flag them.
  expect_true(all(pa$significant[1:2]))

  expect_output(print(pa), "Parallel Analysis")
  expect_output(print(pa), "exceed the chance threshold")
})

test_that("reusing the observed EREC leaves the PA bit-for-bit identical", {
  X <- morgana_example_data()
  R <- stats::cor(X)

  set.seed(99)
  fresh <- morganafa:::.parallel_analysis_erec(X, 2, 2L, R, 10, 95)

  obs <- LargestDoubletsEREC(X, 2, 2, R = R)
  set.seed(99)
  reused <- morganafa:::.parallel_analysis_erec(
    X, 2, 2L, R, 10, 95,
    obs_result = list(doublets = obs$doublets, nd = obs$nd)
  )

  expect_identical(reused, fresh)
})

test_that("parallel execution gives results identical to sequential", {
  skip_on_cran()
  # Workers load the installed morganafa; skip if it is not this version.
  installed <- tryCatch(
    read.dcf(file.path(
      find.package("morganafa", lib.loc = .libPaths()), "DESCRIPTION"
    ))[, "Version"],
    error = function(e) NA_character_
  )
  skip_if_not(
    identical(unname(installed), as.character(utils::packageVersion("morganafa"))),
    "installed morganafa differs from the version under test"
  )
  X <- morgana_example_data()

  set.seed(11)
  pa_seq <- ParallelAnalysisEREC(X, r = 2, pa_reps = 12, cores = 1)
  after_seq <- stats::runif(1)
  set.seed(11)
  pa_par <- ParallelAnalysisEREC(X, r = 2, pa_reps = 12, cores = 2)
  after_par <- stats::runif(1)
  expect_identical(pa_par, pa_seq)
  expect_identical(after_par, after_seq)

  set.seed(12)
  efa_seq <- solutionMEFA(X, r = 2, montecarlo = 60, pa_reps = 6, cores = 1)
  set.seed(12)
  efa_par <- solutionMEFA(X, r = 2, montecarlo = 60, pa_reps = 6, cores = 2)
  expect_identical(efa_par, efa_seq)
})

test_that("cores is validated", {
  X <- morgana_example_data()
  expect_error(ParallelAnalysisEREC(X, r = 2, pa_reps = 2, cores = 0), "cores")
  expect_error(ParallelAnalysisEREC(X, r = 2, pa_reps = 2, cores = 1.5), "cores")
})
