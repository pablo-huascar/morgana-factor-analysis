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
