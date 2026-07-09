test_that("LargestDoubletsEREC recovers the known doublets in the example data", {
  X <- morgana_example_data()
  res <- LargestDoubletsEREC(X, r = 2, cor = 2)

  expect_equal(res$nd, 4)
  expect_equal(dim(res$doublets), c(4, 3))
  expect_equal(colnames(res$doublets), c("P1", "P2", "EREC"))

  # The example data contain doublets on (1, 2) and (6, 7); they must rank
  # first, and the EREC values are a deterministic regression check.
  expect_equal(unname(res$doublets[1, 1:2]), c(1, 2))
  expect_equal(unname(res$doublets[2, 1:2]), c(6, 7))
  expect_equal(unname(res$doublets[1, 3]), 0.5008042, tolerance = 1e-5)
  expect_equal(unname(res$doublets[2, 3]), 0.4929888, tolerance = 1e-5)

  # The full EREC matrix is symmetric
  expect_true(isSymmetric(res$EREC, tol = 1e-10))
})

test_that("DoubletsEREC is consistent with LargestDoubletsEREC", {
  X <- morgana_example_data()
  largest <- LargestDoubletsEREC(X, r = 2, cor = 2)
  user <- DoubletsEREC(X, r = 2, cor = 2,
                       PairsDoublets = largest$doublets[, 1:2])
  expect_equal(user$doublets[, "EREC"], largest$doublets[, "EREC"],
    tolerance = 1e-10
  )
})

test_that("DoubletsEREC validates its pairs", {
  X <- morgana_example_data()
  expect_error(DoubletsEREC(X, 2, 2, PairsDoublets = cbind(1, 1)), "distinct")
  expect_error(DoubletsEREC(X, 2, 2, PairsDoublets = cbind(0, 2)), "distinct|between")
  expect_warning(
    DoubletsEREC(X, 2, 2, PairsDoublets = cbind(1:6, c(2:6, 8))),
    "can be freed"
  )
})

test_that("LargestDoubletsEREC accepts a precomputed R and string cor codes", {
  X <- morgana_example_data()
  R <- cor(X)
  res1 <- LargestDoubletsEREC(X, r = 2, cor = 2)
  res2 <- LargestDoubletsEREC(X, r = 2, cor = "pearson", R = R)
  expect_equal(res1$doublets, res2$doublets)
})

test_that("LargestDoubletsEREC errors when no doublet slot is available", {
  X <- morgana_example_data()[, 1:4]
  expect_error(LargestDoubletsEREC(X, r = 2, cor = 2), "Ledermann")
})
