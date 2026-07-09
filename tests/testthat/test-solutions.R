test_that("solutionMEFA returns a tidy, labelled solution", {
  X <- morgana_example_data()
  labels <- paste("Item text", seq_len(ncol(X)))

  set.seed(123)
  res <- solutionMEFA(X, r = 2, montecarlo = 200, pa_reps = 0,
                      item_labels = labels)

  expect_s3_class(res, "morgana_efa")
  expect_s3_class(res, "morgana_solution")
  expect_s3_class(res$loadings, "tbl_df")
  expect_s3_class(res$doublets, "tbl_df")
  expect_s3_class(res$fit, "tbl_df")

  expect_named(
    res$loadings,
    c("item", "label", "F1", "F2", "h2", "u2")
  )
  expect_equal(res$loadings$item, colnames(X))
  expect_equal(res$loadings$label, labels)
  expect_equal(res$loadings$h2 + res$loadings$u2, rep(1, ncol(X)),
    tolerance = 1e-6
  )

  expect_named(
    res$doublets,
    c("rank", "item_a", "item_b", "label_a", "label_b", "erec", "residual_r")
  )
  expect_equal(res$doublets$item_a[1:2], c("V1", "V6"))
  expect_equal(res$doublets$item_b[1:2], c("V2", "V7"))

  expect_equal(dim(res$residual_correlation), c(10, 10))
  expect_equal(rownames(res$residual_correlation), colnames(X))

  # Oblique rotation returns Phi; fit indices are in sensible ranges
  expect_equal(dim(res$Phi), c(2, 2))
  fit <- setNames(res$fit$value, res$fit$measure)
  expect_gt(fit["cfi"], 0.9)
  expect_lt(fit["rmsea"], 0.08)

  expect_output(print(res), "MORGANA Exploratory Factor Analysis")
  expect_output(print(res), "RMSEA")
})

test_that("solutionMCFA accepts index and name pairs and returns a solution", {
  X <- morgana_example_data()

  set.seed(123)
  res_idx <- solutionMCFA(X, r = 2, doublets = rbind(c(1, 2), c(6, 7)),
                          montecarlo = 200)
  set.seed(123)
  res_nm <- solutionMCFA(X, r = 2,
                         doublets = rbind(c("V1", "V2"), c("V6", "V7")),
                         montecarlo = 200)

  expect_s3_class(res_idx, "morgana_cfa")
  expect_equal(res_idx$doublets$erec, res_nm$doublets$erec)
  expect_equal(res_idx$doublets$item_a, c("V1", "V6"))
  expect_output(print(res_idx), "Confirmatory")
})

test_that("solution input validation catches common mistakes", {
  X <- morgana_example_data()

  expect_error(solutionMEFA(data.frame(a = letters), r = 1), "numeric")
  Xna <- X
  Xna[1, 1] <- NA
  expect_error(solutionMEFA(Xna, r = 2), "missing values")

  expect_error(
    solutionMEFA(X, r = 2, item_labels = c("only one label"), pa_reps = 0),
    "one entry per item"
  )
  expect_error(
    solutionMEFA(X, r = 2, item_labels = c(bad_name = "text"), pa_reps = 0),
    "do not match"
  )

  expect_error(solutionMCFA(X, r = 2, doublets = cbind(1, 1)), "distinct")
  expect_error(
    solutionMCFA(X, r = 2, doublets = rbind(c(1, 2), c(2, 1))),
    "duplicated"
  )
  expect_error(
    solutionMCFA(X, r = 2, doublets = rbind(c("V1", "nope"))),
    "Unknown item"
  )
})

test_that("unrotated loadings reproduce the MORGANA extraction", {
  X <- morgana_example_data()
  set.seed(123)
  res <- solutionMEFA(X, r = 2, montecarlo = 200, pa_reps = 0)

  extracted <- MorganaEFA(X, r = 2, conver = 1e-4, maxiter = 500, cor = 2)
  A <- posit(extracted$A)
  A_pkg <- as.matrix(res$loadings_unrotated[, c("F1", "F2")])
  expect_equal(unname(A_pkg), unname(A), tolerance = 1e-10)
})
