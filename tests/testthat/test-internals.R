test_that("posit reflects columns with negative loading sums", {
  A <- cbind(c(-.7, -.8, .1), c(.5, .6, .4))
  out <- posit(A)
  expect_equal(out[, 1], -A[, 1])
  expect_equal(out[, 2], A[, 2])
})

test_that("cor argument accepts codes and names", {
  expect_equal(morganafa:::.morgana_cor_code(1), 1L)
  expect_equal(morganafa:::.morgana_cor_code("pearson"), 2L)
  expect_equal(morganafa:::.morgana_cor_code("Polychoric"), 3L)
  expect_error(morganafa:::.morgana_cor_code(4), "must be")
  expect_error(morganafa:::.morgana_cor_code("spearman"), "must be")
})

test_that("BERGE returns the unit-norm least-squares solution", {
  set.seed(1)
  f <- matrix(rnorm(30), nrow = 10, ncol = 3)
  f[4, ] <- 0 # mimic the zeroed row used by the MINRES update
  phi <- matrix(rnorm(10), ncol = 1)

  t_hat <- morganafa:::BERGE(f, phi)
  expect_equal(sum(t_hat^2), 1, tolerance = 1e-6)

  # No unit-norm vector does better (spot check over random candidates)
  loss <- function(w) sum((f %*% w - phi)^2)
  best_random <- min(vapply(1:500, function(i) {
    w <- rnorm(3)
    loss(w / sqrt(sum(w^2)))
  }, numeric(1)))
  expect_lte(loss(t_hat), best_random + 1e-8)
})

test_that("minres with Heywood protection keeps communalities admissible", {
  X <- morgana_example_data()
  R <- cor(X)
  fit <- morganafa:::my_minres_hc(R, 2, 1e-5, nrow(X))
  h2 <- rowSums(as.matrix(fit$Loadings)^2)
  expect_true(all(h2 <= 1 + 1e-8))
})

test_that("polychoric dispersion does not fork the session in Positron", {
  skip_if_not_installed("psych")
  set.seed(1)
  X <- matrix(sample(1:5, 200 * 4, replace = TRUE), ncol = 4)

  withr_env <- Sys.getenv("POSITRON", unset = NA)
  Sys.setenv(POSITRON = "1")
  on.exit(
    if (is.na(withr_env)) Sys.unsetenv("POSITRON") else Sys.setenv(POSITRON = withr_env),
    add = TRUE
  )
  old <- options(mc.cores = 4L)
  on.exit(options(old), add = TRUE)

  # Inside the call mc.cores must be forced to 1 (no fork)...
  seen <- NULL
  testthat::local_mocked_bindings(
    polychoric = function(x, ...) {
      seen <<- getOption("mc.cores")
      list(rho = stats::cor(x))
    },
    .package = "psych"
  )
  R <- morganafa:::.morgana_dispersion(X, cor = 3)

  expect_equal(seen, 1L)
  expect_equal(dim(R), c(4, 4))
  # ...and the user's option is restored afterwards.
  expect_equal(getOption("mc.cores"), 4L)
})
