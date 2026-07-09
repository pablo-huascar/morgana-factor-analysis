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
