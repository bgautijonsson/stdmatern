test_that("make_AR_prec_matrix returns correct structure", {
  Q <- make_AR_prec_matrix(5L, 0.5)
  expect_true(inherits(Q, "dgCMatrix") || inherits(Q, "dsCMatrix"))
  expect_equal(nrow(Q), 5)
  expect_equal(ncol(Q), 5)
})

test_that("AR(1) precision matrix is symmetric", {
  Q <- make_AR_prec_matrix(10L, 0.7)
  Q_dense <- as.matrix(Q)
  expect_equal(Q_dense, t(Q_dense), tolerance = 1e-12)
})

test_that("AR(1) precision matrix is positive definite", {
  Q <- make_AR_prec_matrix(8L, 0.8)
  eigenvalues <- eigen(as.matrix(Q), symmetric = TRUE, only.values = TRUE)$values
  expect_true(all(eigenvalues > 0))
})

test_that("AR(1) precision matrix has unit marginal variance", {
  # The inverse of the precision matrix (covariance matrix) should have
  # diagonal entries equal to 1
  Q <- make_AR_prec_matrix(6L, 0.6)
  Sigma <- solve(as.matrix(Q))
  expect_equal(diag(Sigma), rep(1, 6), tolerance = 1e-10)
})

test_that("AR(1) precision matrix inverse gives correct AR(1) correlations", {
  rho <- 0.7
  Q <- make_AR_prec_matrix(5L, rho)
  Sigma <- solve(as.matrix(Q))
  # Correlation between adjacent elements should be rho
  expect_equal(Sigma[1, 2], rho, tolerance = 1e-10)
  # Correlation at lag 2 should be rho^2
  expect_equal(Sigma[1, 3], rho^2, tolerance = 1e-10)
})

test_that("AR(1) precision matrix handles rho = 0 (independence)", {
  Q <- make_AR_prec_matrix(5L, 0.0)
  expect_equal(as.matrix(Q), diag(5), tolerance = 1e-12)
})

test_that("AR(1) precision matrix handles small rho", {
  Q <- make_AR_prec_matrix(5L, 0.01)
  expect_true(all(eigen(as.matrix(Q), only.values = TRUE)$values > 0))
})

test_that("AR(1) precision matrix handles rho near upper bound", {
  Q <- make_AR_prec_matrix(5L, 0.99)
  eigenvalues <- eigen(as.matrix(Q), symmetric = TRUE, only.values = TRUE)$values
  expect_true(all(eigenvalues > 0))
})
