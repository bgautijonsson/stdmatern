test_that("make_standardized_matern returns sparse matrix of correct size", {
  Q <- make_standardized_matern(dim = c(4, 5), rho = 0.5, nu = 0)
  expect_s4_class(Q, "dgCMatrix")
  expect_equal(nrow(Q), 20)
  expect_equal(ncol(Q), 20)
})

test_that("make_standardized_matern is symmetric", {
  Q <- make_standardized_matern(dim = c(5, 5), rho = 0.6, nu = 1)
  Q_dense <- as.matrix(Q)
  expect_equal(Q_dense, t(Q_dense), tolerance = 1e-10)
})

test_that("make_standardized_matern is positive definite", {
  Q <- make_standardized_matern(dim = c(4, 4), rho = 0.5, nu = 0)
  eigenvalues <- eigen(as.matrix(Q), symmetric = TRUE, only.values = TRUE)$values
  expect_true(all(eigenvalues > 0))
})

test_that("standardized precision inverse is a correlation matrix", {
  Q <- make_standardized_matern(dim = c(4, 4), rho = 0.5, nu = 0)
  Sigma <- solve(as.matrix(Q))
  # Diagonal of inverse should be 1 (unit marginal variance = correlation matrix)
  expect_equal(diag(Sigma), rep(1, 16), tolerance = 1e-8)
})

test_that("scalar dim creates square grid", {
  Q1 <- make_standardized_matern(dim = 5, rho = 0.5, nu = 0)
  Q2 <- make_standardized_matern(dim = c(5, 5), rho = 0.5, nu = 0)
  expect_equal(as.matrix(Q1), as.matrix(Q2), tolerance = 1e-12)
})

test_that("scalar rho gives isotropic precision matrix", {
  Q1 <- make_standardized_matern(dim = c(5, 5), rho = 0.6, nu = 0)
  Q2 <- make_standardized_matern(dim = c(5, 5), rho = c(0.6, 0.6), nu = 0)
  expect_equal(as.matrix(Q1), as.matrix(Q2), tolerance = 1e-12)
})

test_that("anisotropic rho gives different result from isotropic", {
  Q1 <- make_standardized_matern(dim = c(5, 5), rho = c(0.3, 0.8), nu = 0)
  Q2 <- make_standardized_matern(dim = c(5, 5), rho = 0.5, nu = 0)
  expect_false(isTRUE(all.equal(as.matrix(Q1), as.matrix(Q2))))
})

test_that("higher nu increases sparsity pattern bandwidth", {
  Q0 <- make_standardized_matern(dim = c(5, 5), rho = 0.5, nu = 0)
  Q1 <- make_standardized_matern(dim = c(5, 5), rho = 0.5, nu = 1)
  # Higher nu means more powers of the base matrix, so more nonzeros
  expect_true(Matrix::nnzero(Q1) >= Matrix::nnzero(Q0))
})

test_that("nu = 0, 1, 2 all produce valid precision matrices", {
  for (nu in 0:2) {
    Q <- make_standardized_matern(dim = c(4, 4), rho = 0.5, nu = nu)
    Sigma <- solve(as.matrix(Q))
    expect_equal(diag(Sigma), rep(1, 16), tolerance = 1e-6,
                 label = paste0("nu=", nu, " diagonal"))
    eigenvalues <- eigen(as.matrix(Q), symmetric = TRUE, only.values = TRUE)$values
    expect_true(all(eigenvalues > 0), label = paste0("nu=", nu, " positive definite"))
  }
})
