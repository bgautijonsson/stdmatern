test_that("rmatern_copula returns correct dimensions", {
  Y <- rmatern_copula_eigen(n = 10L, dim_x = 5L, dim_y = 4L, rho1 = 0.5, rho2 = 0.5, nu = 0L)
  expect_equal(nrow(Y), 20) # 5 * 4
  expect_equal(ncol(Y), 10)
})

test_that("rmatern_copula wrapper handles scalar dim and rho", {
  Y <- rmatern_copula(n = 5, dim = 4, rho = 0.5, nu = 0)
  expect_equal(nrow(Y), 16) # 4 * 4
  expect_equal(ncol(Y), 5)
})

test_that("rmatern_copula wrapper handles vector dim and rho", {
  Y <- rmatern_copula(n = 5, dim = c(3, 4), rho = c(0.3, 0.7), nu = 0)
  expect_equal(nrow(Y), 12) # 3 * 4
  expect_equal(ncol(Y), 5)
})

test_that("rmatern_copula samples have approximately zero mean", {
  set.seed(123)
  Y <- rmatern_copula(n = 500, dim = c(5, 5), rho = 0.5, nu = 0)
  col_means <- rowMeans(Y)
  # With 500 samples, means should be close to 0
  expect_true(all(abs(col_means) < 0.3))
})

test_that("rmatern_copula samples have approximately unit marginal variance", {
  set.seed(123)
  Y <- rmatern_copula(n = 1000, dim = c(5, 5), rho = 0.5, nu = 0)
  marginal_vars <- apply(Y, 1, var)
  # Marginal variances should be approximately 1
  expect_true(all(abs(marginal_vars - 1) < 0.3))
})

test_that("rmatern_copula works with nu = 0, 1, 2", {
  for (nu in 0:2) {
    Y <- rmatern_copula(n = 5, dim = c(4, 4), rho = 0.5, nu = nu)
    expect_equal(nrow(Y), 16, label = paste0("nu=", nu))
    expect_equal(ncol(Y), 5, label = paste0("nu=", nu))
    expect_true(all(is.finite(Y)), label = paste0("nu=", nu, " all finite"))
  }
})

test_that("rmatern_copula works with high rho", {
  Y <- rmatern_copula(n = 5, dim = c(4, 4), rho = 0.95, nu = 0)
  expect_true(all(is.finite(Y)))
})

test_that("rmatern_copula works with low rho", {
  Y <- rmatern_copula(n = 5, dim = c(4, 4), rho = 0.01, nu = 0)
  expect_true(all(is.finite(Y)))
})

test_that("rmatern_copula with high rho and nu=2 produces finite values", {
  Y <- rmatern_copula(n = 5, dim = c(4, 4), rho = 0.9, nu = 2)
  expect_true(all(is.finite(Y)))
})

test_that("rmatern (non-copula) returns correct dimensions", {
  Y <- rmatern(n = 5, dim = c(3, 4), rho = 0.5, nu = 0)
  expect_equal(nrow(Y), 12)
  expect_equal(ncol(Y), 5)
  expect_true(all(is.finite(Y)))
})
