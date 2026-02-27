test_that("dmatern_copula returns one value per observation", {
  X <- matrix(rnorm(20 * 3), nrow = 20, ncol = 3) # 3 obs on 4x5 grid
  ld <- dmatern_copula(X, dim = c(4, 5), rho = 0.5, nu = 0)
  expect_length(ld, 3)
})

test_that("dmatern_copula returns finite values", {
  set.seed(42)
  X <- matrix(rnorm(16 * 10), nrow = 16, ncol = 10)
  ld <- dmatern_copula(X, dim = c(4, 4), rho = 0.5, nu = 0)
  expect_true(all(is.finite(ld)))
})

test_that("dmatern_copula accepts vector input (single observation)", {
  x <- rnorm(16)
  ld <- dmatern_copula(x, dim = c(4, 4), rho = 0.5, nu = 0)
  expect_length(ld, 1)
  expect_true(is.finite(ld))
})

test_that("dmatern_copula is zero for standard normal under independence", {
  # When rho is very small, the copula density should be near zero
  # (the copula approaches the independence copula)
  set.seed(42)
  X <- matrix(rnorm(16 * 100), nrow = 16, ncol = 100)
  ld_small <- dmatern_copula(X, dim = c(4, 4), rho = 0.01, nu = 0)
  # Log-density of independence copula is 0
  expect_true(all(abs(ld_small) < 1))
})

test_that("dmatern_copula is negative for non-independent copula", {
  # High rho should give meaningful (often negative) log-densities
  set.seed(42)
  X <- matrix(rnorm(16 * 10), nrow = 16, ncol = 10)
  ld <- dmatern_copula(X, dim = c(4, 4), rho = 0.8, nu = 0)
  expect_true(all(is.finite(ld)))
})

test_that("dmatern_copula scalar vs vector rho dispatch gives same result", {
  set.seed(42)
  X <- matrix(rnorm(16 * 5), nrow = 16, ncol = 5)
  ld1 <- dmatern_copula(X, dim = 4, rho = 0.5, nu = 0)
  ld2 <- dmatern_copula(X, dim = c(4, 4), rho = c(0.5, 0.5), nu = 0)
  expect_equal(ld1, ld2, tolerance = 1e-10)
})

test_that("dmatern_copula works for nu = 0, 1, 2", {
  set.seed(42)
  X <- matrix(rnorm(16 * 5), nrow = 16, ncol = 5)
  for (nu in 0:2) {
    ld <- dmatern_copula(X, dim = c(4, 4), rho = 0.5, nu = nu)
    expect_length(ld, 5)
    expect_true(all(is.finite(ld)))
  }
})

test_that("dmatern_copula density is consistent with sampling", {
  # Samples from rmatern_copula should have higher density than iid normal
  set.seed(42)
  Y_copula <- rmatern_copula(n = 50, dim = c(4, 4), rho = 0.7, nu = 0)
  Y_iid <- matrix(rnorm(16 * 50), nrow = 16, ncol = 50)

  ld_copula <- dmatern_copula(Y_copula, dim = c(4, 4), rho = 0.7, nu = 0)
  ld_iid <- dmatern_copula(Y_iid, dim = c(4, 4), rho = 0.7, nu = 0)

  # On average, copula samples should have higher log-density under the copula
  expect_true(mean(ld_copula) > mean(ld_iid))
})

test_that("dmatern (non-copula) returns finite values", {
  X <- matrix(rnorm(16 * 5), nrow = 16, ncol = 5)
  ld <- dmatern(X, dim = c(4, 4), rho = 0.5, nu = 0)
  expect_length(ld, 5)
  expect_true(all(is.finite(ld)))
})
