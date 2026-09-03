# Input validation at the R boundary (2026-09-02). Before the guards: a wrong nrow returned a finite
# number from an out-of-bounds read, NA/Inf in Y "converged" on the flat 1e15 guard, rho >= 1 gave an
# all-NaN precision silently, nu < 0 a wrong one, and dim = c(0, 4) segfaulted the R session.
Y_ok <- function() rmatern_copula(6L, c(4L, 4L), 0.5, 0L, seed = 1L)

test_that("dim must be positive integers of length 1 or 2", {
  expect_error(make_standardized_matern(c(4, 4, 4), 0.5, 0L), "dim")
  expect_error(make_standardized_matern(c(4.5, 4), 0.5, 0L), "dim")
  expect_error(rmatern_copula_eigen(2L, 0L, 4L, 0.5, 0.5, 0L), "dim")
  expect_error(dmatern_copula(matrix(0, 16, 1), c(-4, 4), 0.5, 0L), "dim")
})

test_that("dim = c(0, 4) is refused with a message instead of aborting the session", {
  skip_if_not_installed("callr")
  skip_if_not_installed("pkgload")
  root <- system.file(package = "stdmatern")
  dev <- dir.exists(file.path(root, "src"))
  expect_error(
    callr::r(function(root, dev, libs) {
      if (dev) pkgload::load_all(root, quiet = TRUE) else library(stdmatern, lib.loc = libs)
      make_standardized_matern(c(0L, 4L), 0.5, 0L)
    }, args = list(root = root, dev = dev, libs = .libPaths()),
    env = c(callr::rcmd_safe_env(), OMP_NUM_THREADS = "1")),
    "dim must be positive"
  )
})

test_that("rho must lie in (-1, 1)", {
  expect_error(make_standardized_matern(c(4, 4), 1, 0L), "rho")
  expect_error(make_standardized_matern(c(4, 4), 1.5, 0L), "rho")
  expect_error(rmatern_copula_eigen(2L, 4L, 4L, 1, 0.5, 0L), "rho")
  expect_error(dmatern_copula(matrix(0, 16, 1), c(4, 4), c(0.5, -1), 0L), "rho")
  expect_silent(make_standardized_matern(c(4, 4), -0.5, 0L))
})

test_that("nu must be a non-negative integer", {
  expect_error(make_standardized_matern(c(4, 4), 0.5, -1L), "nu")
  expect_error(make_standardized_matern(c(4, 4), 0.5, 1.7), "nu")
  expect_error(rmatern_copula_eigen(2L, 4L, 4L, 0.5, 0.5, -1L), "nu")
})

test_that("the number of rows must equal the number of grid cells", {
  expect_error(dmatern_copula(matrix(rnorm(10), 10, 1), c(4, 4), 0.5, 0L), "rows")
  expect_error(dmatern_copula_eigen(matrix(rnorm(20), 20, 1), 4L, 4L, 0.5, 0.5, 0L), "rows")
  expect_error(fit_matern_copula(matrix(rnorm(50), 10, 5), c(4, 4), 0L, isotropic = TRUE), "rows")
  expect_error(fit_matern_copula(Y_ok(), c(5, 5), 0L, isotropic = TRUE), "rows")
})

test_that("non-finite observations are refused", {
  Y <- Y_ok()
  Y[3, 2] <- NA
  expect_error(fit_matern_copula(Y, c(4, 4), 0L, isotropic = TRUE), "finite")
  Y[3, 2] <- Inf
  expect_error(dmatern_copula(Y, c(4, 4), 0.5, 0L), "finite")
})

test_that("seed must be a whole number and n non-negative", {
  expect_error(rmatern_copula_eigen(2L, 4L, 4L, 0.5, 0.5, 0L, seed = 3.7), "seed")
  expect_error(rmatern_copula(2L, 4L, 0.5, 0L, seed = "a"), "seed")
  expect_error(rmatern_copula_eigen(-1L, 4L, 4L, 0.5, 0.5, 0L), "n")
  expect_identical(dim(rmatern_copula_eigen(0L, 4L, 4L, 0.5, 0.5, 0L)), c(16L, 0L))
})

test_that("a start outside the search bounds is refused", {
  expect_error(fit_matern_copula(Y_ok(), c(4, 4), 0L, isotropic = TRUE, start = 1.5), "start")
})
