test_that("fit_matern_copula isotropic recovers rho", {
  set.seed(42)
  true_rho <- 0.6
  Y <- rmatern_copula(n = 200, dim = c(8, 8), rho = true_rho, nu = 0)
  fit <- fit_matern_copula(Y, dim = c(8, 8), nu = 0, isotropic = TRUE, start = 0.3)

  expect_true(fit$convergence > 0)
  expect_true(abs(fit$rho - true_rho) < 0.15)
})

test_that("fit_matern_copula anisotropic recovers rho1 and rho2", {
  set.seed(42)
  rho1 <- 0.7
  rho2 <- 0.3
  Y <- rmatern_copula(n = 300, dim = c(8, 8), rho = c(rho1, rho2), nu = 0)
  fit <- fit_matern_copula(Y, dim = c(8, 8), nu = 0, start = c(0.5, 0.5))

  expect_true(fit$convergence > 0)
  expect_true(abs(fit$rho1 - rho1) < 0.2)
  expect_true(abs(fit$rho2 - rho2) < 0.2)
})

test_that("fit_matern_copula returns expected list elements (isotropic)", {
  set.seed(42)
  Y <- rmatern_copula(n = 50, dim = c(5, 5), rho = 0.5, nu = 0)
  fit <- fit_matern_copula(Y, dim = c(5, 5), nu = 0, isotropic = TRUE)

  expect_true("rho" %in% names(fit))
  expect_true("loglik" %in% names(fit))
  expect_true("convergence" %in% names(fit))
  expect_true("message" %in% names(fit))
  expect_true(is.numeric(fit$rho))
  expect_true(is.numeric(fit$loglik))
})

test_that("fit_matern_copula returns expected list elements (anisotropic)", {
  set.seed(42)
  Y <- rmatern_copula(n = 50, dim = c(5, 5), rho = 0.5, nu = 0)
  fit <- fit_matern_copula(Y, dim = c(5, 5), nu = 0, isotropic = FALSE)

  expect_true("rho1" %in% names(fit))
  expect_true("rho2" %in% names(fit))
  expect_true("rho" %in% names(fit))
  expect_length(fit$rho, 2)
})

test_that("fit_matern_copula respects bounds", {
  set.seed(42)
  Y <- rmatern_copula(n = 50, dim = c(5, 5), rho = 0.5, nu = 0)

  fit_iso <- fit_matern_copula(Y, dim = c(5, 5), nu = 0, isotropic = TRUE)
  expect_true(fit_iso$rho >= 0.001 && fit_iso$rho <= 0.999)

  fit_aniso <- fit_matern_copula(Y, dim = c(5, 5), nu = 0, isotropic = FALSE)
  expect_true(fit_aniso$rho1 >= 0.001 && fit_aniso$rho1 <= 0.999)
  expect_true(fit_aniso$rho2 >= 0.001 && fit_aniso$rho2 <= 0.999)
})

test_that("fit_matern_copula works with nu = 1", {
  set.seed(42)
  Y <- rmatern_copula(n = 100, dim = c(6, 6), rho = 0.5, nu = 1)
  fit <- fit_matern_copula(Y, dim = c(6, 6), nu = 1, isotropic = TRUE, start = 0.3)

  expect_true(fit$convergence > 0)
  expect_true(is.finite(fit$rho))
  expect_true(is.finite(fit$loglik))
})

test_that("fit_matern_copula works with nu = 2", {
  set.seed(42)
  Y <- rmatern_copula(n = 100, dim = c(6, 6), rho = 0.5, nu = 2)
  fit <- fit_matern_copula(Y, dim = c(6, 6), nu = 2, isotropic = TRUE, start = 0.3)

  expect_true(fit$convergence > 0)
  expect_true(is.finite(fit$rho))
  expect_true(is.finite(fit$loglik))
})

test_that("fit_matern_copula handles high rho with nu = 2 without hanging", {
  set.seed(42)
  # This is the regime that previously caused hangs
  Y <- rmatern_copula(n = 50, dim = c(6, 6), rho = 0.9, nu = 2)
  fit <- fit_matern_copula(Y, dim = c(6, 6), nu = 2, isotropic = TRUE, start = 0.5)

  # Should complete (not hang) and return finite values
  expect_true(is.finite(fit$rho))
  expect_true(is.finite(fit$loglik))
})

test_that("fit_matern_copula uses scalar dim correctly", {
  set.seed(42)
  Y <- rmatern_copula(n = 50, dim = 5, rho = 0.5, nu = 0)
  fit <- fit_matern_copula(Y, dim = 5, nu = 0, isotropic = TRUE)
  expect_true(fit$convergence > 0)
})

test_that("NaN guard prevents infinite loop", {
  # Directly call the C++ function with extreme values that could produce NaN
  # This tests the NaN guard we added
  set.seed(42)
  # Create data with extreme values
  X <- matrix(rnorm(36 * 10, sd = 10), nrow = 36, ncol = 10)
  # Should not hang; may return large penalty value
  ld <- dmatern_copula(X, dim = c(6, 6), rho = 0.98, nu = 2)
  expect_length(ld, 10)
  # Values should be finite (even if very negative)
  expect_true(all(is.finite(ld)))
})
