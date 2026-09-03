# fit_matern_copula's contract (2026-09-02): a status the caller can act on, an evaluation count, a
# deterministic evaluation budget instead of a wall-clock cap, and an opt-in O(N (n1 + n2) n) objective.
Yiso <- function(n = 60L, dims = c(12L, 9L), nu = 1L, rho = 0.6) rmatern_copula(n, dims, rho, nu, seed = 21L)

test_that("the fit reports evaluations, convergence and a stop reason", {
  f <- fit_matern_copula(Yiso(), c(12L, 9L), 1L, isotropic = TRUE)
  expect_true(is.integer(f$n_evals) && f$n_evals > 0L)
  expect_true(isTRUE(f$converged))
  expect_identical(f$stop_reason, "converged")
  expect_true(f$convergence %in% 1:4)
})

test_that("the anisotropic fit reports the same fields", {
  f <- fit_matern_copula(Yiso(), c(12L, 9L), 1L, isotropic = FALSE)
  expect_true(is.integer(f$n_evals) && f$n_evals > 0L)
  expect_true(isTRUE(f$converged))
  expect_length(f$rho, 2L)
})

test_that("maxeval is honoured, counted and reported as the stop reason", {
  f <- fit_matern_copula(Yiso(), c(12L, 9L), 1L, isotropic = TRUE, maxeval = 5L)
  expect_lte(f$n_evals, 5L)
  expect_false(f$converged)
  expect_identical(f$stop_reason, "maxeval")
  expect_identical(f$convergence, 5L)
  expect_true(is.finite(f$rho))
})

test_that("there is no wall-clock cap in the interface", {
  expect_false("maxtime" %in% names(formals(fit_matern_copula)))
})

test_that("a warm start is used and reduces the evaluation count", {
  cold <- fit_matern_copula(Yiso(), c(12L, 9L), 1L, isotropic = TRUE)
  warm <- fit_matern_copula(Yiso(), c(12L, 9L), 1L, isotropic = TRUE, start = cold$rho)
  expect_lt(abs(warm$rho - cold$rho), 1e-6)
  expect_lt(warm$n_evals, cold$n_evals)
})

test_that("the kronecker objective agrees with the eigenpair objective and the dense reference", {
  set.seed(4)
  X <- matrix(rnorm(12 * 9 * 5), 12 * 9, 5)
  a <- dmatern_copula(X, c(12L, 9L), 0.6, 1L)
  b <- dmatern_copula(X, c(12L, 9L), 0.6, 1L, method = "kronecker")
  expect_lt(max(abs(a - b)), 1e-9)
  X5 <- matrix(rnorm(20 * 4), 20, 4)
  for (nu in 0:2) {
    ref <- ref_copula_logdens(X5, ref_Qtilde(5L, 4L, 0.7, 0.4, nu))
    expect_lt(max(abs(dmatern_copula(X5, c(5L, 4L), c(0.7, 0.4), nu, method = "kronecker") - ref)), 1e-9,
      label = paste0("nu = ", nu))
  }
})

test_that("fit_matern_copula(method = 'kronecker') reproduces the default fit", {
  Y <- Yiso()
  a <- fit_matern_copula(Y, c(12L, 9L), 1L, isotropic = TRUE)
  b <- fit_matern_copula(Y, c(12L, 9L), 1L, isotropic = TRUE, method = "kronecker")
  expect_lt(abs(a$rho - b$rho), 1e-8)
  # the objectives agree to ~1e-10, which can move BOBYQA's stopping test by a few evaluations
  # (22 vs 24 seen in an optimised build); the estimate, not the path, is the contract
  expect_lte(abs(a$n_evals - b$n_evals), 6L)
  a2 <- fit_matern_copula(Y, c(12L, 9L), 1L, isotropic = FALSE)
  b2 <- fit_matern_copula(Y, c(12L, 9L), 1L, isotropic = FALSE, method = "kronecker")
  expect_lt(max(abs(a2$rho - b2$rho)), 1e-8)
})

test_that("the default method is the eigenpair objective", {
  expect_identical(eval(formals(fit_matern_copula)$method)[1], "eigenpair")
})
