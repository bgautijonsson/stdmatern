# Agreement with the independent dense reference (helper-reference.R) on a 5 x 4 grid.
n1 <- 5L
n2 <- 4L
rho1 <- 0.7
rho2 <- 0.4

test_that("make_standardized_matern matches the dense reference at nu = 0, 1, 2", {
  for (nu in 0:2) {
    Qt <- as.matrix(make_standardized_matern(c(n1, n2), c(rho1, rho2), nu))
    ref <- ref_Qtilde(n1, n2, rho1, rho2, nu)
    expect_lt(max(abs(Qt - ref)) / max(abs(ref)), 1e-12, label = paste0("nu = ", nu))
  }
})

test_that("station order is x-major, k = (a - 1) * n2 + b", {
  Qt <- as.matrix(make_standardized_matern(c(n1, n2), c(rho1, rho2), 0L))
  ref <- ref_Qtilde(n1, n2, rho1, rho2, 0L)
  k <- (2 - 1) * n2 + 2 # cell (a = 2, b = 2)
  expect_true(Qt[k, k + 1] != 0) # y-neighbour (a, b + 1)
  expect_true(Qt[k, k + n2] != 0) # x-neighbour (a + 1, b)
  expect_equal(Qt[k, k + 2], 0)
  expect_equal(Qt[k, k + 1], ref[k, k + 1], tolerance = 1e-12)
  expect_equal(Qt[k, k + n2], ref[k, k + n2], tolerance = 1e-12)
})

test_that("marginal_sd_eigen matches sqrt(diag(solve(Q)))", {
  for (nu in 0:2) {
    e1 <- eigen(ref_ar1_prec(n1, rho1), symmetric = TRUE)
    e2 <- eigen(ref_ar1_prec(n2, rho2), symmetric = TRUE)
    sd <- stdmatern:::marginal_sd_eigen(e1$values, e1$vectors, n1, e2$values, e2$vectors, n2, nu)
    expect_lt(max(abs(sd - ref_marginal_sd(n1, n2, rho1, rho2, nu))), 1e-12, label = paste0("nu = ", nu))
  }
})

test_that("dmatern_copula matches the dense reference at nu = 0, 1, 2", {
  set.seed(7)
  X <- matrix(rnorm(n1 * n2 * 6), n1 * n2, 6)
  for (nu in 0:2) {
    ld <- dmatern_copula(X, c(n1, n2), c(rho1, rho2), nu)
    ref <- ref_copula_logdens(X, ref_Qtilde(n1, n2, rho1, rho2, nu))
    expect_lt(max(abs(ld - ref)), 1e-9, label = paste0("nu = ", nu))
  }
})

test_that("dmatern (unscaled) matches the dense reference", {
  set.seed(8)
  X <- matrix(rnorm(n1 * n2 * 4), n1 * n2, 4)
  for (nu in 0:2) {
    ld <- dmatern(X, c(n1, n2), c(rho1, rho2), nu)
    ref <- ref_normal_logdens(X, ref_Q_power(n1, n2, rho1, rho2, nu))
    expect_lt(max(abs(ld - ref)), 1e-9, label = paste0("nu = ", nu))
  }
})

# The Cholesky variant was verified against this reference after its F3 fix and then archived
# (archive/retired-2026-09-02/), so no block for it remains here.

test_that("rmatern_copula_eigen draws have the covariance solve(Qtilde) and unit margins", {
  nS <- 50000L
  Z <- rmatern_copula_eigen(nS, n1, n2, rho1, rho2, 1L, seed = 99L)
  Sig <- solve(ref_Qtilde(n1, n2, rho1, rho2, 1L))
  # MC sd of a covariance entry ~ 1/sqrt(nS) = 0.0045; 0.03 is ~6 sd, far below any structural error
  expect_lt(max(abs(cov(t(Z)) - Sig)), 0.03)
  expect_lt(max(abs(apply(Z, 1, var) - 1)), 0.05)
  expect_lt(max(abs(rowMeans(Z))), 0.03)
})
