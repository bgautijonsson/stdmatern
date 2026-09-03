# matern_spectral(): the standardised precision, its half log-determinant, the marginal SDs and the
# two AR(1) eigen systems from one call (2026-09-02). Qtilde must be the make_standardized_matern
# object bit for bit -- the S05/U04/U05 pins in maxandsmooth depend on it.
test_that("matern_spectral returns make_standardized_matern's Qtilde bit for bit", {
  for (nu in 0:2) {
    sp <- matern_spectral(c(9L, 7L), c(0.7, 0.4), nu)
    expect_identical(sp$Q, make_standardized_matern(c(9L, 7L), c(0.7, 0.4), nu))
    expect_identical(sp$x, sp$Q@x)
    expect_identical(sp$i, sp$Q@i)
    expect_identical(sp$p, sp$Q@p)
  }
})

test_that("matern_spectral's half log-determinant and marginal sd match the dense reference", {
  for (nu in 0:2) {
    sp <- matern_spectral(c(5L, 4L), c(0.7, 0.4), nu)
    Qt <- ref_Qtilde(5L, 4L, 0.7, 0.4, nu)
    expect_lt(abs(sp$hld - ref_half_logdet(Qt)) / abs(ref_half_logdet(Qt)), 1e-10, label = paste0("nu = ", nu))
    expect_lt(max(abs(sp$sd - ref_marginal_sd(5L, 4L, 0.7, 0.4, nu))), 1e-12, label = paste0("nu = ", nu))
  }
})

test_that("matern_spectral returns the two AR(1) eigen systems consistent with hld", {
  sp <- matern_spectral(c(5L, 4L), 0.6, 1L)
  expect_named(sp$e1, c("values", "vectors"))
  expect_length(sp$e1$values, 5L)
  expect_length(sp$e2$values, 4L)
  lam <- outer(sp$e1$values, sp$e2$values, `+`)^2
  expect_lt(abs(sp$hld - 0.5 * (sum(log(lam)) + 2 * sum(log(sp$sd)))) / abs(sp$hld), 1e-12)
  expect_identical(sp$dims, c(5L, 4L))
  expect_identical(sp$nu, 1L)
})
