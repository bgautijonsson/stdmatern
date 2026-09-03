test_that("marginal_sd_eigen returns correct length", {
  Q1 <- make_AR_prec_matrix(5L, 0.5)
  Q2 <- make_AR_prec_matrix(4L, 0.5)
  solver1 <- eigen(as.matrix(Q1), symmetric = TRUE)
  solver2 <- eigen(as.matrix(Q2), symmetric = TRUE)

  sds <- stdmatern:::marginal_sd_eigen(
    solver1$values, solver1$vectors, 5L,
    solver2$values, solver2$vectors, 4L,
    0L
  )
  expect_length(sds, 20)
})

test_that("marginal_sd_eigen returns positive values", {
  Q1 <- make_AR_prec_matrix(5L, 0.5)
  Q2 <- make_AR_prec_matrix(5L, 0.5)
  solver1 <- eigen(as.matrix(Q1), symmetric = TRUE)
  solver2 <- eigen(as.matrix(Q2), symmetric = TRUE)

  sds <- stdmatern:::marginal_sd_eigen(
    solver1$values, solver1$vectors, 5L,
    solver2$values, solver2$vectors, 5L,
    0L
  )
  expect_true(all(sds > 0))
})

test_that("marginal SDs are consistent with precision matrix inverse", {
  dim_x <- 4L
  dim_y <- 4L
  rho <- 0.5
  nu <- 0L

  Q <- make_standardized_matern(dim = c(dim_x, dim_y), rho = rho, nu = nu)
  Sigma <- solve(as.matrix(Q))
  marginal_sds_from_Q <- sqrt(diag(Sigma))

  # After standardization, marginal SDs should be 1
  expect_equal(marginal_sds_from_Q, rep(1, dim_x * dim_y), tolerance = 1e-8)
})
