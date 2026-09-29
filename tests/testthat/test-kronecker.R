# The Kronecker forms (perf-backend, 2026-09-29): marginal sds as S2 = W1 M W2', the
# sampler as V2 (lambda^(-(nu+1)/2) o E) V1' per draw. Same mathematics as the
# eigenpair forms, O(N (n1 + n2)) instead of O(N^2); agreement is to rounding.

grids <- list(c(6L, 5L), c(1L, 7L), c(9L, 1L), c(13L, 11L))
cases <- expand.grid(g = seq_along(grids), nu = 0:2, aniso = c(FALSE, TRUE))

test_that("Kronecker marginal sds match the eigenpair sds and the dense inverse", {
  for (r in seq_len(nrow(cases))) {
    d <- grids[[cases$g[r]]]
    rho <- if (cases$aniso[r]) c(0.83, 0.41) else 0.7
    nu <- cases$nu[r]
    Qe <- make_standardized_matern(d, rho, nu)
    Qk <- make_standardized_matern(d, rho, nu, method = "kronecker")
    lab <- sprintf("dims %s nu %d rho %s", paste(d, collapse = "x"), nu, paste(rho, collapse = ","))
    expect_lt(max(abs(Qk - Qe)) / max(abs(Qe)), 1e-13, label = lab)
    expect_equal(Qk@i, Qe@i, label = lab)
    expect_equal(Qk@p, Qe@p, label = lab)
    expect_lt(max(abs(diag(solve(as.matrix(Qk))) - 1)), 1e-10, label = lab)
  }
})

test_that("matern_spectral agrees across methods, half log-determinant included", {
  for (nu in 0:2) {
    se <- matern_spectral(c(12L, 9L), c(0.9, 0.6), nu)
    sk <- matern_spectral(c(12L, 9L), c(0.9, 0.6), nu, method = "kronecker")
    expect_lt(abs(sk$hld - se$hld), 1e-10 * abs(se$hld))
    expect_lt(max(abs(sk$sd - se$sd)), 1e-14)
    dense <- 0.5 * as.numeric(determinant(as.matrix(sk$Q), logarithm = TRUE)$modulus)
    expect_lt(abs(sk$hld - dense), 1e-9 * abs(dense))
  }
})

test_that("the Kronecker sampler reproduces the eigenpair sampler's draws from the same seed", {
  for (r in seq_len(nrow(cases))) {
    d <- grids[[cases$g[r]]]
    rho <- if (cases$aniso[r]) c(0.83, 0.41) else 0.7
    nu <- cases$nu[r]
    Ze <- rmatern_copula(7L, d, rho, nu, seed = 11L)
    Zk <- rmatern_copula(7L, d, rho, nu, seed = 11L, method = "kronecker")
    expect_equal(dim(Zk), dim(Ze))
    expect_lt(max(abs(Zk - Ze)), 1e-12 * max(1, max(abs(Ze))),
      label = sprintf("dims %s nu %d", paste(d, collapse = "x"), nu)
    )
  }
})

test_that("the Kronecker sampler is reproducible, handles n = 0, and follows set.seed without a seed", {
  a <- rmatern_copula(4L, c(8L, 6L), 0.8, 1L, seed = 3L, method = "kronecker")
  expect_identical(a, rmatern_copula(4L, c(8L, 6L), 0.8, 1L, seed = 3L, method = "kronecker"))
  expect_equal(dim(rmatern_copula(0L, c(8L, 6L), 0.8, 1L, seed = 3L, method = "kronecker")), c(48L, 0L))
  set.seed(9)
  b1 <- rmatern_copula(2L, c(5L, 5L), 0.5, 0L, method = "kronecker")
  set.seed(9)
  b2 <- rmatern_copula(2L, c(5L, 5L), 0.5, 0L, method = "kronecker")
  expect_identical(b1, b2)
})

test_that("the Kronecker sampler's linear map has covariance Q~^-1 exactly", {
  d <- c(7L, 5L)
  N <- prod(d)
  sp <- matern_spectral(d, c(0.75, 0.55), 1L, method = "kronecker")
  lam <- outer(sp$e1$values, sp$e2$values, "+")^2
  L <- vapply(seq_len(N), function(k) {
    e <- numeric(N)
    e[k] <- 1
    E <- matrix(e, d[2], d[1]) / sqrt(t(lam))
    as.vector(sp$e2$vectors %*% E %*% t(sp$e1$vectors)) / sp$sd
  }, numeric(N))
  expect_lt(max(abs(tcrossprod(L) - solve(as.matrix(sp$Q)))), 1e-12)
})
