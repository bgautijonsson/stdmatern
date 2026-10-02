# The Student-t copula (0.6.0): gaussian_to_t, dtmatern_copula, rtmatern_copula, fit_tmatern_copula.
# Dense references in helper-reference.R (ref_t_copula_logdens, ref_t_to_gaussian).

# Run `code` with the caller's RNG kind and .Random.seed restored afterwards.
local_rng_guard <- function(env = parent.frame()) {
  withr::local_preserve_seed(.local_envir = env)
  old <- RNGkind()
  withr::defer(suppressWarnings(RNGkind(old[1], old[2], old[3])), envir = env, priority = "last")
}

drop_w <- function(x) {
  attr(x, "W") <- NULL
  x
}

# ---- gaussian_to_t --------------------------------------------------------------------------------

test_that("gaussian_to_t keeps dimensions, is exactly odd and matches qt(pnorm(z)) in the body", {
  local_rng_guard()
  set.seed(1)
  Z <- matrix(rnorm(30), 6, 5)
  X <- gaussian_to_t(Z, 4)
  expect_identical(dim(X), dim(Z))
  expect_identical(gaussian_to_t(-Z, 4), -X)
  expect_equal(X, qt(pnorm(Z), 4), tolerance = 1e-12)
  expect_identical(gaussian_to_t(0, 3), 0)
})

test_that("gaussian_to_t keeps precision in both tails", {
  x8 <- gaussian_to_t(8, 5)
  expect_true(is.finite(x8) && x8 > 0)
  expect_identical(x8, -gaussian_to_t(-8, 5))
  tails <- gaussian_to_t(c(10, 20, 30), 3)
  expect_true(all(is.finite(tails)) && all(diff(tails) > 0))
  # the upper tail of the t matches the Gaussian tail it came from
  expect_equal(pt(x8, 5, lower.tail = FALSE, log.p = TRUE), pnorm(8, lower.tail = FALSE, log.p = TRUE),
    tolerance = 1e-10
  )
  # df -> infinity returns the Gaussian scores
  expect_equal(gaussian_to_t(c(-2, 0.3, 1.7), 1e8), c(-2, 0.3, 1.7), tolerance = 1e-6)
})

# ---- the density ----------------------------------------------------------------------------------

test_that("dtmatern_copula matches the dense reference on 5 x 4 grids", {
  local_rng_guard()
  set.seed(2)
  for (df in c(3, 12)) {
    X <- matrix(rt(20 * 4, df), 20, 4)
    for (nu in 0:2) {
      ref <- ref_t_copula_logdens(X, ref_Qtilde(5L, 4L, 0.7, 0.4, nu), df)
      got <- dtmatern_copula(X, c(5L, 4L), c(0.7, 0.4), nu, df = df)
      expect_lt(max(abs(got - ref)), 1e-10, label = sprintf("df = %g, nu = %d", df, nu))
    }
  }
})

test_that("dtmatern_copula tends to dmatern_copula as df grows", {
  Y <- rmatern_copula(4, c(5L, 4L), c(0.7, 0.4), 1L, seed = 2L)
  X <- gaussian_to_t(Y, 1e8)
  got <- dtmatern_copula(X, c(5L, 4L), c(0.7, 0.4), 1L, df = 1e8)
  expect_lt(max(abs(got - dmatern_copula(Y, c(5L, 4L), c(0.7, 0.4), 1L))), 1e-5)
})

test_that("dtmatern_copula's kronecker and eigenpair forms agree", {
  X <- rtmatern_copula(5, c(12L, 9L), 0.6, 1L, df = 4, seed = 3L)
  a <- dtmatern_copula(X, c(12L, 9L), 0.6, 1L, df = 4, method = "eigenpair")
  b <- dtmatern_copula(X, c(12L, 9L), 0.6, 1L, df = 4, method = "kronecker")
  expect_lt(max(abs(a - b)), 1e-9)
})

# ---- the sampler ----------------------------------------------------------------------------------

test_that("a seeded draw is bit-identical and set.seed governs seed = NULL", {
  local_rng_guard()
  a <- rtmatern_copula(5, c(6L, 5L), 0.7, 1L, df = 4, seed = 11L)
  b <- rtmatern_copula(5, c(6L, 5L), 0.7, 1L, df = 4, seed = 11L)
  expect_identical(a, b)
  set.seed(3)
  c1 <- rtmatern_copula(5, c(6L, 5L), 0.7, 1L, df = 4)
  set.seed(3)
  c2 <- rtmatern_copula(5, c(6L, 5L), 0.7, 1L, df = 4)
  expect_identical(c1, c2)
  # seed = NULL: Z from rmatern_copula(seed = NULL), then W from the session stream after it
  set.seed(3)
  z <- rmatern_copula(5, c(6L, 5L), 0.7, 1L)
  w <- rgamma(5, shape = 2, rate = 2)
  expect_identical(attr(c1, "W"), w)
  expect_equal(drop_w(c1) * rep(sqrt(w), each = 30L), z, tolerance = 1e-14)
})

test_that("a seeded draw leaves the caller's RNG state and kind untouched", {
  local_rng_guard()
  suppressWarnings(RNGkind("L'Ecuyer-CMRG", "Box-Muller", "Rejection"))
  set.seed(99)
  before <- get(".Random.seed", envir = globalenv())
  kind <- RNGkind()
  x <- rtmatern_copula(3, c(6L, 5L), 0.7, 1L, df = 4, seed = 5L)
  expect_identical(get(".Random.seed", envir = globalenv()), before)
  expect_identical(RNGkind(), kind)
  # and a missing .Random.seed stays missing
  rm(".Random.seed", envir = globalenv())
  y <- rtmatern_copula(3, c(6L, 5L), 0.7, 1L, df = 4, seed = 5L)
  expect_false(exists(".Random.seed", envir = globalenv(), inherits = FALSE))
  expect_identical(x, y)
})

test_that("the seeded W is rgamma after set.seed(seed) under Mersenne-Twister", {
  local_rng_guard()
  x <- rtmatern_copula(7, c(6L, 5L), 0.7, 1L, df = 6, seed = 17L)
  set.seed(17L, kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
  expect_identical(attr(x, "W"), rgamma(7, shape = 3, rate = 3))
})

test_that("attr W is a positive vector of length n and X * sqrt(W) is the rmatern_copula draw", {
  for (method in c("kronecker", "eigenpair")) {
    x <- rtmatern_copula(6, c(6L, 5L), 0.7, 1L, df = 4, seed = 8L, method = method)
    W <- attr(x, "W")
    expect_true(is.numeric(W) && length(W) == 6L && all(W > 0))
    expect_identical(dim(x), c(30L, 6L))
    z <- rmatern_copula(6, c(6L, 5L), 0.7, 1L, seed = 8L, method = method)
    expect_equal(drop_w(x) * rep(sqrt(W), each = 30L), z, tolerance = 1e-14)
  }
  expect_identical(dim(rtmatern_copula(0, c(6L, 5L), 0.7, 1L, df = 4, seed = 1L)), c(30L, 0L))
})

test_that("W is recovered from x by E[W | x] = (df + P) / (df + q) on a 40 x 50 grid", {
  df <- 5
  x <- rtmatern_copula(6, c(40L, 50L), 0.7, 1L, df = df, seed = 7L)
  sp <- matern_spectral(c(40L, 50L), 0.7, 1L)
  q <- colSums(x * as.matrix(sp$Q %*% x))
  ew <- (df + 2000) / (df + q)
  expect_lt(max(abs(ew / attr(x, "W") - 1)), 0.1)
})

# ---- the fit --------------------------------------------------------------------------------------

t_fit_data <- function() {
  x <- rtmatern_copula(50, 20L, 0.7, 1L, df = 5, seed = 2026L)
  ref_t_to_gaussian(drop_w(x), 5)
}

test_that("fit_tmatern_copula recovers rho and df", {
  Y <- t_fit_data()
  f <- fit_tmatern_copula(Y, 20L, 1L, isotropic = TRUE)
  expect_lt(abs(f$rho - 0.7), 0.03)
  expect_gt(f$df, 3.5)
  expect_lt(f$df, 7.5)
  expect_true(isTRUE(f$converged))
  expect_identical(f$stop_reason, "converged")
  expect_true(is.integer(f$n_evals) && f$n_evals > 0L)
  expect_false(f$df_fixed)
  expect_false(f$df_at_bound)
  expect_equal(f$loglik, sum(dtmatern_copula(gaussian_to_t(Y, f$df), 20L, f$rho, 1L, f$df)), tolerance = 1e-10)
})

test_that("a fixed df is returned exactly and only rho is searched", {
  Y <- t_fit_data()
  f <- fit_tmatern_copula(Y, 20L, 1L, df = 5, isotropic = TRUE)
  expect_identical(f$df, 5)
  expect_true(f$df_fixed)
  expect_false(f$df_at_bound)
  expect_lt(abs(f$rho - 0.7), 0.03)
  fa <- fit_tmatern_copula(Y, 20L, 1L, df = 5, isotropic = FALSE)
  expect_length(fa$rho, 2L)
  expect_identical(fa$rho, c(fa$rho1, fa$rho2))
  expect_true(isTRUE(fa$converged))
})

test_that("Gaussian data push df to the upper bound", {
  Y <- rmatern_copula(50, 20L, 0.7, 1L, seed = 2026L)
  f <- fit_tmatern_copula(Y, 20L, 1L, df_bounds = c(1, 50), isotropic = TRUE)
  expect_true(f$df_at_bound)
  expect_equal(f$df, 50, tolerance = 1e-6)
})

test_that("maxeval is honoured and reported", {
  f <- fit_tmatern_copula(t_fit_data(), 20L, 1L, isotropic = TRUE, maxeval = 5L)
  expect_lte(f$n_evals, 5L)
  expect_false(f$converged)
  expect_identical(f$stop_reason, "maxeval")
  expect_identical(f$convergence, 5L)
})

# ---- validation -----------------------------------------------------------------------------------

test_that("bad df, df_bounds and df_start are refused at the R boundary", {
  X <- matrix(0.1, 20, 2)
  for (bad in list(0, -1, NA_real_, Inf, c(3, 4), "5", NULL)) {
    expect_error(dtmatern_copula(X, c(5L, 4L), 0.5, 1L, df = bad), "df must be")
    expect_error(gaussian_to_t(X, bad), "df must be")
    expect_error(rtmatern_copula(2, c(5L, 4L), 0.5, 1L, df = bad, seed = 1L), "df must be")
  }
  expect_error(gaussian_to_t("a", 3), "Z must be numeric")
  expect_error(fit_tmatern_copula(X, c(5L, 4L), 1L, df = -2), "df must be")
  for (bad in list(c(5, 1), c(0, 10), c(1, Inf), c(-1, 3), c(2, 2), 3, c(1, 2, 3), c(NA, 4))) {
    expect_error(fit_tmatern_copula(X, c(5L, 4L), 1L, df_bounds = bad), "df_bounds must be")
  }
  expect_error(fit_tmatern_copula(X, c(5L, 4L), 1L, df_start = 0.5), "df_start must lie within")
  expect_error(fit_tmatern_copula(X, c(5L, 4L), 1L, df_start = 600), "df_start must lie within")
  expect_error(fit_tmatern_copula(X, c(5L, 4L), 1L, df_start = 0), "df_start must be")
  expect_error(fit_tmatern_copula(X, c(5L, 4L), 1L, start = 1.2), "start must lie")
})

test_that("review fixes: stable constant at huge df, no overflow of q, failed fits say so, W underflow refused", {
  Y <- rmatern_copula(3, c(12, 10), 0.7, 1, seed = 2)
  for (dfv in c(1e8, 1e12, 1e16)) {
    expect_lt(
      max(abs(dtmatern_copula(gaussian_to_t(Y, dfv), c(12, 10), 0.7, 1, dfv) -
        dmatern_copula(Y, c(12, 10), 0.7, 1))),
      1e-5
    )
  }
  X1 <- gaussian_to_t(c(-30, rep(0.1, 29)), 1)
  expect_true(is.finite(X1[1]))
  expect_true(is.finite(dtmatern_copula(X1, c(6, 5), 0.5, 1, df = 1)))
  Yb <- rmatern_copula(10, c(6, 5), 0.5, 1, seed = 4)
  Yb[1, 1] <- NA
  expect_error(fit_tmatern_copula(Yb, c(6, 5), 1, df = 5, isotropic = TRUE), "finite")
  expect_error(rtmatern_copula(500, c(6, 5), 0.7, 1, df = 0.01, seed = 1), "too small")
})
