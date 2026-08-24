# Reproducibility of the samplers (added 2026-08-24).
#
# Before the seed argument existed, every sampler seeded std::mt19937 from
# std::random_device, so set.seed() could not reach it and repeated calls
# differed. The eigen samplers additionally shared one std::normal_distribution
# across OpenMP threads, making the draw depend on the thread count.

test_that("an explicit seed reproduces the draw", {
  a <- rmatern_copula_eigen(5L, 6L, 5L, 0.5, 0.5, 1L, seed = 42L)
  b <- rmatern_copula_eigen(5L, 6L, 5L, 0.5, 0.5, 1L, seed = 42L)
  expect_identical(a, b)
  c <- rmatern_copula_eigen(5L, 6L, 5L, 0.5, 0.5, 1L, seed = 43L)
  expect_false(isTRUE(all.equal(a, c)))
})

test_that("set.seed() governs the draw when no seed is passed", {
  set.seed(99)
  a <- rmatern_copula_eigen(5L, 6L, 5L, 0.5, 0.5, 1L)
  set.seed(99)
  b <- rmatern_copula_eigen(5L, 6L, 5L, 0.5, 0.5, 1L)
  expect_identical(a, b)
})

test_that("the draw does not depend on the OpenMP thread count", {
  old <- Sys.getenv("OMP_NUM_THREADS", unset = NA)
  on.exit(if (is.na(old)) Sys.unsetenv("OMP_NUM_THREADS") else
    Sys.setenv(OMP_NUM_THREADS = old), add = TRUE)
  Sys.setenv(OMP_NUM_THREADS = "1")
  a <- rmatern_copula_eigen(4L, 6L, 5L, 0.6, 0.6, 1L, seed = 7L)
  Sys.setenv(OMP_NUM_THREADS = "4")
  b <- rmatern_copula_eigen(4L, 6L, 5L, 0.6, 0.6, 1L, seed = 7L)
  expect_identical(a, b)
})

test_that("the other samplers take a seed too", {
  expect_identical(rmatern_eigen(3L, 5L, 4L, 0.5, 0.5, 1L, seed = 11L),
                   rmatern_eigen(3L, 5L, 4L, 0.5, 0.5, 1L, seed = 11L))
  expect_identical(rmatern_copula_circulant(3L, 6L, 6L, 0.5, 0.5, 1L, seed = 11L),
                   rmatern_copula_circulant(3L, 6L, 6L, 0.5, 0.5, 1L, seed = 11L))
  expect_identical(rmatern_copula_folded_full(3L, 6L, 6L, 0.5, 0.5, 1L, seed = 11L),
                   rmatern_copula_folded_full(3L, 6L, 6L, 0.5, 0.5, 1L, seed = 11L))
})

test_that("the R-level wrappers pass the seed through", {
  expect_identical(rmatern_copula(3L, c(5L, 4L), 0.5, 1L, seed = 5L),
                   rmatern_copula(3L, c(5L, 4L), 0.5, 1L, seed = 5L))
})
