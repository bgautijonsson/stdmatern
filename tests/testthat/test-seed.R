# Reproducibility of the samplers (added 2026-08-24; thread-count invariance moved to
# test-determinism.R on 2026-09-02, where it runs one subprocess per OMP_NUM_THREADS).
#
# Before the seed argument existed, every sampler seeded std::mt19937 from
# std::random_device, so set.seed() could not reach it and repeated calls differed.

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

test_that("the unstandardised sampler takes a seed too", {
  expect_identical(
    rmatern(3L, c(5L, 4L), 0.5, 1L, seed = 11L),
    rmatern(3L, c(5L, 4L), 0.5, 1L, seed = 11L)
  )
})

test_that("the R-level wrappers pass the seed through", {
  expect_identical(
    rmatern_copula(3L, c(5L, 4L), 0.5, 1L, seed = 5L),
    rmatern_copula(3L, c(5L, 4L), 0.5, 1L, seed = 5L)
  )
})
