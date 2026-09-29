# The seeding fix of 0.5.0 (2026-09-29): each eigenpair's generator is seeded from the full 64-bit
# key instead of its low 32 bits, which let two eigenpairs of one draw share a stream; "legacy32"
# reproduces the 0.4.0 draws.

test_that("the 32-bit legacy seeding lets two eigenpairs share a stream; the full seeding does not", {
  # found by search (tools: pair_first_normals_cpp over seeds 1..120 on 150 x 150): seeds 16, 41, 83
  first_legacy <- pair_first_normals_cpp(150L, 150L, 16L, TRUE)
  shared <- which(duplicated(first_legacy) | duplicated(first_legacy, fromLast = TRUE))
  expect_identical(shared, c(14772L, 21069L))
  for (s in c(16L, 41L, 83L, 1:20)) {
    expect_false(anyDuplicated(pair_first_normals_cpp(150L, 150L, s, FALSE)) > 0, label = paste("seed", s))
  }
})

test_that("seeding = \"legacy32\" reproduces the 0.4.0 draws exactly", {
  # values printed by stdmatern 0.4.0 on 2026-09-29 (17 significant digits)
  a <- rmatern_copula(3L, c(8L, 6L), c(0.7, 0.4), 1L, seed = 11L, seeding = "legacy32")
  expect_identical(
    sprintf("%.17g", c(sum(a), a[1], a[length(a)])),
    c("83.206532551887008", "0.93389256508829965", "0.34169429290692716")
  )
  b <- rmatern_copula(2L, c(5L, 7L), 0.5, 0L, seed = 2026L, method = "eigenpair", seeding = "legacy32")
  expect_identical(
    sprintf("%.17g", c(sum(b), b[1], b[length(b)])),
    c("2.8129298917786767", "0.81365147017147321", "-0.31478818364126104")
  )
  u <- rmatern(2L, c(6L, 5L), 0.6, 1L, seed = 5L, seeding = "legacy32")
  expect_identical(
    sprintf("%.17g", c(sum(u), u[1], u[length(u)])),
    c("-11.731285464176759", "-0.088721120524879918", "0.010019082253163501")
  )
})

test_that("the full seeding is the default, reproducible, and shared by both samplers", {
  expect_identical(eval(formals(rmatern_copula)$seeding), "full")
  a <- rmatern_copula(4L, c(9L, 7L), 0.6, 1L, seed = 3L)
  expect_identical(a, rmatern_copula(4L, c(9L, 7L), 0.6, 1L, seed = 3L, seeding = "full"))
  expect_false(isTRUE(all.equal(a, rmatern_copula(4L, c(9L, 7L), 0.6, 1L, seed = 3L, seeding = "legacy32"))))
  e <- rmatern_copula(4L, c(9L, 7L), 0.6, 1L, seed = 3L, method = "eigenpair")
  expect_lt(max(abs(a - e)), 1e-12)
})

test_that("the option and the environment variable choose the default; bad values are refused", {
  ref_legacy <- rmatern_copula(2L, c(6L, 5L), 0.5, 1L, seed = 9L, seeding = "legacy32")
  ref_full <- rmatern_copula(2L, c(6L, 5L), 0.5, 1L, seed = 9L, seeding = "full")
  withr::with_options(list(stdmatern.seeding = "legacy32"), {
    expect_identical(rmatern_copula(2L, c(6L, 5L), 0.5, 1L, seed = 9L), ref_legacy)
    expect_identical(rmatern_copula_eigen(2L, 6L, 5L, 0.5, 0.5, 1L, seed = 9L), ref_legacy)
  })
  withr::with_envvar(c(STDMATERN_SEEDING = "legacy32"), {
    expect_identical(rmatern_copula(2L, c(6L, 5L), 0.5, 1L, seed = 9L), ref_legacy)
    withr::with_options(list(stdmatern.seeding = "full"), {
      expect_identical(rmatern_copula(2L, c(6L, 5L), 0.5, 1L, seed = 9L), ref_full)
    })
  })
  expect_error(rmatern_copula(2L, c(6L, 5L), 0.5, 1L, seed = 9L, seeding = "legacy"), "seeding must be")
  expect_error(rmatern(2L, c(6L, 5L), 0.5, 1L, seed = 9L, seeding = NA), "seeding must be")
})
