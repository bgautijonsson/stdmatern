# The three production paths must be bit-for-bit reproducible regardless of the OpenMP thread count
# the process was launched with. Before 2026-09-02 the samplers, the precision builder and the density
# ran OpenMP regions with a critical-section merge, so results differed run to run at 4 or 10 threads
# and depended on OMP_NUM_THREADS -- which libgomp reads once, at load; Sys.setenv() inside a session
# cannot change it, hence one fresh subprocess per setting.
skip_if_not_installed("callr")
skip_if_not_installed("pkgload")

pkg_root <- system.file(package = "stdmatern")
dev_tree <- dir.exists(file.path(pkg_root, "src"))

run_with_threads <- function(omp) {
  env <- c(callr::rcmd_safe_env(), VECLIB_MAXIMUM_THREADS = "1")
  if (!is.null(omp)) env <- c(env, OMP_NUM_THREADS = omp)
  callr::r(
    function(root, dev, libs) {
      if (dev) pkgload::load_all(root, quiet = TRUE) else library(stdmatern, lib.loc = libs)
      set.seed(11)
      list(
        draw = rmatern_copula_eigen(4L, 60L, 50L, 0.8, 0.8, 1L, seed = 7L),
        draw_setseed = rmatern_copula_eigen(2L, 60L, 50L, 0.8, 0.8, 1L),
        Q = make_standardized_matern(c(60L, 50L), 0.8, 1L)@x,
        dens = dmatern_copula_eigen(matrix(sin(seq_len(3000 * 3)), 3000, 3), 60L, 50L, 0.7, 0.7, 1L)
      )
    },
    args = list(root = pkg_root, dev = dev_tree, libs = .libPaths()),
    env = env
  )
}

test_that("draw, precision and density are identical at 1, 4 and the default thread count", {
  one <- run_with_threads("1")
  four <- run_with_threads("4")
  dflt <- run_with_threads(NULL)
  for (nm in names(one)) {
    expect_identical(four[[nm]], one[[nm]], label = paste0(nm, " at 4 threads vs 1"))
    expect_identical(dflt[[nm]], one[[nm]], label = paste0(nm, " at the default thread count vs 1"))
  }
})

test_that("two launches at the default thread count agree bit for bit", {
  a <- run_with_threads(NULL)
  b <- run_with_threads(NULL)
  expect_identical(a, b)
})
