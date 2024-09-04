library(stdmatern)

n <- 1
dim_x <- 256
dim_y <- 256
rho1 <- 0.5
rho2 <- 0.5
nu <- 2




bench::mark(
  "copula" = rmatern_copula_eigen(
    n = n,
    dim_x = dim_x,
    dim_y = dim_y,
    rho1 = rho1,
    rho2 = rho2,
    nu = nu
  ),
  "unscaled" = rmatern_eigen(
    n = n,
    dim_x = dim_x,
    dim_y = dim_y,
    rho1 = rho1,
    rho2 = rho2,
    nu = nu
  ),
  "circulant" = rmatern_copula_circulant(
    n = n,
    dim1 = dim_x,
    dim2 = dim_y,
    rho1 = rho1,
    rho2 = rho2,
    nu = nu
  ),
  "folded" = rmatern_copula_folded_full(
    n = n,
    dim1 = dim_x,
    dim2 = dim_y,
    rho1 = rho1,
    rho2 = rho2,
    nu = nu
  ),
  iterations = 10,
  filter_gc = FALSE,
  check = FALSE
)
