library(stdmatern)

bench::mark(
  "copula" = rmatern_copula_eigen(
    n = 50,
    dim_x = 100,
    dim_y = 100,
    rho1 = 0.5,
    rho2 = 0.5,
    nu = 2
  ),
  "unscaled" = rmatern_eigen(
    n = 50,
    dim_x = 100,
    dim_y = 100,
    rho1 = 0.5,
    rho2 = 0.5,
    nu = 2
  ),
  "circulant" = rmatern_copula_circulant(
    n = 50,
    dim1 = 100,
    dim2 = 100,
    rho1 = 0.5,
    rho2 = 0.5,
    nu = 2
  ),
  "folded" = rmatern_copula_folded_full(
    n = 50,
    dim1 = 100,
    dim2 = 100,
    rho1 = 0.5,
    rho2 = 0.5,
    nu = 2
  ),
  iterations = 10,
  filter_gc = FALSE,
  check = FALSE
)
