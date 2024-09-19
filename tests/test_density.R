library(tidyverse)
library(stdmatern)
dim1 <- 40
dim2 <- 40
rho1 <- 0.5
rho2 <- 0.3
nu <- 2
n <- 10

Z <- rmatern_copula_eigen(n, dim1, dim2, rho1, rho2, nu)

tibble(
  eigen = dmatern_copula_eigen(Z, dim1, dim2, rho1, rho2, nu),
  cholesky = dmatern_copula_cholesky(Z, dim1, dim2, rho1, rho2, nu),
  circulant = dmatern_copula_circulant(Z, dim1, dim2, rho1, rho2, nu),
  folded = dmatern_copula_folded(Z, dim1, dim2, rho1, rho2, nu),
  eigen_unsc = dmatern_eigen(Z, dim1, dim2, rho1, rho2, nu),
  cholesky_unsc = dmatern_cholesky(Z, dim1, dim2, rho1, rho2, nu)
) |>
  cor() |>
  round(3)
