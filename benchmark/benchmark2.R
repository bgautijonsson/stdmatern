library(dplyr)
library(sparseMVN)
library(ggplot2)
library(stdmatern)
library(purrr)
library(patchwork)




f1 <- function() {
  C <- matern_mvn_density(X, grid_dim, rho, nu) |> 
    sum()
  return(1)
}

f2 <- function() {
  R <- dmvn.sparse(
    X |> t(),
    mu = rep(0, grid_dim^2),
    CH = chol_Q
  ) |> 
    sum()
  return(1)
}

grid_dim <- 20
n_replicates <- 20
rho <- 0.5
nu <- 0
Q <- make_standardized_matern_eigen(grid_dim, rho, nu)
chol_Q <- Cholesky(Q)

X <- rmvn.sparse(
  n = n_replicates,
  mu = rep(0, nrow(Q)),
  CH = chol_Q
) |> t()

bench::mark(
  "Eigen" = f1(),
  "Basic" = f2()
) |> 
  mutate(
    relative = as.numeric(median / min(median)),
    .before = min
  )


