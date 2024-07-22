library(dplyr)
library(sparseMVN)
library(ggplot2)
library(stdmatern)
library(purrr)
library(patchwork)

grid_dim <- 10
n_replicates <- 20
rho <- 0.5
nu <- 1
Q <- make_standardized_matern_eigen(grid_dim, rho, nu)
chol_Q <- Cholesky(Q)



my_fun <- function(iter) {
  
  X <- rmvn.sparse(
    n = n_replicates,
    mu = rep(0, nrow(Q)),
    CH = chol_Q
  ) |> t()

  C <- matern_mvn_density(X, grid_dim, rho, nu) |> 
    sum()

  R <- dmvn.sparse(
    X |> t(),
    mu = rep(0, grid_dim^2),
    CH = chol_Q
  ) |> 
    sum()


  tibble(
    C = C,
    R = R,
    iter = iter
  )
}

d <- map(seq_len(200), my_fun) |> 
  list_rbind()



p1 <- d |> 
  ggplot(aes(C, R)) +
  geom_point()



p2 <- d |> 
  mutate(
    diff = C / R
  ) |> 
  ggplot(aes(C, diff)) +
  geom_point()

p1 + p2

d |> 
  select(-iter) |> 
  cor()