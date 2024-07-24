library(dplyr)
library(sparseMVN)
library(ggplot2)
library(stdmatern)
library(purrr)
library(patchwork)

grid_dim <- 20
n_replicates <- 10
rho <- 0.5
nu <- 0



my_fun <- function(iter) {
  
  X <- sample_standardized_matern(grid_dim, rho, nu, n_replicates)

  Eigen <- matern_mvn_density_eigen(X, grid_dim, rho, nu) |> 
    sum()

  Choll <- matern_mvn_density_cholesky(X, grid_dim, rho, nu) |> 
    sum()


  tibble(
    Eigen = Eigen,
    Cholesky = Choll,
    iter = iter
  )
}

d <- map(seq_len(100), my_fun) |> 
  list_rbind()



p1 <- d |> 
  ggplot(aes(Eigen, Cholesky)) +
  geom_point()



p2 <- d |> 
  mutate(
    diff = Eigen / Cholesky
  ) |> 
  ggplot(aes(Eigen, diff)) +
  geom_point()

p1 + p2

d |> 
  select(-iter) |> 
  cor()
