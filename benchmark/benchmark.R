library(stdmatern)
library(Matrix)
library(purrr)
library(tidyverse)
library(glue)

x <- sample_standardized_matern(5, 0.5, 0, 1)

matern_mvn_density_cholesky(x, 5, 0.5, 0)
matern_mvn_density_eigen_whitened(x, 5, 0.5, 0)

my_fun <- function(dim) {
  x <- sample_standardized_matern(dim, 0.5, 0, 1)
  bench::mark(
    "Cholesky" = matern_mvn_density_cholesky(x, dim, 0.5, 0),
    "Eigen" = matern_mvn_density_eigen_whitened(x, dim, 0.5, 0),
    filter_gc = FALSE,
    iterations = 10,
    check = FALSE
  ) |> 
    mutate(
      dim = dim
    )
}

my_fun(100)


results <- map(c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 200), my_fun)

results |> 
  list_rbind() |> 
  select(Q_size = dim, time = median, memory = mem_alloc) |> 
  mutate(
    Q_size = glue("{Q_size^2}x{Q_size^2}")
  ) |> 
  select(-memory)
