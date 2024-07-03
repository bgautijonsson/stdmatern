library(stdmatern)
library(Matrix)
library(purrr)
library(tidyverse)
library(glue)


my_fun <- function(dim) {
  x <- rnorm(dim^2)
  bench::mark(
    matern_mvn_density(x, dim, 0.5, 0),
    filter_gc = FALSE,
    iterations = 10,
    check = FALSE
  ) |> 
    mutate(
      dim = dim
    )
}


results <- map(c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100), my_fun)

results |> 
  list_rbind() |> 
  select(Q_size = dim, time = median, memory = mem_alloc) |> 
  mutate(
    Q_size = glue("{Q_size^2}x{Q_size^2}")
  ) |> 
  select(-memory)
