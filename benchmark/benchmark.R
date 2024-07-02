library(cholesky)
library(Matrix)
library(purrr)
library(tidyverse)
library(glue)

dim <- 30
bench::mark(
  make_standardized_matern(dim, 0.5),
  filter_gc = FALSE,
  iterations = 10,
  check = FALSE
)

my_fun <- function(dim) {
  bench::mark(
    make_standardized_matern(dim, 0.5),
    filter_gc = FALSE,
    iterations = 10,
    check = FALSE
  ) |> 
    mutate(
      dim = dim
    )
}


results <- map(c(10, 20, 40, 80), my_fun)

results |> 
  list_rbind() |> 
  select(Q_size = dim, time = median, memory = mem_alloc) |> 
  mutate(
    Q_size = glue("{Q_size^2}x{Q_size^2}")
  )
