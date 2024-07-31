library(stdmatern)
library(Matrix)
library(purrr)
library(tidyverse)
library(glue)

my_fun <- function(dim) {
  x_copula <- rmatern_copula(1, dim, 0.5, 0)
  x <- rmatern(1, dim, 0.5, 0)
  bench::mark(
    "Copula" = dmatern_copula(x, dim, 0.5, 0),
    "Regular" = dmatern(x, dim, 0.5, 0),
    filter_gc = FALSE,
    iterations = 10,
    check = FALSE
  ) |> 
    mutate(
      dim = dim
    )
}

my_fun(200)



results <- map(seq(10, 200, by = 10), my_fun)

results |> 
  list_rbind() |> 
  select(Q_size = dim, Type = expression, time = median, memory = mem_alloc) |> 
  mutate(
    Q_size = glue("{Q_size^2}x{Q_size^2}")
  ) |> 
  select(-memory) |> 
  pivot_wider(names_from = Type, values_from = time)




results |> 
  list_rbind() |> 
  select(dim, type = expression, time = median, memory = mem_alloc) |> 
  mutate(
    locations = dim^2,
    type = as.character(type)
  ) |> 
  ggplot(aes(locations, time)) +
  geom_line(aes(lty = type))

