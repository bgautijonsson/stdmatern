library(stdmatern)
library(tidyverse)
dim1 <- 20
dim2 <- 30
rho1 <- 0.9
rho2 <- 0.9
nu <- 2

Z <- rmatern_copula_eigen(1, dim1, dim2, rho1, rho2, nu)

matrix(
  Z,
  nrow = dim2,
  ncol = dim1,
  byrow = FALSE
) |> 
  as.data.frame() |> 
  as_tibble() |> 
  mutate(
    row = row_number()
  ) |> 
  pivot_longer(c(-row), names_to = "column", names_transform = parse_number) |> 
  ggplot(aes(column, row, fill = value)) +
  geom_raster() +
  scale_y_reverse()



fold_data(Z, dim1, dim2) |> 
  matrix(
    nrow = 2 * dim2,
    ncol = 2 * dim1
  ) |> 
  as.data.frame() |> 
  as_tibble() |> 
  mutate(
    row = row_number()
  ) |> 
  pivot_longer(c(-row), names_to = "column", names_transform = parse_number) |> 
  mutate(
    value = na_if(value, 0)
  ) |> 
  ggplot(aes(column, row, fill = value)) +
  geom_raster() +
  scale_y_reverse()


matrix(
  c(
    11, 21, 31, 41,
    12, 22, 32, 42,
    13, 23, 33, 43
  ),
  byrow = FALSE,
  nrow = 4,
  ncol = 3
) |> 
  as.numeric() |> 
  fold_data(dim1 = 3, dim2 = 4) |> 
  matrix(
    nrow = 8,
    ncol = 6
  ) |> 
  as.data.frame() |> 
  as_tibble() |> 
  mutate(
    row = row_number()
  ) |> 
  pivot_longer(c(-row), names_to = "column", names_transform = parse_number) |> 
  mutate(
    value = na_if(value, 0)
  ) |> 
  ggplot(aes(column, row, fill = value)) +
  geom_raster() +
  geom_text(aes(label = value, col = value)) +
  scale_colour_distiller(palette = "Greys", direction = 1) +
  scale_fill_distiller(palette = "Greys") +
  scale_x_continuous(position = "top") +
  scale_y_reverse()
