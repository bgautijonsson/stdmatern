library(stdmatern)
library(tidyverse)
library(ggh4x)
theme_set(bggjphd::theme_bggj())

dim <- 20
rho <- 0.8
nu <- 2
n_obs <- 4000


exact <- cor(t(rmatern_copula(n_obs, dim, rho, nu)))[1, ]
circulant <- cor(t(rmatern_copula_circulant(n_obs, dim, rho, nu)))[1, ]
folded <- cor(t(rmatern_copula_folded_full(n_obs, dim, rho, nu)))[1, ]

tibble(
  Exact = exact,
  Circulant = circulant,
  Folded = folded
) |> 
  mutate(
    index = row_number()
  ) |> 
  pivot_longer(c(-index)) |> 
  ggplot(aes(index, value, col = name)) +
  geom_line() +
  scale_x_continuous(
    guide = guide_axis_truncated()
  ) +
  scale_y_continuous(
    breaks = seq(0, 1, by = 0.1),
    guide = guide_axis_truncated()
  ) +
  scale_colour_brewer(
    palette = "Set1"
  ) +
  labs(
    title = "Comparing the first line in the correlation matrices for each method",
    subtitle = "Shown for rho = 0.8, nu = 2 and dim = 20",
    x = NULL,
    y = "Correlation",
    col = NULL
  ) +
  theme(
    legend.position = "top"
  )





  dim <- 40
  rho <- 0.8
  nu <- 1

rmatern_copula(1, dim, rho, nu) |> 
  pnorm() |> 
  evd::qgev(loc = 10, scale = 5, shape = 0.1) |> 
  as.numeric() |> 
  tibble() |> 
  rename(
    y = 1
  ) |> 
  mutate(
    id = row_number(),
    lat = (id - 1) %% dim,
    lon = cumsum(lat == 0)
  ) |> 
  ggplot(aes(lat, lon)) +
  geom_raster(aes(fill = y))






