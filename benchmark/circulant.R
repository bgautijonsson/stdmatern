library(stdmatern)
library(tidyverse)
library(ggh4x)
theme_set(bggjphd::theme_bggj())

dim <- 30
rho <- 0.8
nu <- 2
n_obs <- 8000


exact <- cor(t(rmatern_copula_eigen(n_obs, dim, dim, rho, rho, nu)))[1, ]
circulant <- cor(t(rmatern_copula_circulant(n_obs, dim, dim, rho, rho, nu)))[1, ]
folded <- cor(t(rmatern_copula_folded_full(n_obs, dim, dim, rho, rho, nu)))[1, ]

tibble(
  Exact = exact,
  Circulant = circulant,
  Folded = folded
) |> 
  mutate(
    index = row_number()
  ) |> 
  pivot_longer(c(-index)) |> 
  mutate(
    name = fct_relevel(name, "Exact")
  ) |> 
  ggplot(aes(index, value, col = name)) +
  geom_smooth(
    se = 0,
    span = 0.021,
    n = 400
  ) +
  scale_x_continuous(
    guide = guide_axis_truncated(),
    breaks = breaks_extended(9)
  ) +
  scale_y_continuous(
    breaks = seq(0, 1, by = 0.1),
    guide = guide_axis_truncated()
  ) +
  scale_colour_manual(
    values = c(
      "#525252",
      "#e41a1c",
      "#377eb8"
    )
  ) +
  scale_fill_manual(
    values = c(
      "#525252",
      "#e41a1c",
      "#377eb8"
    )
  ) +
  labs(
    title = "Comparing the first line in the correlation matrices for each method",
    subtitle = "Shown for rho = 0.8 and nu = 2 on a 30x30 grid",
    x = "Column index",
    y = "Correlation",
    col = NULL
  ) +
  theme(
    legend.position = "top"
  )






