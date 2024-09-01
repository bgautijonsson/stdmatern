library(stdmatern)
library(tidyverse)
library(purrr)
library(scales)
library(ggh4x)
theme_set(bggjphd::theme_bggj())


test_fun <- function(dim1, dim2, rho1, rho2, nu) {
  X <- rmatern_copula_eigen(n, dim1, dim2, rho1, rho2, nu)

  start <- Sys.time()
  folded <- optimize(
    function(par) {
      -sum(dmatern_copula_folded(X, dim1, dim2, par, par, nu))
    },
    interval = c(0, 1)
  )$minimum
  folded_time <- Sys.time() - start

  circulant <- optimize(
    function(par) {
      -sum(dmatern_copula_circulant(X, dim1, dim2, par, par, nu))
    },
    interval = c(0, 1)
  )$minimum
  circulant_time <- (Sys.time() - start) -  folded_time

  eigen <- optimize(
    function(par) {
      -sum(dmatern_copula_eigen(X, dim1, dim2, par, par, nu))
    },
    interval = c(0, 1)
  )$minimum
  eigen_time <- (Sys.time() - start) - circulant_time


  i <<- i + 1
  print(i / nrow(d))

  tibble(
    par = c("rho"),
    folded,
    eigen,
    circulant
  ) |>
    pivot_longer(c(-par), names_to = "model") |>
    pivot_wider(names_from = par) |>
    rename(rho_hat = rho) |> 
    arrange(model) |>
    mutate(
      time = as.numeric(c(circulant_time, eigen_time, folded_time))
    )
}

i <- 0

d <- crossing(
  dim = c(10, 20, 30),
  nu = c(0, 1, 2),
  rho = c(0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95),
  replicate = 1:5
) |>
  mutate(
    iter = row_number()
  )
d <- d |>
  group_by(iter) |>
  reframe(
    dim = dim,
    rho = rho,
    nu = nu,
    test_fun(dim, dim, rho, rho, nu)
  )

write_csv(d, "tests/ml_bench.csv")
d <- read_csv("tests/ml_bench.csv")
# Define colors
exact_color <- "#e41a1c"
circulant_color <- "#969696"
folded_color <- "#4daf4a"

d |>
  filter(dim %in% c(10, 20, 40, 80)) |> 
  mutate(
    grid_size = glue::glue("{dim} x {dim}"),
    model = fct_recode(
      model,
      "Exact" = "eigen",
      "Circulant" = "circulant",
      "Folded" = "folded"
    ) |> 
      fct_relevel("Exact")
  ) |>
  ggplot(aes(rho, rho_hat - rho)) +
  geom_abline(intercept = 0, slope = 0, lty = 2) +
  geom_smooth(aes(col = model, fill = model)) +
  scale_x_continuous(
    guide = guide_axis_truncated(),
    breaks = seq(0, 1, by = 0.2),
    limits = c(0, 1)
  ) +
  scale_y_continuous(
    guide = guide_axis_truncated()
  ) +
  scale_colour_manual(
    values = c(
      exact_color,
      circulant_color,
      folded_color
    )
  ) +
  scale_fill_manual(
    values = c(
      exact_color,
      circulant_color,
      folded_color
    )
  ) +
  facet_grid(
    cols = vars(grid_size),
    rows = vars(nu),
    labeller = label_both
  ) +
  coord_cartesian(ylim = c(-0.1, 0.1)) +
  labs(
    col = NULL,
    fill = NULL,
    x = expression(rho),
    y = "Error in ML estimate",
    title = "Benchmarking Performance of Maximum Likelihood Estimates",
    subtitle = "The folded approximation does better for larger grids and dependence"
  ) +
  theme(legend.position = "top")

ggsave(
  "bench_ml_bias.png",
  scale = 1.3, width = 8, height = 0.8 * 8
)


d |> 
  mutate(
    model = fct_recode(
      model,
      "Exact" = "eigen",
      "Circulant" = "circulant",
      "Folded" = "folded"
    ) |> 
      fct_relevel("Exact")
  ) |> 
  ggplot(aes(dim, time)) +
  geom_smooth(aes(col = model, fill = model)) +
    scale_x_continuous(
      guide = guide_axis_truncated(),
      breaks = breaks_extended(8),
      labels = \(x) glue::glue("{x}x{x}")
    ) +
    scale_y_continuous(
      guide = guide_axis_truncated(),
      labels = label_timespan(),
      breaks = breaks_extended(10)
    ) +
    scale_colour_manual(
      values = c(
        exact_color,
        circulant_color,
        folded_color
      )
    ) +
    scale_fill_manual(
      values = c(
        exact_color,
        circulant_color,
        folded_color
      )
    ) +
  labs(
    col = NULL,
    fill = NULL,
    x = "Grid Size",
    y = "Computation time",
    title = "Benchmarking Computation Time for Maximum Likelihood Estimation",
    subtitle = "The approximations scale better than the exact method"
  ) +
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.5, 0.9),
    legend.direction = "horizontal"
  )

  ggsave(
    "bench_ml_speed.png",
    scale = 1.3, width = 8, height = 0.621 * 8
  )
  

d |> 
  mutate(
    model = fct_recode(
      model,
      "Exact" = "eigen",
      "Circulant" = "circulant",
      "Folded" = "folded"
    ) |> 
      fct_relevel("Exact")
  ) |> 
  summarise(
    median = median(time),
    .by = c(model, dim)
  ) |> 
  pivot_wider(names_from = model, values_from = median) |> 
  select(dim, Exact, Folded, Circulant) |> 
  mutate(
    dim = glue::glue("{dim} x {dim}"),
    d_1 = Exact / Folded,
    d_2 = Exact / Circulant
  ) |> 
  mutate_at(
    vars(Exact, Folded, Circulant),
    bench::as_bench_time
  ) |> 
  mutate_at(
    vars(d_1, d_2),
    \(x) as.numeric(x) |> round(2) |> paste0("x"),
  ) |> 
  select(dim, Exact, Folded, d_1, Circulant, d_2) 
  write_csv("ml_bench_speed.csv")
