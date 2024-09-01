library(stdmatern)
library(Matrix)
library(purrr)
library(tidyverse)
library(glue)
library(gt)
library(ggh4x)
library(scales)
library(geomtextpath)
theme_set(bggjphd::theme_bggj())

rho <- 0.5


my_fun <- function(dim, nu) {
  X <- rmatern_eigen(1, dim, dim, rho, rho, nu)
  bench::mark(
    chol = dmatern_cholesky(X, dim, dim, rho, rho, nu),
    eig = dmatern_eigen(X, dim, dim, rho, rho, nu),
    filter_gc = FALSE,
    iterations = 20,
    check = FALSE
  ) 
}

results <- crossing(
  dim = seq(10, 200, by = 10),
  nu = 0:2
) |> 
  mutate(
    results = map2(dim, nu, my_fun)
  )

results |> 
  unnest(results) |> 
  mutate(
    expression = as.character(expression),
    nu = factor(nu),
    median = as.numeric(median)
  ) |> 
  select(dim, nu, expression, median) |>  
  pivot_wider(names_from = expression, values_from = median) |> 
  mutate(
    diff = eig / chol
  ) |> 
  ggplot(aes(dim, diff, col = nu)) +
  geom_hline(yintercept = 1, lty = 2) +
  geom_line() +
  scale_x_continuous(
    breaks = c(seq(10, 160, by = 30), 200),
    labels = function(x) glue("{x^2} x {x^2}"),
    guide = guide_axis_truncated()
  ) +
  scale_y_continuous(
    breaks = c(1, 0.25, 0.5, 0.75, 1.25, 1.5, 2, 4),
    labels = function(x) {
      x <- x - 1
      case_when(
        x == 0 ~ "Equal",
        x < 0 ~ glue("{percent(abs(x))} faster"),
        x > 0 ~ glue("{percent(abs(x))} slower")
      )
    },
    trans = "log10",
    guide = guide_axis_truncated()
  ) +
  scale_colour_brewer(
    palette = "Set1",
    guide = guide_legend(
      override.aes = list(linewidth = 2)
    )
  ) +
  labs(
    title = "How fast is the eigen version compared to the cholesky version?",
    x = "Size of Q",
    y = NULL,
    col = expression(nu)
  ) +
  theme(
    legend.position = "top",
    plot.margin = margin(5, 25, 5, 5)
  )

  results |> 
    unnest(results) |> 
    mutate(
      expression = as.character(expression) |> 
        fct_recode(
          "Cholesky" = "chol",
          "Eigen" = "eig"
        ),
      nu = factor(nu),
      median = as.numeric(median)
    ) |> 
    select(dim, nu, expression, median) |> 
    ggplot(aes(dim, median)) +
    geom_line(aes(lty = expression, group = paste(expression, nu))) +
    geom_textline(
      data = ~filter(.x, nu == 0),
      aes(label = expression, group = expression, hjust = expression, lty = expression),
      vjust = -0.2,
      text_only = TRUE,
      size = 6
    ) +
    scale_x_continuous(
      breaks = c(10, 50, 100, 150, 200),
      labels = function(x) glue("{x^2} x {x^2}"),
      guide = guide_axis_truncated()
    ) +
    scale_y_continuous(
      labels = label_timespan(),
      breaks = breaks_extended(10),
      guide = guide_axis_truncated()
    ) +
    scale_linetype_discrete(
      guide = "none"
    ) +
    scale_hjust_manual(
      values = c(0.8, 0.8)
    ) +
    labs(
      title = "Benchmarking density evaluations",
      x = "Size of Q",
      y = NULL,
      col = expression(nu),
      lty = NULL
    ) +
    theme(
      legend.position.inside = c(0.5, 0.8),
      plot.margin = margin(5, 25, 5, 5)
    )

ggsave(
  "unscaled_benchmark.png",
  scale = 0.9, width = 8, height = 8
)

results |> 
  unnest(results) |> 
  mutate(
    expression = as.character(expression)
  ) |> 
  filter(nu == 0) |> 
  select(dim, expression, median) |> 
  pivot_wider(names_from = expression, values_from = median) |> 
  mutate(
    diff = as.numeric(chol / eig) |> round(2) |> paste0("x"),
    dim = glue("{dim^2}x{dim^2}")
  ) |> 
  write_csv("unscaled_benchmark.csv")

results |> 
  unnest(results) |> 
  mutate(
    expression = as.character(expression)
  ) |> 
  filter(nu == 0) |> 
  select(dim, expression, median) |> 
  pivot_wider(names_from = expression, values_from = median) |> 
  mutate(
    diff = as.numeric(chol / eig) |> round(2) |> paste0("x"),
    dim = glue("{dim^2}x{dim^2}")
  ) |> 
  gt() |> 
  cols_label(
    dim = "Q-size",
    chol = "Cholesky",
    eig = "Eigen",
    diff = "Cholesky / Eigen"
  ) |> 
  tab_caption(
    md("Benchmarking how long it takes to evaluate the density of a Mátern($\\nu$)-like field with correlation parameter $\\rho$, unscaled version")
  )
