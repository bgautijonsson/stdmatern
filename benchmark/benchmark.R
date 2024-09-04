library(stdmatern)
library(Matrix)
library(purrr)
library(tidyverse)
library(glue)
library(gt)

rho <- 0.5
nu <- 0

my_fun <- function(dim) {
  X <- rmatern_copula_eigen(1, dim, dim, rho, rho, nu)
  bench::mark(
    "Cholesky (Unscaled)" = dmatern_cholesky(X, dim, dim, rho, rho, nu),
    "Eigen (Unscaled)" = dmatern_eigen(X, dim, dim, rho, rho, nu),
    "Eigen" = dmatern_copula_eigen(X, dim, dim, rho, rho, nu),
    "Circulant" = dmatern_copula_circulant(X, dim, dim, rho, rho, nu),
    "Folded" = dmatern_copula_folded(X, dim, dim, rho, rho, nu),
    filter_gc = FALSE,
    iterations = 20,
    check = FALSE
  ) |>
    mutate(
      dim = dim
    )
}


results <- map(
  c(seq(20, 240, by = 20)),
  my_fun,
  .progress = TRUE
)

results |>
  list_rbind() |>
  select(Grid = dim, Type = expression, time = median, memory = mem_alloc) |>
  mutate(
    Grid = glue("{Grid}x{Grid}")
  ) |>
  select(-memory) |>
  pivot_wider(names_from = Type, values_from = time) |>
  mutate(
    sp_1 = as.numeric(Circulant / Eigen),
    sp_2 = as.numeric(Folded / Eigen),
    sp_3 = as.numeric(`Eigen (Unscaled)` / `Cholesky (Unscaled)`)
  ) |>
  mutate_at(
    vars(starts_with("sp_")),
    \(x) scales::percent(x, accuracy = 0.1, decimal_mark = ".", big.mark = ",")
  ) |>
  select(
    Grid,
    "Cholesky (Unscaled)",
    "Eigen (Unscaled)",
    sp_3,
    eig = Eigen,
    circ = Circulant,
    sp_1,
    fol = Folded,
    sp_2
  ) |>
  # write_csv("benchmark_all.csv")
  gt() |>
  cols_label(
    `Cholesky (Unscaled)` = "Cholesky",
    `Eigen (Unscaled)` = "Time",
    sp_3 = "Relative",
    eig = "Eigen",
    circ = "Time",
    sp_1 = "Relative",
    fol = "Time",
    sp_2 = "Relative"
  ) |>
  tab_spanner(
    label = "Circulant",
    columns = 6:7
  ) |>
  tab_spanner(
    label = "Folded",
    columns = 8:9
  ) |>
  tab_spanner(
    label = "Eigen",
    columns = 3:4
  ) |>
  tab_spanner(
    label = "Unscaled",
    2:4
  ) |>
  tab_spanner(
    label = "Scaled",
    columns = 5:9
  ) |>
  tab_caption(
    md("Benchmarking how long it takes to evaluate the density of a Mátern($\\nu$)-like field with correlation parameter $\\rho$, either unscaled or scaled to have unit marginal variance")
  )






results |>
  list_rbind() |>
  select(dim, expression, time = median) |>
  mutate(
    time = as.numeric(time),
    expression = as.character(expression),
    dim = dim^2
  ) |>
  reframe(
    lm(log(time) ~ log(dim)) |>
      broom::tidy(),
    .by = expression
  ) |>
  filter(term == "log(dim)")

results |>
  list_rbind() |>
  select(dim, expression, time = median) |>
  mutate(
    time = as.numeric(time),
    expression = as.character(expression),
    dim = dim^2
  ) |>
  reframe(
    lm(log(time) ~ log(dim)) |>
      broom::augment(),
    .by = expression
  ) |>
  janitor::clean_names() |>
  ggplot(aes(exp(log_dim), exp(log_time), col = expression)) +
  geom_point() +
  geom_line() +
  scale_y_log10(
    labels = scales::label_timespan()
  )
