my_fun <- function(n, dim) {
  Z <- rmatern_copula_eigen(n, dim, dim, rho1, rho2, nu)
  out <- dmatern_copula_eigen(Z, dim, dim, rho1, rho2, nu)
  tibble(dim = dim, value = out, name = names(out))
}

library(purrr)

d <- crossing(
  n = c(1),
  dim = c(10, 20, 30, 40, 60, 80, 100, 120, 140, 160, 180, 200)
) |>
  mutate(
    iter = row_number()
  ) |>
  group_by(iter) |>
  reframe(
    n = n,
    dim = dim,
    my_fun(n, dim)
  )

d |>
  mutate(
    value = value / 1e9,
    value = value - lag(value, default = 0),
    perc = value / sum(value),
    .by = iter
  ) |> 
  filter(dim == max(dim))


d |>
  mutate(
    value = value / 1e9,
    value = value - lag(value, default = 0),
    value = value / sum(value),
    .by = iter
  ) |>
  ggplot(aes(x = dim, y = value, color = name)) +
  geom_smooth() +
  geom_point() +
  scale_y_continuous(
    # trans = "log10",
    labels = scales::label_percent()
  )


d |>
  mutate(
    value = value / 1e9,
    value = value - lag(value, default = 0),
    .by = iter
  ) |>
  ggplot(aes(x = dim, y = value, color = name)) +
  geom_smooth() +
  geom_point() +
  scale_y_continuous(
    trans = "log10",
    labels = scales::label_timespan()
  )
