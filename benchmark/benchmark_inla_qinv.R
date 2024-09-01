library(stdmatern)
library(INLA)
library(purrr)
library(tidyverse)
library(broom)

dim <- 50
rho <- 0.5
nu <- 0



msd <- function(dim) {
  Q1 <- make_AR_prec_matrix(dim, rho)
  Q2 <- make_AR_prec_matrix(dim, rho)

  E1 <- eigen(Q1)
  E2 <- eigen(Q2)

  marginal_sd_eigen(
    E1$values,
    E1$vectors,
    dim,
    E2$values,
    E2$vectors,
    dim,
    nu
  ) |> sort()
}

inla_fun <- function(dim) {
  Q1 <- make_AR_prec_matrix(dim, rho)
  Q2 <- make_AR_prec_matrix(dim, rho)

  I1 <- Matrix::Diagonal(dim)
  I2 <- Matrix::Diagonal(dim)

  Q <- temp <- kronecker(Q1, I2) + kronecker(I1, Q2)
  for (i in seq_len(nu)) Q <- Q %*% temp
  inla.qinv(Q) |>
    diag() |>
    sqrt() |>
    sort()
}


my_fun <- function(dim) {
  i <<- i + 1
  print(i / length(my_seq))
  bench::mark(
    "inla" = inla_fun(dim),
    "stdmatern" = msd(dim),
    iterations = 1,
    filter_gc = FALSE
  )
}

i <- 0
my_seq <- seq(5, 205, by = 10)
d <- map(
  my_seq,
  my_fun
)


d |>
  list_rbind() |>
  select(expression, median) |>
  mutate(
    dim = my_seq,
    expression = as.character(expression),
    median = as.numeric(median),
    .by = expression
  ) |>
  reframe(
    lm(log(median) ~ log(dim)) |>
      tidy(),
    .by = expression
  )

d |>
  list_rbind() |>
  select(expression, median) |>
  mutate(
    dim = my_seq,
    expression = as.character(expression),
    median = as.numeric(median),
    .by = expression
  ) |>
  filter(
    dim == max(dim)
  )



d |>
  list_rbind() |>
  select(expression, median) |>
  mutate(
    dim = my_seq,
    expression = as.character(expression),
    median = as.numeric(median),
    .by = expression
  ) |>
  reframe(
    lm(log(median) ~ log(dim)) |>
      augment(),
    .by = expression
  ) |>
  janitor::clean_names() |>
  ggplot(aes(exp(log_dim), exp(log_median), col = expression)) +
  geom_point() +
  geom_line()
