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
    "Eigen" = dmatern_copula_eigen(X, dim, dim, rho, rho, nu),
    "Circulant" = dmatern_copula_circulant(X, dim, rho, nu),
    "Folded" = dmatern_copula_folded(X, dim, rho, nu),
    filter_gc = FALSE,
    iterations = 20,
    check = FALSE
  ) |> 
    mutate(
      dim = dim
    )
}


results <- map(seq(10, 200, by = 10), my_fun)


results |> 
  list_rbind() |> 
  select(Q_size = dim, Type = expression, time = median, memory = mem_alloc) |> 
  mutate(
    Q_size = glue("{Q_size^2}x{Q_size^2}")
  ) |> 
  select(-memory) |> 
  pivot_wider(names_from = Type, values_from = time) |> 
  mutate(
    sp_1 = as.numeric(Eigen/Circulant) |> round(2) |> paste0("x"),
    sp_2 = as.numeric(Eigen/Folded) |> round(2) |> paste0("x")
  ) |> 
  select(
    Q_size, Exact = Eigen, circ = Circulant, sp_1, fol = Folded, sp_2
  ) |> 
  gt() |> 
  cols_label(
    circ = "Time",
    sp_1 = "Speed-Up",
    fol = "Time",
    sp_2 = "Speed-Up"
  ) |> 
  tab_spanner(
    label = "Circulant",
    columns = 3:4
  ) |> 
    tab_spanner(
      label = "Folded",
      columns = 5:6
    ) |> 
  tab_caption(
    md("Benchmarking how long it takes to evaluate the density of a Mátern($\\nu$)-like field with correlation parameter $\\rho$, scaled to have unit marginal variance")
  )


Tdim <- 20
rho <- 0.9
nu <- 2
X <- rmatern_copula_circulant(1, dim, rho, nu)

tibble(
  Z = as.numeric(X)
) |> 
  mutate(
    id = row_number(),
    lat = (id - 1) %% dim[1],
    lon = cumsum(lat == 0),
  ) |> 
  ggplot(aes(lat, lon, fill = Z)) +
  geom_raster() +
  scale_fill_viridis_c() +
  coord_fixed(expand = FALSE)

dmatern_eigen(X, dim, dim, rho, rho, nu)
dmatern_copula_circulant(X, dim, rho, nu)

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
  select(Q_size = dim, Type = expression, time = median, memory = mem_alloc) |> 
  mutate(
    time = as.numeric(time),
    Q_size = Q_size^2
  ) |> 
  group_by(Type) |> 
  reframe(
    lm(log(time) ~ log(Q_size)) |> 
      broom::tidy()
  )


  results |> 
    list_rbind() |> 
    select(dim, type = expression, time = median, memory = mem_alloc) |> 
    mutate(
      locations = dim^2,
      type = as.character(type),
      time = as.numeric(time)
    ) |> 
    ggplot(aes(locations, time)) +
    geom_line(aes(lty = type)) +
    scale_y_log10()
  
