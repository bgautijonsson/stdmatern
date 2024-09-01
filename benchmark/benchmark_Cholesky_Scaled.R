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
    "Cholesky" = dmatern_copula_cholesky(X, dim, dim, rho, rho, nu),
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

results <- list()
sizes <- seq(10, 50, by = 10)

for (i in seq_along(sizes)) {
  results[[i]] <- my_fun(sizes[i])
}

results[[6]] <- my_fun(60)
results[[7]] <- my_fun(70)
results[[8]] <- my_fun(80)
results[[9]] <- my_fun(90)
results[[10]] <- my_fun(100)


tab <- results |> 
  list_rbind() |> 
  select(Grid = dim, Type = expression, time = median, memory = mem_alloc) |> 
  mutate(
    Grid = glue("{Grid}x{Grid}")
  ) |> 
  select(-memory) |> 
  pivot_wider(names_from = Type, values_from = time) |> 
  mutate(
    sp_1 = as.numeric(Cholesky/Eigen) |> round(2) |> paste0("x"),
    sp_2 = as.numeric(Cholesky/Circulant) |> round(2) |> paste0("x"),
    sp_3 = as.numeric(Cholesky / Folded) |> round(2) |> paste0("x")
  ) |> 
  select(
    Grid, 
    Cholesky,
    eig = Eigen, 
    sp_1,
    circ = Circulant, 
    sp_2, 
    fol = Folded, 
    sp_3
  ) |> 
  gt() |> 
  cols_label(
    Grid = "Grid Size",
    eig = "Time",
    sp_1 = "Speed-up",
    circ = "Time",
    sp_2 = "Speed-Up",
    fol = "Time",
    sp_3 = "Speed-Up"
  ) |> 
  tab_spanner(
    label = "Eigen",
    columns = 3:4
  ) |> 
  tab_spanner(
    label = "Circulant",
    columns = 5:6
  ) |> 
  tab_spanner(
    label = "Folded",
    columns = 7:8
  ) |> 
  tab_caption(
    md("Benchmarking how long it takes to evaluate the density of a Mátern($\\nu$)-like field with correlation parameter $\\rho$ scaled to have unit marginal variance")
  )

tab

tab |> 
  gtsave("benchmark_scaled.html")



results |> 
  list_rbind() |> 
  select(Grid = dim, Type = expression, time = median, memory = mem_alloc) |> 
  mutate(
    Grid = glue("{Grid}x{Grid}")
  ) |> 
  select(-memory) |> 
  pivot_wider(names_from = Type, values_from = time) |> 
  mutate(
    sp_1 = as.numeric(Cholesky/Eigen) |> round(2) |> paste0("x"),
    sp_2 = as.numeric(Cholesky/Circulant) |> round(2) |> paste0("x"),
    sp_3 = as.numeric(Cholesky / Folded) |> round(2) |> paste0("x")
  ) |> 
  select(
    Grid, 
    Cholesky,
    eig = Eigen, 
    sp_1,
    circ = Circulant, 
    sp_2, 
    fol = Folded, 
    sp_3
  ) |> 
  write_csv("benchmark_scaled.csv")

