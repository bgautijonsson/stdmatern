library(tidyverse)
library(stdmatern)
dim1 <- 20
dim2 <- 20
rho1 <- 0.5
rho2 <- 0.3
nu <- 1
n <- 100

eigen <- rmatern_copula_eigen(n, dim1, dim2, rho1, rho2, nu) |> apply(1, sd)
folded <- rmatern_copula_folded_full(n, dim1, dim2, rho1, rho2, nu) |> apply(1, sd)
circulant <- rmatern_copula_circulant(n, dim1, dim2, rho1, rho2, nu) |> apply(1, sd)

tibble(
  eigen,
  folded,
  circulant
) |> 
  pivot_longer(c(everything())) |> 
  ggplot(aes(value, fill = name)) +
  geom_histogram(position = "identity", alpha = 0.6)


eigen <- rmatern_copula_eigen(n, dim1, dim2, rho1, rho2, nu) |> apply(1, mean)
folded <- rmatern_copula_folded_full(n, dim1, dim2, rho1, rho2, nu) |> apply(1, mean)
circulant <- rmatern_copula_circulant(n, dim1, dim2, rho1, rho2, nu) |> apply(1, mean)

tibble(
  eigen,
  folded,
  circulant
) |> 
  pivot_longer(c(everything())) |> 
  ggplot(aes(value, fill = name)) +
  geom_histogram(position = "identity", alpha = 0.6)