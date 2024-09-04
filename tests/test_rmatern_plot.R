library(stdmatern)
library(tidyverse)

dim1 <- 80
dim2 <- 100
rho1 <- 0.9
rho2 <- 0.8
nu <- 2
n <- 1

eigen <- rmatern_copula_eigen(n, dim1, dim2, rho1, rho2, nu)
folded <- rmatern_copula_folded_full(n, dim1, dim2, rho1, rho2, nu)
circulant <- rmatern_copula_circulant(n, dim1, dim2, rho1, rho2, nu)

tibble(
  eigen = as.numeric(eigen),
  folded = as.numeric(folded),
  circulant = as.numeric(circulant)
) |>
  pivot_longer(c(everything())) |>
  mutate(
    id = row_number(),
    lat = rep(seq_len(dim1), each = dim2),
    lon = rep(seq_len(dim2), times = dim1),
    value = (value - mean(value)) / sd(value),
    .by = name
  ) |>
  ggplot(aes(lat, lon, fill = value)) +
  geom_raster() +
  scale_fill_viridis_c() +
  facet_wrap("name") +
  coord_fixed(expand = FALSE)
