library(stdmatern)
library(tidyverse)
library(purrr)
library(scales)
n <- 10
dim <- 10
rho <- 0.5
nu <- 0

test_fun <- function(iter) {
  nu <- 0
  X <- rmatern_copula(n, dim, rho, nu)
  
  folded <- optimise(
    function(par) {
      dmatern_copula_folded(X, dim, par, nu) |> sum()
    },
    interval = c(0, 1), 
    maximum = TRUE
  )$maximum

  circulant <- optimise(
    function(par) {
      dmatern_copula_circulant(X, dim, par, nu) |> sum()
    },
    interval = c(0, 1), 
    maximum = TRUE
  )$maximum

  eigen <- optimise(
    function(par) {
      dmatern_copula_eigen(X, dim, dim, par, par, nu) |> sum()
    },
    interval = c(0, 1), 
    maximum = TRUE
  )$maximum

  tibble(
    iter = iter,
    nu = nu,
    folded,
    circulant, 
    eigen
  )
}


results <- map(1:100, test_fun) |> 
  list_rbind() 

results |> 
  rename(
    Circulant = circulant,
    Folded = folded,
    True = eigen
  ) |> 
  pivot_longer(c(-iter, -nu)) |> 
  mutate(
    diff = value / rho,
    name = fct_relevel(name, "Circulant", "True"),
    nu = paste0("nu = ", nu)
  ) |> 
  ggplot(aes(diff)) +
  geom_histogram(aes(fill = name), col = "white", bins = 50, position = "identity", alpha = 0.6) +
  geom_vline(xintercept = 1, lty = 2) +
  scale_x_continuous(
    labels = function(x) percent(x - 1),
    trans = "log10",
    breaks = breaks_pretty(8)
  ) +
  scale_y_continuous(
    expand = expansion()
  ) +
  scale_fill_brewer(
    palette = "Set1"
  ) +
  # facet_wrap("nu", ncol = 1) +
  labs(
    title = latex2exp::TeX("Bias and Variance in Maximum Likelihood estimates of $\\rho$"),
    x = "Difference between estimate and true value",
    y = NULL,
    fill = NULL
  ) +
  bggjphd::theme_bggj() +
  theme(
    axis.line.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "top"
  )
