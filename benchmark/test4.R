library(stdmatern)
library(tidyverse)
library(purrr)
library(scales)
n <- 10
dim1 <- 30
dim2 <- 30
rho1 <- 0.8
rho2 <- 0.4
nu <- 0


test_fun <- function(iter) {
  X <- rmatern_copula_eigen(n, dim1, dim2, rho1, rho2, nu)


  folded <- optim(
    c(0.5, 0.5),
    function(par) {
      -sum(dmatern_copula_folded(X, dim1, dim2, par[1], par[2], nu))
    }
  )$par

  circulant <- optim(
    c(0.5, 0.5),
    function(par) {
      -sum(dmatern_copula_circulant(X, dim1, dim2, par[1], par[2], nu))
    }
  )$par

  eigen <- optim(
    c(0.5, 0.5),
    function(par) {
      -sum(dmatern_copula_eigen(X, dim1, dim2, par[1], par[2], nu))
    }
  )$par

  chol <- optim(
    c(0.5, 0.5),
    function(par) {
      -sum(dmatern_cholesky(X, dim1, dim2, par[1], par[2], nu))
    }
  )$par

  tibble(
    iter = iter,
    nu = nu,
    par = c("rho_1", "rho_2"),
    folded,
    eigen,
    circulant,
    chol = chol
  )
}


results <- map(1:10, test_fun) |> 
  list_rbind() 

results |> 
  rename(
    Folded = folded,
    Exact = eigen,
    Circulant = circulant,
    Cholesky = chol
  ) 
  # pivot_longer(c(-iter, -nu, -par)) |> 
  # pivot_wider(names_from = par) |> 
  # mutate(
  #   rho_1 = rho_1 / rho1,
  #   rho_2 = rho_2 / rho2,
  #   name = fct_relevel(name, "Exact", "Folded"),
  #   nu = paste0("nu = ", nu)
  # ) |> 
  # pivot_longer(c(rho_1, rho_2), names_to = "param") |> 
  # ggplot(aes(value)) +
  # geom_histogram(aes(fill = name), col = "white", bins = 50, position = "identity", alpha = 0.6) +
  # geom_vline(xintercept = 1, lty = 2) +
  # facet_wrap("param") +
  # scale_x_continuous(
  #   labels = function(x) percent(x - 1),
  #   trans = "log10",
  #   breaks = breaks_pretty(8)
  # ) +
  # scale_y_continuous(
  #   expand = expansion()
  # ) +
  # scale_fill_brewer(
  #   palette = "Set1"
  # ) +
  # # facet_wrap("nu", ncol = 1) +
  # labs(
  #   title = latex2exp::TeX("Bias and Variance in Maximum Likelihood estimates of $\\rho$"),
  #   x = "Difference between estimate and true value",
  #   y = NULL,
  #   fill = NULL
  # ) +
  # bggjphd::theme_bggj() +
  # theme(
  #   axis.line.y = element_blank(),
  #   axis.text.y = element_blank(),
  #   axis.ticks.y = element_blank(),
  #   legend.position = "top"
  # )
