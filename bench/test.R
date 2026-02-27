library(stdmatern)
library(Matrix)
library(purrr)
library(tidyverse)
library(glue)
library(gt)

rho <- 0.5
nu <- 0
n_x <- 40
n_y <- 40
n_rep <- 5


get_optim <- function(i) {
  X <- rmatern_copula(n = n_rep, dim = c(n_x, n_y), rho = c(rho,rho), nu = nu)
  optim(
    0.5,
    function(par) {
      -sum(dmatern_copula(X, dim = c(n_x, n_y), rho = c(par, par), nu = nu))
    }, 
    method = "Brent",
    lower = 0,
    upper = 1
  )$par
}



map_dbl(1:10, get_optim)