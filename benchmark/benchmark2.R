library(stdmatern)
library(dplyr)
library(sparseMVN)
library(ggplot2)
grid_dim <- 30
rho <- 0.8
nu <- 2
Q <- make_standardized_matern(grid_dim, rho, nu)
L <- make_standardized_matern_cholesky(grid_dim, rho, nu)
chol_Q <- Matrix::Cholesky(Q, perm = FALSE)

comp_fun <- function(iter) {
  Z <- rmvn.sparse(
      n = 1,
      mu = rep(0, nrow(Q)),
      CH = chol_Q
  )

  R <- dmvn.sparse(
    Z,
    mu = rep(0, nrow(Q)),
    CH = chol_Q,
    log = TRUE
  )

  C <- dmvn_chol(
    as.numeric(Z), 
    L
  )

  tibble(
    R = R,
    C = C,
    iter = iter
  )
}

res <- purrr::map(1:100, comp_fun) |> 
  purrr::list_rbind()

res |> select(-iter) |> cor()

res |> 
  ggplot(aes(R, C)) +
  geom_abline(
    intercept = 0,
    slope = 1,
    lty = 2
  ) +
  geom_point()

#res |> 
#  ggplot(aes(C/R)) +
#  geom_density()