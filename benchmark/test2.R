library(stdmatern)
library(sparseMVN)

grid_dim <- 3
rho <- 0.1
n_replicate <- 10
nu <- 0

Q <- make_standardized_matern_eigen(grid_dim, rho, nu = nu)

Z <- rmvn.sparse(
    n = n_replicate,
    mu = rep(0, nrow(Q)),
    CH = Matrix::Cholesky(Q)
) |> t()


neg_ll_function <- function(params) {
    rho <- params[1]
    - sum(matern_mvn_density(Z, grid_dim, rho, nu))
}

res <- optimize(
    f = neg_ll_function,
    interval = c(0, 1)
)

res



