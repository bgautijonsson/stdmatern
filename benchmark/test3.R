library(stdmatern)
library(sparseMVN)

grid_dim <- 10
rho <- 0.5
nu <- 0



Q1 <- make_AR_prec_matrix(grid_dim, rho)
eig1 <- eigen(Q1)
A1 <- eig1$values
V1 <- eig1$vectors

Q0 <- make_matern_prec_matrix(grid_dim, rho, 0)
Q <- Q0

for (i in seq_len(nu)) Q <- Q %*% Q0

mds1 <- fast_marginal_standard_deviations(A1, V1, grid_dim, nu)

mds1 <- numeric(grid_dim^2)

for (i in seq_len(grid_dim)) {
  for (j in seq_len(grid_dim)) {
    v <- Matrix::kronecker(V1[, j], V1[, i])
    a <- A1[i] + A1[j]
    a <- a^(nu + 1)
    
    mds1 <- mds1 + v^2 / a
  }
}
mds1 <- sqrt(mds1)

mds2 <- Q |> solve() |> diag() |> sqrt()

mean((mds1 - mds2))

plot(mds1, mds2)
abline(a = 0, b = 1)
