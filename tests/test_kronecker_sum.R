library(stdmatern)
library(tidyverse)

dim1 <- 3
dim2 <- 4
rho1 <- 0.4
rho2 <- 0.5
nu <- 0

Q1 <- make_AR_prec_matrix(dim1, rho1)
I1 <- Diagonal(dim1, 1)
Q2 <- make_AR_prec_matrix(dim2, rho2)
I2 <- Diagonal(dim2, 1)

Q <- kronecker(Q1, I2) + kronecker(I1, Q2)
temp <- Q
for (i in seq_len(nu)) Q <- Q %*% temp

E1 <- eigen(Q1)
E2 <- eigen(Q2)


A <- numeric(dim1 * dim2)
V <- matrix(
  0,
  nrow = dim1 * dim2,
  ncol = dim1 * dim2
)

idx <- 1
for (i in seq_len(dim1)) {
  for (j in seq_len(dim2)) {
    A[idx] <- (E1$values[i] + E2$values[j])^(nu + 1)
    V[ , idx] <- kronecker(E1$vectors[, i], E2$vectors[, j])
    idx <- idx + 1
  }
}

zapsmall(Q - (V %*% diag(A) %*% t(V))) |> round(3)
zapsmall(make_matern_prec_matrix(dim1, dim2, rho1, rho2, nu) - (V %*% diag(A) %*% t(V))) |> round(3)
