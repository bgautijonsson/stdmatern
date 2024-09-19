library(stdmatern)

dim1 <- 3
dim2 <- 3
rho1 <- 0.5
rho2 <- 0.5
nu <- 0

Q1 <- make_AR_prec_matrix(3, 0.5)
Q2 <- make_AR_prec_matrix(3, 0.5)

I1 <- Diagonal(dim1, 1)
I2 <- Diagonal(dim2, 1)

Q <- kronecker(Q1, I2) + kronecker(I1, Q2)

