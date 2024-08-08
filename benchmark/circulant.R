dim <- 5
rho <- 0.5

base_Q1 <- numeric(dim)
base_Q1[1] <- 1 + rho^2
base_Q1[c(2, dim)] <- -rho

Q1 <- matrix(
  0,
  nrow = dim,
  ncol = dim
)

for (i in 1:dim) {
  Q1[i, (i + 1:dim - 2) %% dim + 1] <- base_Q1
}

I <- diag(1, nrow = dim, ncol = dim)

Q <- kronecker(Q1, I) + kronecker(I, Q1)



base_Q <- matrix(
  0,
  nrow = dim,
  ncol = dim
)

base_Q[1, ] <- base_Q1 + c(1 + rho^2, rep(0, dim - 1)) 
base_Q[c(2, dim), 1] <- -rho

fft(base_Q)

mean(abs(sort(Re(fft(base_Q)), decreasing = TRUE) - eigen(Q)$values))

Q_inv <- solve(Q)

temp <- base_Q

for (i in seq_len(dim)) {
  temp[i, ] <- Q_inv[(i - 1) * dim + 1, seq_len(dim)]
}

#temp
Re(fft(1 / fft(base_Q), inverse = TRUE) / dim^2)

mean(abs(temp - Re(fft(1 / fft(base_Q), inverse = TRUE) / dim^2)))

base_C <- base_Q

eigenvalues <- fft(base_C)

base_C_inverse <- fft(1 / eigenvalues, inverse = TRUE) / dim^2

mvar <- base_C_inverse[1, 1]

updated_eigenvalues <- mvar * eigenvalues


