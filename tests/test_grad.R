library(stdmatern)
dim <- 50
rho <- 0.9
Q <- make_AR_prec_matrix(dim, rho)

L <- chol(Q)

z <- rnorm(dim)
x <- solve(L, z)


plot(x, type = "l")

ll <- function(x, rho) {
  Q <- make_AR_prec_matrix(dim, rho)
  logdet <- log(det(Q))
  quadform <- t(x) %*% Q %*% x
  as.numeric((logdet + quadform) / 2)
}

grad <- function(y, rho) {
  n <- length(y)
  y2 <- y^2
  y_lag <- y[-n] * y[-1]
  dlogdet <- (n - 1) * rho / (1 - rho^2)
  dquadform <- 2 * rho * (y2[1] + y2[n] + 2 * sum(y2[-c(1, n)]) - sum(y_lag))
  dquadform <- dquadform - 2 * sum(y_lag)
  dquadform <- dquadform / (1 - rho^2)^2
  dlogdet - dquadform / 2
}

rho_start <- 0.1
stepsize <- 1.2e-4
n_iter <- 400
values <- numeric(n_iter)
grads <- numeric(n_iter)
loglik <- numeric(n_iter)
values[1] <- rho_start

for (i in seq_len(n_iter - 1)) {
  grads[i] <- grad(x, values[i])
  loglik[i] <- ll(x, values[i])
  values[i + 1] <- values[i] + stepsize * mean(grads[seq(pmax(1, i - 4), i)])
}

loglik[n_iter] <- ll(x, values[n_iter])
grads[n_iter] <- grad(x, values[n_iter])


values
plot(values, type = "l")

grads
plot(grads, type = "l")

loglik
plot(loglik, type = "l")

plot(values, loglik, type = "l")
