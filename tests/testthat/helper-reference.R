# Dense reference for the standardised Matern copula on a small grid -- an independent transcription
# of the maths that shares no code with src/ (vendored from
# maxandsmooth/tests/reference/reference_density.R, 2026-09-02). Conventions: AR(1) precision with
# unit marginal variance (off-diagonal -rho/(1-rho^2), interior diagonal (1+rho^2)/(1-rho^2),
# corner 1/(1-rho^2)); Q = (Q1 (x) I_n2 + I_n1 (x) Q2)^(nu+1); station order k = (a-1)*n2 + b.
ref_ar1_prec <- function(n, rho) {
  s <- 1 / (1 - rho^2)
  d <- rep((1 + rho^2) * s, n)
  d[1] <- s
  d[n] <- s
  Q <- diag(d, n)
  for (i in seq_len(n - 1)) {
    Q[i, i + 1] <- -rho * s
    Q[i + 1, i] <- -rho * s
  }
  Q
}
ref_Q_power <- function(n1, n2, rho1, rho2, nu) {
  Qs <- kronecker(ref_ar1_prec(n1, rho1), diag(n2)) + kronecker(diag(n1), ref_ar1_prec(n2, rho2))
  Q <- Qs
  if (nu > 0) for (i in seq_len(nu)) Q <- Q %*% Qs
  Q
}
ref_marginal_sd <- function(n1, n2, rho1, rho2, nu) {
  sqrt(diag(solve(ref_Q_power(n1, n2, rho1, rho2, nu))))
}
ref_Qtilde <- function(n1, n2, rho1, rho2, nu) {
  Q <- ref_Q_power(n1, n2, rho1, rho2, nu)
  sd <- sqrt(diag(solve(Q)))
  diag(sd) %*% Q %*% diag(sd)
}
ref_half_logdet <- function(Qt) 0.5 * as.numeric(determinant(Qt, logarithm = TRUE)$modulus)
ref_copula_logdens <- function(X, Qt) {
  -0.5 * (colSums(X * (Qt %*% X)) - 2 * ref_half_logdet(Qt) - colSums(X^2))
}
ref_normal_logdens <- function(X, Q) {
  N <- nrow(Q)
  -0.5 * (N * log(2 * pi) - as.numeric(determinant(Q, logarithm = TRUE)$modulus) + colSums(X * (Q %*% X)))
}
