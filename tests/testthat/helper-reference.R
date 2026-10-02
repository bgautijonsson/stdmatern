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

# Dense reference for the Student-t copula (0.6.0): correlation R = solve(Qt), multivariate t density
# with df degrees of freedom at t scores X (one observation per column) minus the univariate t log
# densities, each written out from the gamma-function formula (no dt(), no sparse algebra).
ref_t_copula_logdens <- function(X, Qt, df) {
  R <- solve(Qt)
  P <- nrow(R)
  q <- colSums(X * solve(R, X))
  logdet_R <- as.numeric(determinant(R, logarithm = TRUE)$modulus)
  joint <- lgamma((df + P) / 2) - lgamma(df / 2) - (P / 2) * log(df * pi) - 0.5 * logdet_R -
    (df + P) / 2 * log(1 + q / df)
  margin <- lgamma((df + 1) / 2) - lgamma(df / 2) - 0.5 * log(df * pi) - (df + 1) / 2 * log(1 + X^2 / df)
  joint - colSums(margin)
}
# Inverse of gaussian_to_t written independently: Gaussian scores from t scores through the nearer tail.
ref_t_to_gaussian <- function(X, df) {
  -sign(X) * qnorm(pt(-abs(X), df, log.p = TRUE), log.p = TRUE)
}
