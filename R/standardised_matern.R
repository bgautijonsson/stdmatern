#' @export
matern_mvn_density <- function(X, dim, rho, nu) {
  if (is.vector(X) | (is.atomic(X) & NCOL(X) == 1)) {
    X <- matrix(X, ncol = 1)
  }

  matern_mvn_density_eigen(X, dim, rho, nu)
}


#' @export
dmvn_chol <- function(X, L) {
  if (is.vector(X) | (is.atomic(X) & NCOL(X) == 1)) {
    X <- matrix(X, ncol = 1)
  }

  dmvn_chol_eigen(X, L)
}