#' @export
matern_mvn_density <- function(X, dim, rho, nu) {
  if (is.vector(X) | (is.atomic(X) & NCOL(X) == 1)) {
    X <- matrix(X, ncol = 1)
  }

  matern_mvn_density_eigen_whitened(X, dim, rho, nu)
}

