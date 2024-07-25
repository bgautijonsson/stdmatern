#' @export
dmatern_copula <- function(X, dim, rho, nu) {
  if (is.vector(X) | (is.atomic(X) & NCOL(X) == 1)) {
    X <- matrix(X, ncol = 1)
  }

  dmatern_copula_eigen(X, dim, rho, nu)
}

#' @export
dmatern <- function(X, dim, rho, nu) {
  if (is.vector(X) | (is.atomic(X) & NCOL(X) == 1)) {
    X <- matrix(X, ncol = 1)
  }

  dmatern_eigen(X, dim, rho, nu)
}

#' @export
make_standardized_matern <- function(dim, rho, nu) {
  make_standardized_matern_eigen(dim, rho, nu)
}