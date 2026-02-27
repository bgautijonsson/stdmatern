#' @keywords internal
parse_dim_rho <- function(dim, rho) {
  if (length(dim) == 1) {
    dim_x <- dim
    dim_y <- dim
  } else {
    dim_x <- dim[1]
    dim_y <- dim[2]
  }
  if (length(rho) == 1) {
    rho1 <- rho
    rho2 <- rho
  } else {
    rho1 <- rho[1]
    rho2 <- rho[2]
  }
  list(dim_x = dim_x, dim_y = dim_y, rho1 = rho1, rho2 = rho2)
}

#' @export
dmatern_copula <- function(X, dim, rho, nu) {
  if (is.vector(X) | (is.atomic(X) & NCOL(X) == 1)) {
    X <- matrix(X, ncol = 1)
  }
  p <- parse_dim_rho(dim, rho)
  dmatern_copula_eigen(X, p$dim_x, p$dim_y, p$rho1, p$rho2, nu)
}

#' @export
dmatern <- function(X, dim, rho, nu) {
  if (is.vector(X) | (is.atomic(X) & NCOL(X) == 1)) {
    X <- matrix(X, ncol = 1)
  }
  p <- parse_dim_rho(dim, rho)
  dmatern_eigen(X, p$dim_x, p$dim_y, p$rho1, p$rho2, nu)
}

#' @export
rmatern_copula <- function(n, dim, rho, nu) {
  p <- parse_dim_rho(dim, rho)
  rmatern_copula_eigen(n, p$dim_x, p$dim_y, p$rho1, p$rho2, nu)
}

#' @export
rmatern <- function(n, dim, rho, nu) {
  p <- parse_dim_rho(dim, rho)
  rmatern_eigen(n, p$dim_x, p$dim_y, p$rho1, p$rho2, nu)
}

#' @export
make_standardized_matern <- function(dim, rho, nu) {
  p <- parse_dim_rho(dim, rho)
  make_standardized_matern_eigen(p$dim_x, p$dim_y, p$rho1, p$rho2, nu)
}

#' Fit Matern copula parameters via MLE
#'
#' Estimates the correlation parameters of a standardized Matern copula by
#' maximum likelihood estimation using the BOBYQA derivative-free optimization
#' algorithm.
#'
#' @param Y A matrix of observations, where each column is a single observation
#'   of dimension dim_x * dim_y (or dim^2 if dim is scalar).
#' @param dim Either a single integer (for square grids) or a vector of two
#'   integers c(dim_x, dim_y) specifying the grid dimensions.
#' @param nu The smoothness parameter for the Matern field (non-negative integer).
#' @param start Optional starting values for optimization. For anisotropic case,
#'   either a single value (used for both rho1 and rho2) or a vector c(rho1, rho2).
#'   For isotropic case, a single value. Defaults to 0.5.
#' @param isotropic Logical. If TRUE, estimates a single rho parameter that
#'   applies to both spatial directions. If FALSE (default), estimates separate
#'   rho1 and rho2 parameters.
#'
#' @return A list containing:
#' \describe{
#'   \item{rho}{Estimated correlation parameter(s). A scalar if isotropic=TRUE,
#'     or a vector c(rho1, rho2) if isotropic=FALSE.}
#'   \item{rho1}{Estimated rho1 parameter (only if isotropic=FALSE)}
#'   \item{rho2}{Estimated rho2 parameter (only if isotropic=FALSE)}
#'   \item{loglik}{Log-likelihood at the MLE}
#'   \item{convergence}{Integer convergence code from NLopt (positive = success)}
#'   \item{message}{Convergence message from the optimizer}
#' }
#'
#' @examples
#' # Anisotropic case: estimate separate rho1 and rho2
#' Y <- rmatern_copula(n = 100, dim = 10, rho = c(0.7, 0.3), nu = 1)
#' fit <- fit_matern_copula(Y, dim = 10, nu = 1)
#' fit$rho
#'
#' # Isotropic case: estimate single rho for both directions
#' Y <- rmatern_copula(n = 100, dim = 10, rho = 0.6, nu = 1)
#' fit <- fit_matern_copula(Y, dim = 10, nu = 1, isotropic = TRUE)
#' fit$rho
#'
#' @export
fit_matern_copula <- function(Y, dim, nu, start = NULL, isotropic = FALSE) {
  # Handle input dimensions
  if (is.vector(Y) || (is.atomic(Y) && NCOL(Y) == 1)) {
    Y <- matrix(Y, ncol = 1)
  }

  if (length(dim) == 1) {
    dim_x <- dim
    dim_y <- dim
  } else {
    dim_x <- dim[1]
    dim_y <- dim[2]
  }

  if (isotropic) {
    # Isotropic case: single rho for both directions
    if (!is.null(start) && length(start) > 1) {
      start <- start[1] # Use first value if vector provided
    }

    result <- fit_matern_copula_eigen_isotropic(Y, dim_x, dim_y, nu, start)
  } else {
    # Anisotropic case: separate rho1 and rho2
    if (!is.null(start)) {
      if (length(start) == 1) {
        start <- c(start, start)
      }
    }

    result <- fit_matern_copula_eigen(Y, dim_x, dim_y, nu, start)

    # Add combined rho vector for convenience
    result$rho <- c(result$rho1, result$rho2)
  }

  result
}
