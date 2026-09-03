#' stdmatern: standardised Matern Gaussian fields and copulas on rectangular grids
#'
#' Precision matrices, densities, samplers and a maximum-likelihood fit for the
#' standardised anisotropic Matern-like field
#' \eqn{Q = (Q_{\rho_1} \otimes I + I \otimes Q_{\rho_2})^{\nu + 1}} on an
#' \eqn{n_x \times n_y} grid, scaled so that its inverse is a correlation
#' matrix. Station order is x-major: cell \eqn{(a, b)} is row
#' \eqn{k = (a - 1) n_y + b}. All numerical work is single-threaded C++
#' (RcppEigen); the copula fit uses nlopt's BOBYQA through \pkg{nloptr}.
#'
#' @useDynLib stdmatern, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @importFrom nloptr nloptr
#' @importClassesFrom Matrix dgCMatrix
#' @keywords internal
"_PACKAGE"

.onLoad <- function(libname, pkgname) {
  # The fits reach nlopt through nloptr's registered C callables (nloptrAPI.h);
  # loading the namespace here guarantees they are registered before the first call.
  requireNamespace("nloptr", quietly = TRUE)
  invisible(NULL)
}
