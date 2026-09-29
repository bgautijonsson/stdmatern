# ---- argument checks (the R boundary; the C++ never throws) -------------------------------------

check_dim <- function(dim) {
  if (!is.numeric(dim) || !(length(dim) %in% 1:2) || anyNA(dim) ||
    any(dim != round(dim)) || any(dim <= 0)) {
    stop("dim must be positive integers, a scalar (square grid) or c(dim_x, dim_y)", call. = FALSE)
  }
  dim <- as.integer(dim)
  if (length(dim) == 1L) c(dim, dim) else dim
}

check_rho <- function(rho) {
  if (!is.numeric(rho) || !(length(rho) %in% 1:2) || anyNA(rho) || any(abs(rho) >= 1)) {
    stop("rho must lie in (-1, 1), a scalar (isotropic) or c(rho1, rho2)", call. = FALSE)
  }
  rho <- as.numeric(rho)
  if (length(rho) == 1L) c(rho, rho) else rho
}

check_nu <- function(nu) {
  if (!is.numeric(nu) || length(nu) != 1L || is.na(nu) || nu != round(nu) || nu < 0) {
    stop("nu must be a non-negative integer", call. = FALSE)
  }
  as.integer(nu)
}

check_n <- function(n) {
  if (!is.numeric(n) || length(n) != 1L || is.na(n) || n != round(n) || n < 0) {
    stop("n must be a non-negative integer", call. = FALSE)
  }
  as.integer(n)
}

check_seed <- function(seed) {
  if (is.null(seed)) {
    return(NULL)
  }
  if (!is.numeric(seed) || length(seed) != 1L || is.na(seed) || seed != round(seed) ||
    abs(seed) > .Machine$integer.max) {
    stop("seed must be a single whole number within the integer range, or NULL", call. = FALSE)
  }
  as.integer(seed)
}

check_obs <- function(X, dims, what = "X") {
  if (is.vector(X) || (is.atomic(X) && NCOL(X) == 1L)) {
    X <- matrix(X, ncol = 1L)
  }
  if (!is.numeric(X) || !is.matrix(X)) {
    stop(what, " must be a numeric matrix (one observation per column)", call. = FALSE)
  }
  N <- as.integer(dims[1]) * as.integer(dims[2])
  if (nrow(X) != N) {
    stop(sprintf("%s has %d rows but the %d x %d grid has %d cells", what, nrow(X), dims[1], dims[2], N),
      call. = FALSE
    )
  }
  if (!all(is.finite(X))) stop(what, " must be finite (no NA, NaN or Inf)", call. = FALSE)
  storage.mode(X) <- "double"
  X
}

check_start <- function(start, k) {
  if (is.null(start)) {
    return(NULL)
  }
  if (!is.numeric(start) || anyNA(start) || !(length(start) %in% 1:2)) {
    stop("start must be numeric, a scalar or c(rho1, rho2)", call. = FALSE)
  }
  if (any(start < 0.001 | start > 0.999)) {
    stop("start must lie within the search bounds [0.001, 0.999]", call. = FALSE)
  }
  start <- as.numeric(start)
  if (k == 1L) start[1] else if (length(start) == 1L) c(start, start) else start
}

check_count <- function(x, what) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || x != round(x) || x < 1) {
    stop(what, " must be a positive integer", call. = FALSE)
  }
  as.integer(x)
}

check_tol <- function(x, what) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x <= 0) {
    stop(what, " must be a positive number", call. = FALSE)
  }
  as.numeric(x)
}

#' @keywords internal
parse_dim_rho <- function(dim, rho) {
  d <- check_dim(dim)
  r <- check_rho(rho)
  list(dim_x = d[1], dim_y = d[2], rho1 = r[1], rho2 = r[2])
}

# ---- the standardised precision --------------------------------------------------------------------

#' Standardised Matern precision matrix
#'
#' The sparse precision \eqn{\tilde Q = D Q^{\nu+1} D} with
#' \eqn{Q = Q_{\rho_1} \otimes I + I \otimes Q_{\rho_2}} and
#' \eqn{D} the diagonal of marginal standard deviations, so that
#' \eqn{\tilde Q^{-1}} is a correlation matrix. Rows are in x-major station order
#' \eqn{k = (a - 1) n_y + b}.
#'
#' @param dim Grid size: a scalar (square) or `c(dim_x, dim_y)`, positive integers.
#' @param rho Correlation parameter(s) in (-1, 1): a scalar (isotropic) or `c(rho1, rho2)`.
#' @param nu Smoothness, a non-negative integer (the precision is raised to `nu + 1`).
#' @param method How the marginal standard deviations are computed: `"kronecker"`
#'   (default; \eqn{S^2 = W_1 M W_2'}, O(N (n_x + n_y))) or `"eigenpair"` (the
#'   2024 arithmetic, one Kronecker eigenvector at a time, O(N^2)), kept as the
#'   reference form. Equal to rounding (1e-13 relative).
#' @return A `dgCMatrix` of order `dim_x * dim_y`.
#' @seealso [matern_spectral()] for the same matrix together with its half
#'   log-determinant and eigen systems.
#' @examples
#' Q <- make_standardized_matern(c(6, 5), rho = 0.7, nu = 1)
#' range(diag(solve(Q))) # unit marginal variances
#' @export
make_standardized_matern <- function(dim, rho, nu, method = c("kronecker", "eigenpair")) {
  method <- match.arg(method)
  p <- parse_dim_rho(dim, rho)
  make_standardized_matern_eigen_cpp(p$dim_x, p$dim_y, p$rho1, p$rho2, check_nu(nu), method == "kronecker")
}

#' @rdname make_standardized_matern
#' @param dim_x,dim_y Grid dimensions (positive integers).
#' @param rho1,rho2 Correlation parameters along x and y, in (-1, 1).
#' @export
make_standardized_matern_eigen <- function(dim_x, dim_y, rho1, rho2, nu) {
  make_standardized_matern(c(dim_x, dim_y), c(rho1, rho2), nu)
}

#' Standardised precision with its spectral pieces
#'
#' One call returning what the copula likelihood and its consumers need:
#' the precision of [make_standardized_matern()] (bit for bit the same object),
#' its CSC slots, the half log-determinant
#' \eqn{\tfrac12 \log|\tilde Q| = \tfrac12\big((\nu+1)\sum_{ij}\log(\lambda_i + \mu_j) + 2\sum_k \log s_k\big)},
#' the marginal standard deviations \eqn{s_k} and the two AR(1) eigen systems.
#'
#' @inheritParams make_standardized_matern
#' @return A list: `Q` (`dgCMatrix`), `i`, `p`, `x` (its CSC slots), `hld`, `sd`,
#'   `e1` and `e2` (each `list(values, vectors)`), `dims`, `rho`, `nu`.
#' @examples
#' sp <- matern_spectral(c(6, 5), 0.7, 1)
#' sp$hld
#' @export
matern_spectral <- function(dim, rho, nu, method = c("kronecker", "eigenpair")) {
  method <- match.arg(method)
  p <- parse_dim_rho(dim, rho)
  nu <- check_nu(nu)
  cpp <- matern_spectral_cpp(p$dim_x, p$dim_y, p$rho1, p$rho2, nu, method == "kronecker")
  Q <- cpp$Q
  list(
    Q = Q, i = Q@i, p = Q@p, x = Q@x,
    hld = cpp$hld, sd = as.numeric(cpp$sd),
    e1 = list(values = as.numeric(cpp$A1), vectors = cpp$V1),
    e2 = list(values = as.numeric(cpp$A2), vectors = cpp$V2),
    dims = c(p$dim_x, p$dim_y), rho = c(p$rho1, p$rho2), nu = nu
  )
}

#' AR(1) precision matrix with unit marginal variance
#'
#' @param dim Length of the chain (a positive integer).
#' @param rho Correlation parameter in (-1, 1).
#' @return A sparse `dim x dim` tridiagonal precision matrix.
#' @export
make_AR_prec_matrix <- function(dim, rho) {
  if (!is.numeric(dim) || length(dim) != 1L || is.na(dim) || dim != round(dim) || dim <= 0) {
    stop("dim must be a positive integer", call. = FALSE)
  }
  if (!is.numeric(rho) || length(rho) != 1L || is.na(rho) || abs(rho) >= 1) {
    stop("rho must lie in (-1, 1)", call. = FALSE)
  }
  make_AR_prec_matrix_cpp(as.integer(dim), as.numeric(rho))
}

# ---- densities -------------------------------------------------------------------------------------

#' Gaussian copula log-density of the standardised Matern field
#'
#' Per-observation log-density of the Gaussian copula with precision
#' \eqn{\tilde Q} ([make_standardized_matern()]) evaluated at Gaussian scores:
#' \eqn{-\tfrac12 (z'\tilde Q z - \log|\tilde Q| - z'z)}.
#'
#' @param X Numeric matrix of Gaussian scores, one observation per column and
#'   `dim_x * dim_y` rows in x-major station order (a vector is one observation).
#' @inheritParams make_standardized_matern
#' @param method `"kronecker"` (default) evaluates the quadratic form through
#'   \eqn{T = V_2' Z V_1}, O(N (n_x + n_y) n); `"eigenpair"` evaluates one
#'   Kronecker eigenvector per eigenpair against every observation, O(N^2 n), the
#'   2024 reference form. The two agree to roughly 1e-10 per observation
#'   (summation order).
#' @return A numeric vector with one log-density per observation.
#' @examples
#' Z <- rmatern_copula(3, c(6, 5), 0.7, 1, seed = 1)
#' dmatern_copula(Z, c(6, 5), 0.7, 1)
#' @export
dmatern_copula <- function(X, dim, rho, nu, method = c("kronecker", "eigenpair")) {
  method <- match.arg(method)
  p <- parse_dim_rho(dim, rho)
  nu <- check_nu(nu)
  X <- check_obs(X, c(p$dim_x, p$dim_y))
  if (method == "kronecker") {
    as.numeric(dmatern_copula_kron_cpp(X, p$dim_x, p$dim_y, p$rho1, p$rho2, nu))
  } else {
    as.numeric(dmatern_copula_eigen_cpp(X, p$dim_x, p$dim_y, p$rho1, p$rho2, nu))
  }
}

#' @rdname dmatern_copula
#' @param dim_x,dim_y Grid dimensions (positive integers).
#' @param rho1,rho2 Correlation parameters along x and y, in (-1, 1).
#' @export
dmatern_copula_eigen <- function(X, dim_x, dim_y, rho1, rho2, nu) {
  dmatern_copula(X, c(dim_x, dim_y), c(rho1, rho2), nu)
}

#' Log-density of the unstandardised Matern-like field
#'
#' The multivariate normal log-density with precision
#' \eqn{Q = (Q_{\rho_1} \otimes I + I \otimes Q_{\rho_2})^{\nu+1}} (no
#' standardisation), per observation.
#'
#' @inheritParams dmatern_copula
#' @return A numeric vector with one log-density per observation.
#' @export
dmatern <- function(X, dim, rho, nu) {
  p <- parse_dim_rho(dim, rho)
  nu <- check_nu(nu)
  X <- check_obs(X, c(p$dim_x, p$dim_y))
  as.numeric(dmatern_eigen_cpp(X, p$dim_x, p$dim_y, p$rho1, p$rho2, nu))
}

# ---- samplers --------------------------------------------------------------------------------------

#' Sample the standardised Matern copula field
#'
#' Draws Gaussian scores with correlation \eqn{\tilde Q^{-1}}
#' ([make_standardized_matern()]). Each eigenpair's stream is seeded from
#' `seed`, so a draw is bit-for-bit reproducible from the seed alone; with
#' `seed = NULL` the base seed is drawn from R's RNG, so `set.seed()` governs it.
#'
#' Since 0.5.0 each eigenpair's generator is seeded from the full 64-bit hash of
#' (seed, pair index) (`seeding = "full"`). Up to 0.4.0 the hash was cut to 32
#' bits, so two eigenpairs of one draw could share a random stream (about 0.7% of
#' seeds on a 7,776-cell grid, 9% on 28,063 cells); `seeding = "legacy32"`
#' reproduces those draws exactly, for results produced before the fix. The
#' default is `getOption("stdmatern.seeding")`, else the environment variable
#' `STDMATERN_SEEDING`, else `"full"`, so a study can pin its scheme once.
#'
#' @param n Number of draws (a non-negative integer).
#' @inheritParams make_standardized_matern
#' @param seed A whole number, or `NULL` to draw the base seed from R's RNG.
#' @param method `"kronecker"` (default) forms each draw as
#'   \eqn{V_2 (\Lambda^{-(\nu+1)/2} \circ E) V_1'}, O(N (n_x + n_y) n);
#'   `"eigenpair"` adds one Kronecker eigenvector at a time, O(N^2 n), the 2024
#'   reference form. Both use the same per-eigenpair random numbers, so for a
#'   given seed they agree to rounding.
#' @param seeding `"full"` (the 0.5.0 scheme) or `"legacy32"` (the 0.4.0 draws);
#'   see Details.
#' @return A `dim_x * dim_y` by `n` matrix, rows in x-major station order.
#' @examples
#' Z <- rmatern_copula(4, c(6, 5), 0.7, 1, seed = 42)
#' identical(Z, rmatern_copula(4, c(6, 5), 0.7, 1, seed = 42))
#' @export
rmatern_copula <- function(n, dim, rho, nu, seed = NULL, method = c("kronecker", "eigenpair"),
                           seeding = default_seeding()) {
  method <- match.arg(method)
  legacy <- check_seeding(seeding)
  p <- parse_dim_rho(dim, rho)
  draw <- if (method == "kronecker") rmatern_copula_kron_cpp else rmatern_copula_eigen_cpp
  draw(check_n(n), p$dim_x, p$dim_y, p$rho1, p$rho2, check_nu(nu), check_seed(seed), legacy)
}

# The seeding scheme's default: the option, else the environment variable, else "full".
default_seeding <- function() {
  getOption("stdmatern.seeding", Sys.getenv("STDMATERN_SEEDING", "full"))
}

# TRUE for the 0.4.0 scheme, FALSE for the full one; anything else is refused.
check_seeding <- function(seeding) {
  if (!is.character(seeding) || length(seeding) != 1L || !seeding %in% c("full", "legacy32")) {
    stop('seeding must be "full" or "legacy32"', call. = FALSE)
  }
  identical(seeding, "legacy32")
}

#' @rdname rmatern_copula
#' @param dim_x,dim_y Grid dimensions (positive integers).
#' @param rho1,rho2 Correlation parameters along x and y, in (-1, 1).
#' @export
rmatern_copula_eigen <- function(n, dim_x, dim_y, rho1, rho2, nu, seed = NULL,
                                 seeding = default_seeding()) {
  rmatern_copula(n, c(dim_x, dim_y), c(rho1, rho2), nu, seed, seeding = seeding)
}

#' Sample the unstandardised Matern-like field
#'
#' @inheritParams rmatern_copula
#' @return A `dim_x * dim_y` by `n` matrix, rows in x-major station order.
#' @export
rmatern <- function(n, dim, rho, nu, seed = NULL, seeding = default_seeding()) {
  legacy <- check_seeding(seeding)
  p <- parse_dim_rho(dim, rho)
  rmatern_eigen_cpp(check_n(n), p$dim_x, p$dim_y, p$rho1, p$rho2, check_nu(nu), check_seed(seed), legacy)
}

# ---- the fit ---------------------------------------------------------------------------------------

#' Fit the copula correlation parameter(s) by maximum likelihood
#'
#' BOBYQA (nlopt, through \pkg{nloptr}) on the bounded box `[0.001, 0.999]`
#' per parameter, maximising the summed [dmatern_copula()] log-density. There
#' is no wall-clock cap: the deterministic budget is `maxeval`, and the result
#' says how it stopped.
#'
#' @param Y Numeric matrix of Gaussian scores, one observation per column and
#'   `dim_x * dim_y` rows in x-major station order.
#' @inheritParams make_standardized_matern
#' @param start Starting value(s): a scalar, or `c(rho1, rho2)` for the
#'   anisotropic fit; must lie within the search bounds. Default 0.5.
#' @param isotropic If `TRUE`, one `rho` for both directions.
#' @param method Objective form, see [dmatern_copula()].
#' @param maxeval Maximum number of objective evaluations (positive integer).
#' @param xtol_rel,ftol_rel Relative tolerances on the parameter and on the objective.
#' @return A list: `rho` (scalar, or `c(rho1, rho2)` with `rho1`/`rho2` also
#'   given), `loglik`, `n_evals` (integer), `converged` (`TRUE` for nlopt codes
#'   1--4), `stop_reason` (`"converged"`, `"maxeval"`, `"roundoff"` or
#'   `"failed"`), `convergence` (the raw nlopt code: 1 success, 2 stopval,
#'   3 ftol, 4 xtol, 5 maxeval, 6 maxtime, negative failure), `message`, `method`.
#' @examples
#' Y <- rmatern_copula(60, c(8, 8), 0.6, 1, seed = 3)
#' fit_matern_copula(Y, c(8, 8), 1, isotropic = TRUE)[c("rho", "n_evals", "stop_reason")]
#' @export
fit_matern_copula <- function(Y, dim, nu, start = NULL, isotropic = FALSE,
                              method = c("kronecker", "eigenpair"),
                              maxeval = 1000L, xtol_rel = 1e-6, ftol_rel = 1e-8) {
  method <- match.arg(method)
  d <- check_dim(dim)
  nu <- check_nu(nu)
  Y <- check_obs(Y, d, what = "Y")
  maxeval <- check_count(maxeval, "maxeval")
  xtol_rel <- check_tol(xtol_rel, "xtol_rel")
  ftol_rel <- check_tol(ftol_rel, "ftol_rel")
  kron <- method == "kronecker"
  if (isotropic) {
    start <- check_start(start, 1L)
    result <- fit_matern_copula_eigen_isotropic_cpp(Y, d[1], d[2], nu, start, kron, maxeval, xtol_rel, ftol_rel)
  } else {
    start <- check_start(start, 2L)
    result <- fit_matern_copula_eigen_cpp(Y, d[1], d[2], nu, start, kron, maxeval, xtol_rel, ftol_rel)
    result$rho <- c(result$rho1, result$rho2)
  }
  result
}
