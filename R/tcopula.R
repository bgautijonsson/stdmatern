# The Student-t copula of the standardised Matern field (0.6.0). R-level only: the correlation
# matrix is the inverse of the same Qtilde as the Gaussian copula (matern_spectral()), the Gaussian
# part of a draw is rmatern_copula(), and the mixing variable is drawn from R's own RNG.

# ---- argument checks -------------------------------------------------------------------------------

check_df <- function(df, what = "df") {
  if (!is.numeric(df) || length(df) != 1L || !is.finite(df) || df <= 0) {
    stop(what, " must be a single finite number > 0", call. = FALSE)
  }
  as.numeric(df)
}

check_df_bounds <- function(df_bounds) {
  if (!is.numeric(df_bounds) || length(df_bounds) != 2L || !all(is.finite(df_bounds)) ||
    df_bounds[1] <= 0 || df_bounds[1] >= df_bounds[2]) {
    stop("df_bounds must be c(lower, upper) with 0 < lower < upper, both finite", call. = FALSE)
  }
  as.numeric(df_bounds)
}

# ---- internal pieces -------------------------------------------------------------------------------

# Per-observation t-copula log density from t scores X (P x n), the sparse Qtilde, its half
# log-determinant and the summed univariate t log densities of each column (ldt).
tcopula_logdens_core <- function(X, Q, hld, df, ldt) {
  P <- nrow(X)
  q <- colSums(X * as.matrix(Q %*% X))
  K <- lgamma((df + P) / 2) - lgamma(df / 2) - (P / 2) * log(df * pi)
  K + hld - (df + P) / 2 * log1p(q / df) - ldt
}

# Messages and stop reasons as the C++ fits report them (src/matern_copula_eigen.cpp, run_bobyqa).
nlopt_message <- function(code) {
  switch(as.character(code),
    "1" = "Optimization succeeded",
    "2" = "Stopval reached",
    "3" = "Function tolerance reached",
    "4" = "Parameter tolerance reached",
    "5" = "Maximum evaluations reached",
    "6" = "Maximum time reached",
    "-4" = "Roundoff limited",
    "-5" = "Forced stop",
    "-2" = "Optimization failed: invalid arguments",
    "-3" = "Optimization failed: out of memory",
    "Optimization failed"
  )
}

nlopt_stop_reason <- function(code) {
  if (code >= 1L && code <= 4L) {
    "converged"
  } else if (code == 5L) {
    "maxeval"
  } else if (code == 6L) {
    "maxtime"
  } else if (code == -4L) {
    "roundoff"
  } else {
    "failed"
  }
}

# ---- exported functions ----------------------------------------------------------------------------

#' Map Gaussian scores to Student-t scores
#'
#' Elementwise \eqn{x = F_\nu^{-1}(\Phi(z))} for a t distribution with `df`
#' degrees of freedom, computed on the log scale from the nearer tail,
#' \eqn{x = -\mathrm{sign}(z)\, F_\nu^{-1}(\log \Phi(-|z|))}, so the map is exactly
#' odd and keeps full precision in both tails (no \eqn{1 - \Phi(z)} cancellation).
#'
#' @param Z Numeric vector or matrix of Gaussian scores (`NA` propagates).
#' @param df Degrees of freedom, a single finite number > 0.
#' @return `Z` with every entry mapped (dimensions and attributes kept).
#' @examples
#' gaussian_to_t(c(-8, -1, 0, 1, 8), df = 5)
#' @export
gaussian_to_t <- function(Z, df) {
  if (!is.numeric(Z)) stop("Z must be numeric", call. = FALSE)
  df <- check_df(df)
  out <- Z
  storage.mode(out) <- "double"
  out[] <- -sign(Z) * stats::qt(stats::pnorm(-abs(Z), log.p = TRUE), df, log.p = TRUE)
  out
}

#' Student-t copula log-density of the standardised Matern field
#'
#' Per-observation log-density of the t copula with correlation matrix
#' \eqn{R = \tilde Q^{-1}} ([make_standardized_matern()]) and `df` degrees of
#' freedom, evaluated at t scores \eqn{x} (see [gaussian_to_t()]):
#' \deqn{\log c(x) = K + \tfrac12 \log|\tilde Q| - \tfrac{\nu_t + P}{2}
#' \log\big(1 + q / \nu_t\big) - \sum_{k=1}^P \log t_{\nu_t}(x_k),}
#' with \eqn{q = x'\tilde Q x},
#' \eqn{K = \log\Gamma((\nu_t + P)/2) - \log\Gamma(\nu_t/2) - \tfrac P2 \log(\nu_t \pi)}
#' and \eqn{t_{\nu_t}} the univariate t density. \eqn{\tilde Q} and its half
#' log-determinant come from [matern_spectral()]; the quadratic forms are sparse
#' products. As `df` grows the density tends to [dmatern_copula()] at the
#' matching Gaussian scores.
#'
#' @param X Numeric matrix of t scores, one observation per column and
#'   `dim_x * dim_y` rows in x-major station order (a vector is one observation).
#' @inheritParams make_standardized_matern
#' @param df Degrees of freedom of the t copula, a single finite number > 0.
#' @param method How [matern_spectral()] computes the marginal standard
#'   deviations of \eqn{\tilde Q}: `"kronecker"` (default) or `"eigenpair"`.
#' @return A numeric vector with one log-density per observation.
#' @examples
#' X <- rtmatern_copula(3, c(6, 5), 0.7, 1, df = 5, seed = 1)
#' dtmatern_copula(X, c(6, 5), 0.7, 1, df = 5)
#' @export
dtmatern_copula <- function(X, dim, rho, nu, df, method = c("kronecker", "eigenpair")) {
  method <- match.arg(method)
  p <- parse_dim_rho(dim, rho)
  nu <- check_nu(nu)
  df <- check_df(df)
  X <- check_obs(X, c(p$dim_x, p$dim_y))
  sp <- matern_spectral(c(p$dim_x, p$dim_y), c(p$rho1, p$rho2), nu, method)
  ldt <- colSums(stats::dt(X, df, log = TRUE))
  tcopula_logdens_core(X, sp$Q, sp$hld, df, ldt)
}

#' Sample the Student-t copula of the standardised Matern field
#'
#' Draws t scores \eqn{x_i = z_i / \sqrt{W_i}} with \eqn{z_i} a draw of
#' [rmatern_copula()] (correlation \eqn{\tilde Q^{-1}}) and
#' \eqn{W_i \sim \mathrm{Gamma}(\nu_t/2, \text{rate} = \nu_t/2)}, so each column
#' is multivariate t with `df` degrees of freedom and unit-scale t margins.
#'
#' With a `seed`, \eqn{z} is exactly `rmatern_copula(n, dim, rho, nu, seed, method, seeding)`
#' and \eqn{W} is drawn by `rgamma()` after `set.seed(seed)` under
#' `RNGkind("Mersenne-Twister", "Inversion", "Rejection")`; the caller's
#' `.Random.seed` (or its absence) and `RNGkind()` are restored on exit, so a
#' seeded draw is reproducible and leaves the session stream untouched. With
#' `seed = NULL` both the base seed of \eqn{z} and then \eqn{W} come from the
#' session stream, so `set.seed()` governs the draw.
#'
#' @inheritParams rmatern_copula
#' @param df Degrees of freedom, a single finite number > 0.
#' @return A `dim_x * dim_y` by `n` matrix of t scores, rows in x-major station
#'   order, with attribute `"W"` holding the `n` mixing variables.
#' @examples
#' X <- rtmatern_copula(4, c(6, 5), 0.7, 1, df = 5, seed = 42)
#' attr(X, "W")
#' @export
rtmatern_copula <- function(n, dim, rho, nu, df, seed = NULL, method = c("kronecker", "eigenpair"),
                            seeding = default_seeding()) {
  method <- match.arg(method)
  n <- check_n(n)
  df <- check_df(df)
  seed <- check_seed(seed)
  if (!is.null(seed)) {
    env <- globalenv()
    had_seed <- exists(".Random.seed", envir = env, inherits = FALSE)
    old_seed <- if (had_seed) get(".Random.seed", envir = env, inherits = FALSE)
    old_kind <- RNGkind()
    on.exit(
      {
        suppressWarnings(RNGkind(old_kind[1], old_kind[2], old_kind[3]))
        if (had_seed) {
          assign(".Random.seed", old_seed, envir = env)
        } else if (exists(".Random.seed", envir = env, inherits = FALSE)) {
          rm(".Random.seed", envir = env)
        }
      },
      add = TRUE
    )
  }
  Z <- rmatern_copula(n, dim, rho, nu, seed = seed, method = method, seeding = seeding)
  if (!is.null(seed)) {
    set.seed(seed, kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
  }
  W <- stats::rgamma(n, shape = df / 2, rate = df / 2)
  X <- Z / rep(sqrt(W), each = nrow(Z))
  attr(X, "W") <- W
  X
}

#' Fit the t-copula correlation parameter(s) and degrees of freedom
#'
#' Maximum likelihood for the t copula of [dtmatern_copula()] by BOBYQA (nlopt,
#' through [nloptr::nloptr()]) over `rho` in `[0.001, 0.999]` (one value, or
#' two when `isotropic = FALSE`) and \eqn{\log \nu_t} in `log(df_bounds)`, or
#' over `rho` alone when `df` is given. The data are Gaussian scores, as for
#' [fit_matern_copula()]; each trial `df` maps them to t scores by
#' [gaussian_to_t()] (cached for the last `df` value, as is the spectral
#' decomposition for the last `rho`). A non-finite objective is replaced by
#' 1e15. As for [fit_matern_copula()] there is no wall-clock cap.
#'
#' @param Y Numeric matrix of Gaussian scores, one observation per column and
#'   `dim_x * dim_y` rows in x-major station order.
#' @inheritParams make_standardized_matern
#' @param df `NULL` (default) to estimate the degrees of freedom, or a single
#'   finite number > 0 to hold them fixed.
#' @param start Starting value(s) for `rho`: a scalar, or `c(rho1, rho2)` for the
#'   anisotropic fit, within `[0.001, 0.999]`. Default 0.5.
#' @param df_start Starting degrees of freedom, inside `df_bounds`.
#' @param df_bounds Search box for the degrees of freedom, `c(lower, upper)` with
#'   `0 < lower < upper`, both finite (searched on the log scale).
#' @param isotropic If `TRUE`, one `rho` for both directions.
#' @param method Objective form, see [dtmatern_copula()].
#' @param maxeval Maximum number of objective evaluations (positive integer).
#' @param xtol_rel,ftol_rel Relative tolerances on the parameters and on the objective.
#' @return A list: `rho` (scalar, or `c(rho1, rho2)` with `rho1`/`rho2` also
#'   given), `df`, `loglik`, `convergence` (the raw nlopt code), `message`,
#'   `n_evals` (integer), `converged` (`TRUE` for nlopt codes 1--4),
#'   `stop_reason` (`"converged"`, `"maxeval"`, `"maxtime"`, `"roundoff"` or
#'   `"failed"`), `method`, `df_fixed` (`TRUE` when `df` was given) and
#'   `df_at_bound` (`TRUE` when the estimate lies within 1e-6 relative of a
#'   bound of `df_bounds`; always `FALSE` for a fixed `df`).
#' @examples
#' X <- rtmatern_copula(40, c(8, 8), 0.6, 1, df = 5, seed = 3)
#' Y <- qnorm(pt(X, 5))
#' fit_tmatern_copula(Y, c(8, 8), 1, isotropic = TRUE)[c("rho", "df", "stop_reason")]
#' @export
fit_tmatern_copula <- function(Y, dim, nu, df = NULL, start = NULL, df_start = 10,
                               df_bounds = c(1, 500), isotropic = FALSE,
                               method = c("kronecker", "eigenpair"),
                               maxeval = 1000L, xtol_rel = 1e-6, ftol_rel = 1e-8) {
  method <- match.arg(method)
  d <- check_dim(dim)
  nu <- check_nu(nu)
  Y <- check_obs(Y, d, what = "Y")
  df_fixed <- !is.null(df)
  if (df_fixed) df <- check_df(df)
  df_bounds <- check_df_bounds(df_bounds)
  df_start <- check_df(df_start, "df_start")
  if (df_start < df_bounds[1] || df_start > df_bounds[2]) {
    stop("df_start must lie within df_bounds", call. = FALSE)
  }
  maxeval <- check_count(maxeval, "maxeval")
  xtol_rel <- check_tol(xtol_rel, "xtol_rel")
  ftol_rel <- check_tol(ftol_rel, "ftol_rel")
  k <- if (isotropic) 1L else 2L
  start <- check_start(start, k)
  if (is.null(start)) start <- rep(0.5, k)

  # single-entry caches: BOBYQA often moves one coordinate at a time
  t_df <- NA_real_
  t_X <- NULL
  t_ldt <- NULL
  sp_rho <- c(NA_real_, NA_real_)
  sp <- NULL
  loglik_at <- function(rho, dfv) {
    if (!identical(dfv, t_df)) {
      t_X <<- gaussian_to_t(Y, dfv)
      t_ldt <<- colSums(stats::dt(t_X, dfv, log = TRUE))
      t_df <<- dfv
    }
    if (!identical(rho, sp_rho)) {
      sp <<- matern_spectral(d, rho, nu, method)
      sp_rho <<- rho
    }
    sum(tcopula_logdens_core(t_X, sp$Q, sp$hld, dfv, t_ldt))
  }
  # WHY: nloptr::nloptr() calls eval_f twice at x0 before handing it to nlopt; the memo keeps those
  # probes from costing two extra evaluations, and n_evals is nlopt's own count.
  last_x <- NULL
  last_val <- NULL
  objective <- function(x) {
    if (identical(x, last_x)) {
      return(last_val)
    }
    rho <- if (k == 1L) c(x[1], x[1]) else x[1:2]
    dfv <- if (df_fixed) df else exp(x[k + 1L])
    val <- -loglik_at(rho, dfv)
    if (!is.finite(val)) val <- 1e15
    last_x <<- x
    last_val <<- val
    val
  }

  x0 <- start
  lb <- rep(0.001, k)
  ub <- rep(0.999, k)
  if (!df_fixed) {
    x0 <- c(x0, log(df_start))
    lb <- c(lb, log(df_bounds[1]))
    ub <- c(ub, log(df_bounds[2]))
  }
  res <- nloptr::nloptr(
    x0 = x0, eval_f = objective, lb = lb, ub = ub,
    opts = list(algorithm = "NLOPT_LN_BOBYQA", xtol_rel = xtol_rel, ftol_rel = ftol_rel, maxeval = maxeval)
  )
  x <- res$solution
  code <- as.integer(res$status)
  df_hat <- if (df_fixed) df else exp(x[k + 1L])
  df_at_bound <- !df_fixed &&
    (abs(df_hat - df_bounds[1]) <= 1e-6 * df_bounds[1] || abs(df_hat - df_bounds[2]) <= 1e-6 * df_bounds[2])
  rho_part <- if (k == 1L) list(rho = x[1]) else list(rho = x[1:2], rho1 = x[1], rho2 = x[2])
  c(rho_part, list(
    df = df_hat,
    loglik = if (is.finite(res$objective) && res$objective < 1e15) -res$objective else NA_real_,
    convergence = code,
    message = nlopt_message(code),
    n_evals = as.integer(res$iterations),
    converged = code >= 1L && code <= 4L,
    stop_reason = nlopt_stop_reason(code),
    method = method,
    df_fixed = df_fixed,
    df_at_bound = df_at_bound
  ))
}
