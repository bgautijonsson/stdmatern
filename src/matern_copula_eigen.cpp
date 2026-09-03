#include <RcppEigen.h>
#include <random>
#include <string>
#include <vector>
#include <cmath>
#include "ar_matrix.h"
#include "rng_seed.h"
// nlopt through nloptr's registered C callables: no link against a system
// libnlopt, no Homebrew path (review F7, 2026-09-02). The nloptr namespace
// is loaded in .onLoad so the callables are registered before the first fit.
extern "C" {
#include <nloptrAPI.h>
}

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(nloptr)]]

using namespace Rcpp;
using namespace Eigen;

// ---------------------------------------------------------------------------
// Everything here runs single-threaded since 2026-09-02 (review F2): the
// OpenMP regions merged per-thread partial sums in arrival order, so results
// were bit-reproducible only at OMP_NUM_THREADS=1. The loops below keep the
// row-major (i, j) accumulation order of the one-thread execution, so a build
// of this file reproduces the pinned gate results bit for bit.
// Arguments are validated at the R boundary (R/standardised_matern.R); the
// C++ never throws -- a malformed call returns an empty or NA result.
// ---------------------------------------------------------------------------

// Marginal standard deviations of the field with precision
// (Q1 (x) I + I (x) Q2)^(nu + 1), from the two AR(1) eigen systems:
// sd_k^2 = sum_ij (V1_ai V2_bj)^2 / (A1_i + A2_j)^(nu + 1), k = a * dim_y + b.
// [[Rcpp::export]]
Eigen::ArrayXd marginal_sd_eigen(
    const Eigen::VectorXd& A1, const Eigen::MatrixXd& V1, int dim_x,
    const Eigen::VectorXd& A2, const Eigen::MatrixXd& V2, int dim_y,
    int nu
) {
    ArrayXd marginal_sds = ArrayXd::Zero(dim_x * dim_y);
    VectorXd v(dim_x * dim_y);
    for (int i = 0; i < dim_x; ++i) {
        for (int j = 0; j < dim_y; ++j) {
            v = kroneckerProduct(V1.col(i), V2.col(j));
            marginal_sds += v.array().square() / std::pow(A1(i) + A2(j), nu + 1);
        }
    }
    return marginal_sds.sqrt();
}

namespace {

// The standardised precision Q~ = D Q^(nu + 1) D and the pieces it is built
// from. This is the body of make_standardized_matern_eigen (2024) unchanged
// in its statements and their order; matern_spectral reuses it so both
// return the same Q~ bit for bit.
struct StandardizedPieces {
    VectorXd A1, A2;
    MatrixXd V1, V2;
    ArrayXd sd;
    SparseMatrix<double> Q_std;
};

StandardizedPieces standardized_pieces(int dim_x, int dim_y, double rho1, double rho2, int nu) {
    StandardizedPieces out;
    int N = dim_x * dim_y;
    SparseMatrix<double> Q1 = make_AR_prec_matrix_cpp(dim_x, rho1);
    SparseMatrix<double> Q2 = make_AR_prec_matrix_cpp(dim_y, rho2);

    SelfAdjointEigenSolver<SparseMatrix<double>> solver1(Q1);
    out.A1 = solver1.eigenvalues();
    out.V1 = solver1.eigenvectors();

    SelfAdjointEigenSolver<SparseMatrix<double>> solver2(Q2);
    out.A2 = solver2.eigenvalues();
    out.V2 = solver2.eigenvectors();

    out.sd = marginal_sd_eigen(out.A1, out.V1, dim_x, out.A2, out.V2, dim_y, nu);

    DiagonalMatrix<double, Dynamic> D(N);
    D.diagonal() = out.sd;

    SparseMatrix<double> I1(dim_x, dim_x);
    I1.setIdentity();
    SparseMatrix<double> I2(dim_y, dim_y);
    I2.setIdentity();

    SparseMatrix<double> Q = kroneckerProduct(Q1, I2) + kroneckerProduct(I1, Q2);

    if (nu > 0) {
        SparseMatrix<double> temp = Q;
        for (int i = 0; i < nu; ++i) {
            Q = Q * temp;
            Q.makeCompressed();
        }
    }

    out.Q_std = D * Q * D;
    return out;
}

}  // namespace

// Standardised Matern precision Q~ (sparse), whose inverse is a correlation matrix.
// [[Rcpp::export]]
Eigen::SparseMatrix<double> make_standardized_matern_eigen_cpp(int dim_x, int dim_y, double rho1, double rho2, int nu) {
    if (dim_x <= 0 || dim_y <= 0) return SparseMatrix<double>(0, 0);
    return standardized_pieces(dim_x, dim_y, rho1, rho2, nu).Q_std;
}

// Q~ plus its half log-determinant, marginal sds and the two eigen systems
// from one call (review F6). hld = 0.5 * ((nu + 1) sum_ij log(A1_i + A2_j) + 2 sum_k log sd_k).
// [[Rcpp::export]]
Rcpp::List matern_spectral_cpp(int dim_x, int dim_y, double rho1, double rho2, int nu) {
    if (dim_x <= 0 || dim_y <= 0) return Rcpp::List::create();
    StandardizedPieces sp = standardized_pieces(dim_x, dim_y, rho1, rho2, nu);
    double log_lam = 0.0;
    for (int i = 0; i < dim_x; ++i)
        for (int j = 0; j < dim_y; ++j)
            log_lam += std::log(sp.A1(i) + sp.A2(j));
    double hld = 0.5 * ((nu + 1.0) * log_lam + 2.0 * sp.sd.log().sum());
    return Rcpp::List::create(
        Rcpp::Named("Q") = sp.Q_std,
        Rcpp::Named("hld") = hld,
        Rcpp::Named("sd") = sp.sd,
        Rcpp::Named("A1") = sp.A1, Rcpp::Named("V1") = sp.V1,
        Rcpp::Named("A2") = sp.A2, Rcpp::Named("V2") = sp.V2
    );
}

// Gaussian copula log-density, the 2024 "eigenpair" form: one Kronecker
// eigenvector per pair, dotted against every observation -- O(N^2 n) per call.
// Kept unchanged (bar the OpenMP removal) because the production two-stage fit
// and its pins run on it; the O(N (n1 + n2) n) form is dmatern_copula_kron_cpp.
// [[Rcpp::export]]
Eigen::ArrayXd dmatern_copula_eigen_cpp(const Eigen::MatrixXd& X, int dim_x, int dim_y, double rho1, double rho2, int nu) {
    int n_obs = X.cols();
    int N = X.rows();
    if (dim_x <= 0 || dim_y <= 0 || N != dim_x * dim_y) return ArrayXd::Constant(n_obs, NA_REAL);

    MatrixXd Q1 = make_AR_prec_matrix_cpp(dim_x, rho1);
    MatrixXd Q2 = make_AR_prec_matrix_cpp(dim_y, rho2);

    SelfAdjointEigenSolver<MatrixXd> solver1(Q1);
    VectorXd A1 = solver1.eigenvalues();
    MatrixXd V1 = solver1.eigenvectors();

    SelfAdjointEigenSolver<MatrixXd> solver2(Q2);
    VectorXd A2 = solver2.eigenvalues();
    MatrixXd V2 = solver2.eigenvectors();

    ArrayXd marginal_sds = marginal_sd_eigen(A1, V1, dim_x, A2, V2, dim_y, nu);

    MatrixXd X_scaled = (X.array().colwise() * marginal_sds).matrix();

    // log|Q~| = log|Q| + 2 log|D|, D = diag(marginal_sds)
    double log_det_D = marginal_sds.log().sum();
    double log_det_Q = 0.0;

    ArrayXd quadform_sums = ArrayXd::Zero(n_obs);
    VectorXd v(dim_x * dim_y);
    for (int i = 0; i < dim_x; ++i) {
        for (int j = 0; j < dim_y; ++j) {
            v = kroneckerProduct(V1.col(i), V2.col(j));
            double lambda = std::pow(A1(i) + A2(j), nu + 1.0);
            log_det_Q += std::log(lambda);
            for (int k = 0; k < n_obs; ++k) {
                double dot_product = v.dot(X_scaled.col(k));
                quadform_sums(k) += dot_product * dot_product * lambda;
            }
        }
    }

    double log_det = log_det_Q + 2 * log_det_D;
    ArrayXd x_squared = X.colwise().squaredNorm().array();
    ArrayXd log_densities = -0.5 * (quadform_sums - log_det - x_squared);
    return log_densities;
}

// The same log-density through the Kronecker structure (review F1): for each
// observation reshape x (station order k = a * dim_y + b) into the dim_y x dim_x
// matrix Z(b, a), form T = V2' (Z o sd) V1, and x'Q~x = sum_ij lam_ij T_ji^2 with
// lam_ij = (A1_i + A2_j)^(nu + 1). Marginal variances S2 = W1 M W2', W = V o V,
// M = 1/lam. Cost O(N (n1 + n2)) per observation instead of O(N^2).
// [[Rcpp::export]]
Eigen::ArrayXd dmatern_copula_kron_cpp(const Eigen::MatrixXd& X, int dim_x, int dim_y, double rho1, double rho2, int nu) {
    const int n_obs = X.cols();
    const int N = X.rows();
    if (dim_x <= 0 || dim_y <= 0 || N != dim_x * dim_y) return ArrayXd::Constant(n_obs, NA_REAL);

    MatrixXd Q1 = make_AR_prec_matrix_cpp(dim_x, rho1);
    MatrixXd Q2 = make_AR_prec_matrix_cpp(dim_y, rho2);
    SelfAdjointEigenSolver<MatrixXd> solver1(Q1);
    SelfAdjointEigenSolver<MatrixXd> solver2(Q2);
    const VectorXd& A1 = solver1.eigenvalues();
    const MatrixXd& V1 = solver1.eigenvectors();
    const VectorXd& A2 = solver2.eigenvalues();
    const MatrixXd& V2 = solver2.eigenvectors();

    ArrayXXd lam(dim_x, dim_y);
    for (int i = 0; i < dim_x; ++i)
        for (int j = 0; j < dim_y; ++j)
            lam(i, j) = std::pow(A1(i) + A2(j), nu + 1.0);

    MatrixXd W1 = V1.array().square().matrix();
    MatrixXd W2 = V2.array().square().matrix();
    MatrixXd S2 = W1 * lam.inverse().matrix() * W2.transpose();   // (a, b)
    ArrayXd sd(N);
    for (int a = 0; a < dim_x; ++a)
        for (int b = 0; b < dim_y; ++b)
            sd(a * dim_y + b) = std::sqrt(S2(a, b));

    const double log_det = lam.log().sum() + 2.0 * sd.log().sum();

    MatrixXd X_scaled = (X.array().colwise() * sd).matrix();
    Map<const MatrixXd> Xm(X_scaled.data(), dim_y, dim_x * n_obs);   // column a + dim_x * s holds Z_s(., a)
    MatrixXd A = V2.transpose() * Xm;                                // (j, a + dim_x * s)
    ArrayXXd lamT = lam.transpose();                                 // (j, i)
    ArrayXd quad(n_obs);
    for (int s = 0; s < n_obs; ++s) {
        MatrixXd T = A.block(0, s * dim_x, dim_y, dim_x) * V1;      // (j, i)
        quad(s) = (T.array().square() * lamT).sum();
    }
    ArrayXd x_squared = X.colwise().squaredNorm().array();
    return -0.5 * (quad - log_det - x_squared);
}

// Samples from the standardised copula field: deterministic per-eigenpair
// seeding (src/rng_seed.h), single-threaded, row-major (i, j) accumulation.
// [[Rcpp::export]]
Eigen::MatrixXd rmatern_copula_eigen_cpp(int n, int dim_x, int dim_y, double rho1,
                                         double rho2, int nu,
                                         Rcpp::Nullable<int> seed = R_NilValue) {
    if (n < 0 || dim_x <= 0 || dim_y <= 0) return MatrixXd(0, 0);
    SparseMatrix<double> Q1 = make_AR_prec_matrix_cpp(dim_x, rho1);
    SparseMatrix<double> Q2 = make_AR_prec_matrix_cpp(dim_y, rho2);

    SelfAdjointEigenSolver<SparseMatrix<double>> solver1(Q1);
    VectorXd A1 = solver1.eigenvalues();
    MatrixXd V1 = solver1.eigenvectors();

    SelfAdjointEigenSolver<SparseMatrix<double>> solver2(Q2);
    VectorXd A2 = solver2.eigenvalues();
    MatrixXd V2 = solver2.eigenvectors();

    ArrayXd marginal_sds = marginal_sd_eigen(A1, V1, dim_x, A2, V2, dim_y, nu);

    int D = dim_x * dim_y;
    MatrixXd samples = MatrixXd::Zero(D, n);
    const std::uint64_t rng_base = stdmatern::seed_base(seed);

    VectorXd v(D);
    for (int i = 0; i < dim_x; ++i) {
        for (int j = 0; j < dim_y; ++j) {
            std::mt19937 gen = stdmatern::pair_generator(rng_base, i, j, dim_y);
            std::normal_distribution<double> d(0.0, 1.0);
            double lambda = std::pow(A1(i) + A2(j), -(nu + 1.0) / 2.0);
            v = kroneckerProduct(V1.col(i), V2.col(j));
            for (int s = 0; s < n; ++s) {
                samples.col(s) += lambda * d(gen) * v;
            }
        }
    }
    if (n > 0) samples = samples.array().colwise() / marginal_sds;
    return samples;
}

// ============================================================================
// MLE of the copula parameter(s): BOBYQA on [0.001, 0.999] through nlopt's C
// API. No wall-clock cap (the 120 s maxtime of the 2024 code made the returned
// rho load-dependent past ~16k stations, review F1); the deterministic budget
// is maxeval. Evaluations are counted here; the status contract is
// converged = code in 1..4, stop_reason in {converged, maxeval, roundoff, failed}.
// ============================================================================

struct OptData {
    const MatrixXd* Y;
    int dim_x;
    int dim_y;
    int nu;
    bool kron;
    int n_evals;
};

static double nll_at(OptData* d, double rho1, double rho2) {
    d->n_evals += 1;
    ArrayXd ld = d->kron
        ? dmatern_copula_kron_cpp(*(d->Y), d->dim_x, d->dim_y, rho1, rho2, d->nu)
        : dmatern_copula_eigen_cpp(*(d->Y), d->dim_x, d->dim_y, rho1, rho2, d->nu);
    double nll = -ld.sum();
    if (!std::isfinite(nll)) return 1e15;   // keeps BOBYQA off non-finite ground
    return nll;
}

static double nll_objective(unsigned n, const double* x, double* grad, void* data) {
    (void) n; (void) grad;
    return nll_at(static_cast<OptData*>(data), x[0], x[1]);
}

static double nll_objective_isotropic(unsigned n, const double* x, double* grad, void* data) {
    (void) n; (void) grad;
    return nll_at(static_cast<OptData*>(data), x[0], x[0]);
}

struct FitResult {
    int code;
    double min_nll;
    std::string message;
    std::string stop_reason;
    bool converged;
};

static FitResult run_bobyqa(OptData& data, std::vector<double>& x, nlopt_func f,
                            double xtol_rel, double ftol_rel, int maxeval) {
    FitResult r;
    r.min_nll = NA_REAL;
    const unsigned n = static_cast<unsigned>(x.size());
    nlopt_opt opt = nlopt_create(NLOPT_LN_BOBYQA, n);
    if (opt == NULL) {
        r.code = -3; r.message = "Optimization failed: nlopt_create returned NULL";
        r.stop_reason = "failed"; r.converged = false;
        return r;
    }
    nlopt_set_lower_bounds1(opt, 0.001);
    nlopt_set_upper_bounds1(opt, 0.999);
    nlopt_set_min_objective(opt, f, &data);
    nlopt_set_xtol_rel(opt, xtol_rel);
    nlopt_set_ftol_rel(opt, ftol_rel);
    nlopt_set_maxeval(opt, maxeval);
    double min_nll = NA_REAL;
    nlopt_result res = nlopt_optimize(opt, x.data(), &min_nll);
    nlopt_destroy(opt);
    r.code = static_cast<int>(res);
    r.min_nll = min_nll;
    switch (res) {
        case NLOPT_SUCCESS:          r.message = "Optimization succeeded"; break;
        case NLOPT_STOPVAL_REACHED:  r.message = "Stopval reached"; break;
        case NLOPT_FTOL_REACHED:     r.message = "Function tolerance reached"; break;
        case NLOPT_XTOL_REACHED:     r.message = "Parameter tolerance reached"; break;
        case NLOPT_MAXEVAL_REACHED:  r.message = "Maximum evaluations reached"; break;
        case NLOPT_MAXTIME_REACHED:  r.message = "Maximum time reached"; break;
        case NLOPT_ROUNDOFF_LIMITED: r.message = "Roundoff limited"; break;
        case NLOPT_FORCED_STOP:      r.message = "Forced stop"; break;
        case NLOPT_INVALID_ARGS:     r.message = "Optimization failed: invalid arguments"; break;
        case NLOPT_OUT_OF_MEMORY:    r.message = "Optimization failed: out of memory"; break;
        default:                     r.message = "Optimization failed";
    }
    r.converged = (r.code >= 1 && r.code <= 4);
    r.stop_reason = r.converged ? "converged"
                  : r.code == NLOPT_MAXEVAL_REACHED ? "maxeval"
                  : r.code == NLOPT_MAXTIME_REACHED ? "maxtime"
                  : r.code == NLOPT_ROUNDOFF_LIMITED ? "roundoff"
                  : "failed";
    return r;
}

// [[Rcpp::export]]
Rcpp::List fit_matern_copula_eigen_cpp(
    const Eigen::MatrixXd& Y,
    int dim_x, int dim_y, int nu,
    Rcpp::Nullable<Eigen::VectorXd> start = R_NilValue,
    bool kron = false, int maxeval = 1000, double xtol_rel = 1e-6, double ftol_rel = 1e-8
) {
    OptData opt_data;
    opt_data.Y = &Y;
    opt_data.dim_x = dim_x;
    opt_data.dim_y = dim_y;
    opt_data.nu = nu;
    opt_data.kron = kron;
    opt_data.n_evals = 0;

    std::vector<double> x(2, 0.5);
    if (start.isNotNull()) {
        VectorXd start_vec = Rcpp::as<VectorXd>(start);
        x[0] = start_vec(0);
        x[1] = start_vec(1);
    }
    FitResult r = run_bobyqa(opt_data, x, nll_objective, xtol_rel, ftol_rel, maxeval);
    return Rcpp::List::create(
        Rcpp::Named("rho1") = x[0],
        Rcpp::Named("rho2") = x[1],
        Rcpp::Named("loglik") = std::isfinite(r.min_nll) ? -r.min_nll : NA_REAL,
        Rcpp::Named("convergence") = r.code,
        Rcpp::Named("message") = r.message,
        Rcpp::Named("n_evals") = opt_data.n_evals,
        Rcpp::Named("converged") = r.converged,
        Rcpp::Named("stop_reason") = r.stop_reason,
        Rcpp::Named("method") = std::string(kron ? "kronecker" : "eigenpair")
    );
}

// [[Rcpp::export]]
Rcpp::List fit_matern_copula_eigen_isotropic_cpp(
    const Eigen::MatrixXd& Y,
    int dim_x, int dim_y, int nu,
    Rcpp::Nullable<double> start = R_NilValue,
    bool kron = false, int maxeval = 1000, double xtol_rel = 1e-6, double ftol_rel = 1e-8
) {
    OptData opt_data;
    opt_data.Y = &Y;
    opt_data.dim_x = dim_x;
    opt_data.dim_y = dim_y;
    opt_data.nu = nu;
    opt_data.kron = kron;
    opt_data.n_evals = 0;

    std::vector<double> x(1, 0.5);
    if (start.isNotNull()) x[0] = Rcpp::as<double>(start);
    FitResult r = run_bobyqa(opt_data, x, nll_objective_isotropic, xtol_rel, ftol_rel, maxeval);
    return Rcpp::List::create(
        Rcpp::Named("rho") = x[0],
        Rcpp::Named("loglik") = std::isfinite(r.min_nll) ? -r.min_nll : NA_REAL,
        Rcpp::Named("convergence") = r.code,
        Rcpp::Named("message") = r.message,
        Rcpp::Named("n_evals") = opt_data.n_evals,
        Rcpp::Named("converged") = r.converged,
        Rcpp::Named("stop_reason") = r.stop_reason,
        Rcpp::Named("method") = std::string(kron ? "kronecker" : "eigenpair")
    );
}
