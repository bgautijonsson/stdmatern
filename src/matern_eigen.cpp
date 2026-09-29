#include <RcppEigen.h>
#include <random>
#include "ar_matrix.h"
#include "rng_seed.h"

// [[Rcpp::depends(RcppEigen)]]

using namespace Rcpp;
using namespace Eigen;

// Log-density of the (unstandardised) Matern-like field with precision
// Q = (Q1 (x) I + I (x) Q2)^(nu + 1), via the two AR(1) eigen systems.
// Single-threaded since 2026-09-02 (review F2). Inputs are validated in R;
// a row mismatch returns NA.
// [[Rcpp::export]]
Eigen::VectorXd dmatern_eigen_cpp(const Eigen::MatrixXd& X, int dim_x, int dim_y, double rho1, double rho2, int nu) {
    int n_obs = X.cols();
    int N = dim_x * dim_y;
    if (dim_x <= 0 || dim_y <= 0 || X.rows() != N) return VectorXd::Constant(n_obs, NA_REAL);
    ArrayXd quadform_sums = ArrayXd::Zero(n_obs);
    const double C = N * std::log(2 * M_PI);

    SparseMatrix<double> Q1 = make_AR_prec_matrix_cpp(dim_x, rho1);
    SparseMatrix<double> Q2 = make_AR_prec_matrix_cpp(dim_y, rho2);

    SelfAdjointEigenSolver<SparseMatrix<double>> solver1(Q1);
    VectorXd A1 = solver1.eigenvalues();
    MatrixXd V1 = solver1.eigenvectors();

    SelfAdjointEigenSolver<SparseMatrix<double>> solver2(Q2);
    VectorXd A2 = solver2.eigenvalues();
    MatrixXd V2 = solver2.eigenvectors();

    double log_det = 0;
    VectorXd v(N);
    for (int i = 0; i < dim_x; ++i) {
        for (int j = 0; j < dim_y; ++j) {
            v = kroneckerProduct(V1.col(i), V2.col(j));
            double lambda = std::pow(A1(i) + A2(j), nu + 1);
            log_det += std::log(lambda);
            for (int obs = 0; obs < n_obs; ++obs) {
                double u = v.dot(X.col(obs));
                quadform_sums(obs) += u * u * lambda;
            }
        }
    }

    ArrayXd log_densities = -0.5 * (C - log_det + quadform_sums);
    return log_densities;
}

// Samples from the unstandardised field; deterministic per-eigenpair seeding
// (src/rng_seed.h), single-threaded.
// [[Rcpp::export]]
Eigen::MatrixXd rmatern_eigen_cpp(int n, int dim_x, int dim_y, double rho1,
                                  double rho2, int nu,
                                  Rcpp::Nullable<int> seed = R_NilValue, bool legacy32 = false) {
    if (n < 0 || dim_x <= 0 || dim_y <= 0) return MatrixXd(0, 0);
    SparseMatrix<double> Q1 = make_AR_prec_matrix_cpp(dim_x, rho1);
    SparseMatrix<double> Q2 = make_AR_prec_matrix_cpp(dim_y, rho2);

    SelfAdjointEigenSolver<SparseMatrix<double>> solver1(Q1);
    VectorXd A1 = solver1.eigenvalues();
    MatrixXd V1 = solver1.eigenvectors();

    SelfAdjointEigenSolver<SparseMatrix<double>> solver2(Q2);
    VectorXd A2 = solver2.eigenvalues();
    MatrixXd V2 = solver2.eigenvectors();

    int D = dim_x * dim_y;
    MatrixXd samples = MatrixXd::Zero(D, n);
    const std::uint64_t rng_base = stdmatern::seed_base(seed);

    VectorXd v(D);
    for (int i = 0; i < dim_x; ++i) {
        for (int j = 0; j < dim_y; ++j) {
            std::mt19937 gen = stdmatern::pair_generator(rng_base, i, j, dim_y, legacy32);
            std::normal_distribution<double> d(0.0, 1.0);
            double lambda = std::pow(A1(i) + A2(j), -(nu + 1.0) / 2.0);
            v = kroneckerProduct(V1.col(i), V2.col(j));
            for (int s = 0; s < n; ++s) {
                samples.col(s) += lambda * d(gen) * v;
            }
        }
    }
    return samples;
}
