#include <RcppEigen.h>
#include <complex>
#include <random>
#include <fftw3.h>

// [[Rcpp::depends(RcppEigen, fftwtools)]]

using namespace Rcpp;
using namespace Eigen;

#include <fftw3.h>

// Function to create the base matrix c for DCT
// [[Rcpp::export]]
Eigen::MatrixXd create_base_matrix_dct(int dim, double rho) {
    Eigen::MatrixXd c = Eigen::MatrixXd::Zero(dim, dim);
    
    // Set the first row
    c(0, 0) = 2 + 2 * rho * rho;
    c(0, 1) = -2 * rho;
    
    // Set the second row
    c(1, 0) = -rho;
    
    return c;
}

// Function to compute DCT of c, apply nu, and rescale eigenvalues
// [[Rcpp::export]]
Eigen::MatrixXd compute_and_rescale_eigenvalues_dct(const Eigen::MatrixXd& c, int nu) {
    int dim = c.rows();

    double *in = (double*) fftw_malloc(sizeof(double) * dim * dim);
    double *out = (double*) fftw_malloc(sizeof(double) * dim * dim);

    // Copy data from Eigen matrix to FFTW input
    for (int i = 0; i < dim * dim; ++i) {
        in[i] = c(i / dim, i % dim);
    }

    // Forward DCT
    fftw_plan plan_forward = fftw_plan_r2r_2d(dim, dim, in, out, FFTW_REDFT10, FFTW_REDFT10, FFTW_ESTIMATE);
    fftw_execute(plan_forward);

    // Copy eigenvalues to Eigen matrix and apply nu
    Eigen::MatrixXd eigenvalues(dim, dim);
    for (int i = 0; i < dim * dim; ++i) {
        eigenvalues(i / dim, i % dim) = std::pow(out[i], nu + 1);
    }

    // Compute inverse of eigenvalues
    Eigen::MatrixXd inv_eigenvalues = eigenvalues.cwiseInverse();

    // Inverse DCT of inverse eigenvalues
    fftw_plan plan_backward = fftw_plan_r2r_2d(dim, dim, inv_eigenvalues.data(), in, FFTW_REDFT01, FFTW_REDFT01, FFTW_ESTIMATE);
    fftw_execute(plan_backward);

    // Get marginal variance (first element of inverse DCT result)
    double mvar = in[0] / (4 * dim * dim);

    // Scale eigenvalues
    Eigen::MatrixXd updated_eigenvalues = mvar * eigenvalues;

    fftw_destroy_plan(plan_forward);
    fftw_destroy_plan(plan_backward);
    fftw_free(in);
    fftw_free(out);

    return updated_eigenvalues;
}

// Function for matrix-vector product using DCT and eigenvalues
// [[Rcpp::export]]
Eigen::VectorXd matrix_vector_product_dct(const Eigen::MatrixXd& eigenvalues, const Eigen::VectorXd& v) {
    int dim = eigenvalues.rows();

    double *in = (double*) fftw_malloc(sizeof(double) * dim * dim);
    double *out = (double*) fftw_malloc(sizeof(double) * dim * dim);

    // Copy data from Eigen vector to FFTW input
    for (int i = 0; i < dim * dim; ++i) {
        in[i] = v(i);
    }

    // Forward DCT
    fftw_plan plan_forward = fftw_plan_r2r_2d(dim, dim, in, out, FFTW_REDFT10, FFTW_REDFT10, FFTW_ESTIMATE);
    fftw_execute(plan_forward);

    // Multiply with eigenvalues in frequency domain
    for (int i = 0; i < dim * dim; ++i) {
        out[i] *= eigenvalues(i / dim, i % dim);
    }

    // Inverse DCT
    fftw_plan plan_backward = fftw_plan_r2r_2d(dim, dim, out, in, FFTW_REDFT01, FFTW_REDFT01, FFTW_ESTIMATE);
    fftw_execute(plan_backward);

    // Copy result back to Eigen vector and normalize
    Eigen::VectorXd result(dim * dim);
    for (int i = 0; i < dim * dim; ++i) {
        result(i) = in[i] / (4 * dim * dim);
    }

    fftw_destroy_plan(plan_forward);
    fftw_destroy_plan(plan_backward);
    fftw_free(in);
    fftw_free(out);

    return result;
}

// Function for log-density calculation using DCT
// [[Rcpp::export]]
Eigen::VectorXd dmatern_copula_dct(const Eigen::MatrixXd& X, int dim, double rho, int nu) {
    Eigen::MatrixXd c = create_base_matrix_dct(dim, rho);
    Eigen::MatrixXd eigs = compute_and_rescale_eigenvalues_dct(c, nu);
    int N = dim * dim;
    int n_obs = X.cols();
    Eigen::VectorXd quad_forms(n_obs);
    
    // Compute log determinant
    double log_det = eigs.array().log().sum();
    
    // Compute quadratic form x^T Q x
    for (int i = 0; i < n_obs; ++i) {
        Eigen::VectorXd Qx = matrix_vector_product_dct(eigs, X.col(i));
        quad_forms(i) = X.col(i).dot(Qx);
    }
    
    // Subtract x'Ix (sum of squared elements for each observation)
    Eigen::VectorXd x_squared = X.colwise().squaredNorm();

    // Compute log density
    Eigen::VectorXd log_densities = -0.5 * (quad_forms.array() - log_det + x_squared.array());
    
    return log_densities;
}