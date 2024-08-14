#include <RcppEigen.h>
#include <complex>
#include <random>
#include <fftw3.h>
#include "circulant_utils.h"

// [[Rcpp::depends(RcppEigen)]]

using namespace Rcpp;
using namespace Eigen;

// [[Rcpp::export]]
Eigen::VectorXd fold_data(const Eigen::VectorXd& X, int n) {
    Eigen::VectorXd folded(4 * n * n);
    Eigen::Map<const Eigen::MatrixXd> X_mat(X.data(), n, n);
    
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            // First quarter
            folded[i * 2*n + j] = X_mat(i, j);
            
            // Second quarter
            folded[i * 2*n + (2*n - 1 - j)] = X_mat(i, j);
            
            // Third quarter
            folded[(2*n - 1 - i) * 2*n + j] = X_mat(i, j);
            
            // Fourth quarter
            folded[(2*n - 1 - i) * 2*n + (2*n - 1 - j)] = X_mat(i, j);
        }
    }
    
    return folded;
}

// [[Rcpp::export]]
Eigen::VectorXd dmatern_copula_folded(const Eigen::MatrixXd& X, int dim, double rho, int nu) {
    int n = dim;  // original dimension
    int N = 4 * n * n;  // folded dimension
    int n_obs = X.cols();
    Eigen::VectorXd quad_forms(n_obs);
    
    // Create base matrix (unchanged)
    Eigen::MatrixXd c = create_base_matrix(2 * n, rho);
    c /= 2.0;
    
    // Compute eigenvalues (unchanged)
    Eigen::MatrixXcd eigs = compute_and_rescale_eigenvalues(c, nu);
    
    // Compute log determinant (unchanged)
    double log_det = eigs.real().array().log().sum();
    
    // Compute quadratic form x^T Q x for each observation
    for (int i = 0; i < n_obs; ++i) {
        // Fold the data
        Eigen::VectorXd folded_X = fold_data(X.col(i), n);
        
        // Compute Qx using the folded data
        Eigen::VectorXd Qx = matrix_vector_product(eigs, folded_X);
        
        // Compute quadratic form
        quad_forms(i) = (folded_X.dot(Qx) - folded_X.squaredNorm());
    }
    
    // Compute log density
    Eigen::VectorXd log_densities = -0.5 * (quad_forms.array() - log_det);
    
    return log_densities;
}

// [[Rcpp::export]]
Eigen::MatrixXd rmatern_copula_folded_full(int n_samples, int dim, double rho, int nu) {
    int n = dim;  // original dimension
    int N = 2 * n;  // folded dimension
    Eigen::MatrixXd samples(n * n, n_samples);

    // Create base matrix for 2n x 2n grid
    Eigen::MatrixXd c = create_base_matrix(N, rho);
    c /= 2.0;  // Divide by 2 as in dmatern_copula_folded

    // Compute eigenvalues
    Eigen::MatrixXcd eigs = compute_and_rescale_eigenvalues(c, nu);
    Eigen::MatrixXcd lambda_sqrt = eigs.array().sqrt().inverse();

    // Set up random number generation
    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<> d(0, 1);
    double scale = 1.0 / (2 * dim);

    for(int s = 0; s < n_samples; ++s) {
        // Generate complex normal random field z of size 2n x 2n
        Eigen::MatrixXcd z(N, N);
        for(int i = 0; i < N; ++i) {
            for(int j = 0; j < N; ++j) {
                z(i, j) = std::complex<double>(d(gen), d(gen));
            }
        }

        // Perform FFT, multiply by sqrt(inv(eigenvalues)), and inverse FFT
        Eigen::MatrixXcd v = fft2(lambda_sqrt.cwiseProduct(z), scale);

        // Extract the n x n subset corresponding to the original data
        Eigen::MatrixXd sample = v.topLeftCorner(n, n).real();

        // Store the sample
        samples.col(s) = Eigen::Map<Eigen::VectorXd>(sample.data(), n * n);
    }

    return samples;
}