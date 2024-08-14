#include <RcppEigen.h>
#include <complex>
#include <random>
#include <fftw3.h>
#include "circulant_utils.h"

// [[Rcpp::depends(RcppEigen)]]

using namespace Rcpp;
using namespace Eigen;

// [[Rcpp::export]]
Eigen::VectorXd fold_data(const Eigen::VectorXd& X, int n1, int n2) {
    Eigen::VectorXd folded(4 * n1 * n2);
    Eigen::Map<const Eigen::MatrixXd> X_mat(X.data(), n1, n2);
    
    for (int i = 0; i < n1; ++i) {
        for (int j = 0; j < n2; ++j) {
            // First quarter
            folded[i * 2*n2 + j] = X_mat(i, j);
            
            // Second quarter
            folded[i * 2*n2 + (2*n2 - 1 - j)] = X_mat(i, j);
            
            // Third quarter
            folded[(2*n1 - 1 - i) * 2*n2 + j] = X_mat(i, j);
            
            // Fourth quarter
            folded[(2*n1 - 1 - i) * 2*n2 + (2*n2 - 1 - j)] = X_mat(i, j);
        }
    }
    
    return folded;
}

// [[Rcpp::export]]
Eigen::VectorXd dmatern_copula_folded(const Eigen::MatrixXd& X, int dim1, int dim2, double rho1, double rho2, int nu) {
    int N = 4 * dim1 * dim2;  // folded dimension
    int n_obs = X.cols();
    Eigen::VectorXd quad_forms(n_obs);
    
    // Create base matrix
    Eigen::MatrixXd c = create_base_matrix(2 * dim1, 2 * dim2, rho1, rho2);
    
    // Compute eigenvalues
    Eigen::MatrixXcd eigs = compute_and_rescale_eigenvalues(c, nu);
    
    // Compute log determinant
    double log_det = eigs.array().log().sum().real();
    
    // Compute quadratic form x^T Q x for each observation
    for (int i = 0; i < n_obs; ++i) {
        // Fold the data
        Eigen::VectorXd folded_X = fold_data(X.col(i), dim1, dim2);
        
        // Compute Qx using the folded data
        Eigen::VectorXd Qx = matrix_vector_product(eigs, folded_X);
        
        // Compute quadratic form
        quad_forms(i) = (folded_X.dot(Qx) - folded_X.squaredNorm());
    }
    
    // Compute log density
    Eigen::VectorXd log_densities = - (quad_forms.array() - log_det) / 8.0;
    
    return log_densities;
}

// [[Rcpp::export]]
Eigen::MatrixXd rmatern_copula_folded_full(int n_samples, int dim1, int dim2, double rho1, double rho2, int nu) {
    int N1 = 2 * dim1;  // folded dimension
    int N2 = 2 * dim2;  // folded dimension
    Eigen::MatrixXd samples(dim1 * dim2, n_samples);

    // Create base matrix for 2dim1 x 2dim2 grid
    Eigen::MatrixXd c = create_base_matrix(N1, N2, rho1, rho2);

    // Compute eigenvalues
    Eigen::MatrixXcd eigs = compute_and_rescale_eigenvalues(c, nu);
    Eigen::MatrixXcd lambda_sqrt = eigs.array().sqrt().inverse();

    // Set up random number generation
    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<> d(0, 1);
    double scale = 1.0 / std::sqrt(N1 * N2);

    for(int s = 0; s < n_samples; ++s) {
        // Generate complex normal random field z of size N2 x N1
        Eigen::MatrixXcd z(N2, N1);
        for(int i = 0; i < N2; ++i) {
            for(int j = 0; j < N1; ++j) {
                z(i, j) = std::complex<double>(d(gen), d(gen));
            }
        }

        // Perform FFT, multiply by sqrt(inv(eigenvalues)), and inverse FFT
        Eigen::MatrixXcd v = fft2(lambda_sqrt.cwiseProduct(z), scale);

        // Extract the dim1 x dim2 subset corresponding to the original data
        Eigen::MatrixXd sample = v.topLeftCorner(dim2, dim1).real();

        // Store the sample
        samples.col(s) = Eigen::Map<Eigen::VectorXd>(sample.data(), dim1 * dim2);
    }

    return samples;
}