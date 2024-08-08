#include <RcppEigen.h>
#include <complex>
#include <random>
#include <fftw3.h>

// [[Rcpp::depends(RcppEigen, fftwtools)]]

using namespace Rcpp;
using namespace Eigen;

Eigen::MatrixXcd fft2(const Eigen::MatrixXcd& input, double scale) {
    int dim = input.rows();
    Eigen::MatrixXcd output(dim, dim);

    fftw_complex *in, *out;
    in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * dim * dim);
    out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * dim * dim);

    // Copy data from Eigen matrix to FFTW input
    for (int i = 0; i < dim; ++i) {
        for (int j = 0; j < dim; ++j) {
            in[i*dim + j][0] = input(i, j).real();
            in[i*dim + j][1] = input(i, j).imag();
        }
    }

    // Forward FFT
    fftw_plan plan = fftw_plan_dft_2d(dim, dim, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(plan);

    // Copy result back to Eigen matrix and scale by 1/(dim*dim)
    for (int i = 0; i < dim; ++i) {
        for (int j = 0; j < dim; ++j) {
            output(i, j) = std::complex<double>(out[i*dim + j][0], out[i*dim + j][1]) * scale;
        }
    }

    fftw_destroy_plan(plan);
    fftw_free(in);
    fftw_free(out);

    return output;
}

// Function to create the base matrix c
// [[Rcpp::export]]
Eigen::MatrixXd create_base_matrix(int dim, double rho) {
    Eigen::MatrixXd c = Eigen::MatrixXd::Zero(dim, dim);
    
    // Set the first row
    c(0, 0) = 2 + 2 * rho * rho;
    c(0, 1) = -rho;
    c(0, dim - 1) = -rho;
    
    // Set the second and last row
    c(1, 0) = -rho;
    c(dim - 1, 0) = -rho;
    
    return c;
}

// Function to compute FFT of c, apply nu, and rescale eigenvalues using marginal variance
// [[Rcpp::export]]
Eigen::MatrixXcd compute_and_rescale_eigenvalues(const Eigen::MatrixXd& c, int nu) {
    int dim = c.rows();

    fftw_complex *in, *out;
    in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * dim * dim);
    out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * dim * dim);

    // Copy data from Eigen matrix to FFTW input
    for (int i = 0; i < dim; ++i) {
        for (int j = 0; j < dim; ++j) {
            in[i*dim + j][0] = c(i, j);
            in[i*dim + j][1] = 0.0;
        }
    }

    // Forward FFT
    fftw_plan plan_forward = fftw_plan_dft_2d(dim, dim, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(plan_forward);

    // Copy eigenvalues to Eigen matrix and apply nu
    Eigen::MatrixXcd eigenvalues(dim, dim);
    for (int i = 0; i < dim; ++i) {
        for (int j = 0; j < dim; ++j) {
            std::complex<double> eig(out[i*dim + j][0], out[i*dim + j][1]);
            eigenvalues(i, j) = std::pow(eig, nu + 1);
        }
    }

    // Compute inverse of eigenvalues
    Eigen::MatrixXcd inv_eigenvalues = eigenvalues.cwiseInverse();

    // Inverse FFT of inverse eigenvalues
    for (int i = 0; i < dim * dim; ++i) {
        in[i][0] = inv_eigenvalues(i / dim, i % dim).real();
        in[i][1] = inv_eigenvalues(i / dim, i % dim).imag();
    }

    fftw_plan plan_backward = fftw_plan_dft_2d(dim, dim, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_execute(plan_backward);

    // Get marginal variance (first element of inverse FFT result)
    double mvar = out[0][0] / (dim * dim);

    // Scale eigenvalues
    Eigen::MatrixXcd updated_eigenvalues = mvar * eigenvalues;

    fftw_destroy_plan(plan_forward);
    fftw_destroy_plan(plan_backward);
    fftw_free(in);
    fftw_free(out);

    return updated_eigenvalues;
}

// Function for matrix-vector product using FFT and eigenvalues
// [[Rcpp::export]]
Eigen::VectorXd matrix_vector_product(const Eigen::MatrixXcd& eigenvalues, const Eigen::VectorXd& v) {
    int dim = eigenvalues.rows();

    fftw_complex *in, *out;
    in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * dim * dim);
    out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * dim * dim);

    // Copy data from Eigen vector to FFTW input
    for (int i = 0; i < dim * dim; ++i) {
        in[i][0] = v(i);
        in[i][1] = 0.0;
    }

    // Forward FFT
    fftw_plan plan_forward = fftw_plan_dft_2d(dim, dim, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(plan_forward);

    // Multiply with eigenvalues in frequency domain
    for (int i = 0; i < dim * dim; ++i) {
        std::complex<double> temp(out[i][0], out[i][1]);
        temp *= eigenvalues(i / dim, i % dim);
        out[i][0] = temp.real();
        out[i][1] = temp.imag();
    }

    // Inverse FFT
    fftw_plan plan_backward = fftw_plan_dft_2d(dim, dim, out, in, FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_execute(plan_backward);

    // Copy result back to Eigen vector and normalize
    Eigen::VectorXd result(dim * dim);
    for (int i = 0; i < dim * dim; ++i) {
        result(i) = in[i][0] / (dim * dim);
    }

    fftw_destroy_plan(plan_forward);
    fftw_destroy_plan(plan_backward);
    fftw_free(in);
    fftw_free(out);

    return result;
}

// Function for log-density calculation
// [[Rcpp::export]]
Eigen::VectorXd dmatern_copula_circulant(const Eigen::MatrixXd& X, int dim, double rho, int nu) {
    Eigen::MatrixXd c = create_base_matrix(dim, rho);
    Eigen::MatrixXcd eigs = compute_and_rescale_eigenvalues(c, nu);
    int N = dim * dim;
    int n_obs = X.cols();
    Eigen::VectorXd quad_forms(n_obs);
    
    // Compute log determinant
    double log_det = eigs.real().array().log().sum();
    
    // Compute quadratic form x^T Q x
    for (int i = 0; i < n_obs; ++i) {
        Eigen::VectorXd Qx = matrix_vector_product(eigs, X.col(i));
        quad_forms(i) = X.col(i).dot(Qx);
    }
    
    // Subtract x'Ix (sum of squared elements for each observation)
    Eigen::VectorXd x_squared = X.colwise().squaredNorm();

    // Compute log density
    Eigen::VectorXd log_densities = -0.5 * (quad_forms.array() - log_det - x_squared.array());
    
    return log_densities;
}

// [[Rcpp::export]]
Eigen::MatrixXd rmatern_copula_circulant(int n_samples, int dim, double rho, int nu) {
    Eigen::MatrixXd samples(dim * dim, n_samples);

     // Step 1: Compute the (real) eigenvalues, Λ = √(nN) DFT2(θ)
    Eigen::MatrixXd c = create_base_matrix(dim, rho);
    Eigen::MatrixXcd eigs = compute_and_rescale_eigenvalues(c, nu);
    Eigen::MatrixXcd lambda_sqrt = eigs.array().sqrt().inverse();

    // Step 1: Sample z, where Re(z_ij) ~ N(0,1) and Im(z_ij) ~ N(0,1) iid
    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<> d(0, 1);
    Eigen::MatrixXcd z(dim, dim);
    double scale = 1.0 / (dim);
    
    for(int s = 0; s < n_samples; ++s) {

        for(int i = 0; i < dim; ++i) {
            for(int j = 0; j < dim; ++j) {
                z(i, j) = std::complex<double>(d(gen), d(gen));
            }
        }
        
        Eigen::MatrixXcd v = fft2(lambda_sqrt.cwiseProduct(z), scale);
        
        // Step 4: x = Re(v)
        samples.col(s) = v.real().transpose().reshaped(v.size(), 1);

    }
    
    // Step 5: Return x
    return samples;
}


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
        quad_forms(i) = folded_X.dot(Qx) - folded_X.squaredNorm();
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