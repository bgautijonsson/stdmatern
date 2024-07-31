#include <RcppEigen.h>
#include <complex>
#include <unsupported/Eigen/FFT>

// [[Rcpp::depends(RcppEigen)]]

using namespace Rcpp;
using namespace Eigen;

// Function to create a folded circulant AR(1) precision matrix
Eigen::VectorXd make_folded_circulant_AR(int n, double rho) {
    Eigen::VectorXd base = Eigen::VectorXd::Zero(2 * n);
    double scaling = 1.0 / (1.0 - rho * rho);
    base[0] = (1.0 + rho * rho) * scaling;
    base[1] = base[2 * n - 1] = -rho * scaling;
    return base;
}

// Function to compute eigenvalues of a folded circulant matrix
Eigen::VectorXcd compute_folded_circulant_eigenvalues(const Eigen::VectorXd& folded_circulant) {
    Eigen::FFT<double> fft;
    return fft.fwd(folded_circulant);
}

// Function to expand the original data to the folded structure
Eigen::VectorXd expand_data(const Eigen::VectorXd& X, int n) {
    Eigen::VectorXd X_expanded = Eigen::VectorXd::Zero(4 * n * n);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            int idx = i * n + j;
            X_expanded[4 * idx] = X_expanded[4 * idx + 3] = X[i * n + j];
            X_expanded[4 * idx + 1] = X_expanded[4 * idx + 2] = X[(n - 1 - i) * n + j];
        }
    }
    return X_expanded;
}

// [[Rcpp::export]]
List dmatern_copula_folded(const Eigen::MatrixXd& X, int n, double rho, int nu) {
    int N = n * n;
    int n_obs = X.cols();
    
    // Step 1: Create folded circulant matrix C1
    Eigen::VectorXd C1 = make_folded_circulant_AR(n, rho);

    // Step 2: Compute eigenvalues of C1
    Eigen::VectorXcd eig_C1 = compute_folded_circulant_eigenvalues(C1);

    // Step 3: Compute marginal standard deviations
    Eigen::MatrixXcd eig_C(2*n, 2*n);
    for (int i = 0; i < 2*n; ++i) {
        for (int j = 0; j < 2*n; ++j) {
            eig_C(i, j) = eig_C1(i) + eig_C1(j);
        }
    }
    if (nu > 0) {
        eig_C = eig_C.array().pow(nu + 1);
    }
    
    Eigen::VectorXd marginal_sds = Eigen::VectorXd::Zero(4*N);
    for (int i = 0; i < 2*n; ++i) {
        for (int j = 0; j < 2*n; ++j) {
            marginal_sds(i*2*n + j) = 1.0 / std::sqrt(eig_C(i, j).real());
        }
    }

    // Step 4: Update spectral decomposition
    eig_C = eig_C * marginal_sds.array().square().matrix().asDiagonal();

    // Step 5: Compute densities
    Eigen::VectorXd log_densities(n_obs);
    Eigen::FFT<double> fft;

    for (int obs = 0; obs < n_obs; ++obs) {
        // Expand the data
        Eigen::VectorXd x_expanded = expand_data(X.col(obs), n);
        
        // Apply FFT
        Eigen::VectorXcd x_fft = fft.fwd(x_expanded);
        
        // Compute quadratic form and log determinant
        double quad_form = 0;
        double log_det = 0;
        for (int i = 0; i < 2*n; ++i) {
            for (int j = 0; j < 2*n; ++j) {
                int idx = i*2*n + j;
                quad_form += std::norm(x_fft[idx]) * eig_C(i, j).real();
                log_det += std::log(eig_C(i, j).real());
            }
        }
        
        // Compute log density
        log_densities(obs) = -0.5 * (4*N * std::log(2 * M_PI) + log_det + quad_form);
    }

    return List::create(
        Named("C1") = C1,
        Named("marginal_sds") = marginal_sds,
        Named("log_densities") = log_densities
    );
}

// [[Rcpp::export]]
Eigen::MatrixXd rmatern_copula_folded(int n_samples, int n, double rho, int nu) {
    int N = n * n;
    
    // Step 1: Create folded circulant matrix C1
    Eigen::VectorXd C1 = make_folded_circulant_AR(n, rho);

    // Step 2: Compute eigenvalues of C1
    Eigen::VectorXcd eig_C1 = compute_folded_circulant_eigenvalues(C1);

    // Step 3: Compute marginal standard deviations
    Eigen::MatrixXcd eig_C(2*n, 2*n);
    for (int i = 0; i < 2*n; ++i) {
        for (int j = 0; j < 2*n; ++j) {
            eig_C(i, j) = eig_C1(i) + eig_C1(j);
        }
    }
    if (nu > 0) {
        eig_C = eig_C.array().pow(nu + 1);
    }
    
    Eigen::VectorXd marginal_sds = Eigen::VectorXd::Zero(4*N);
    for (int i = 0; i < 2*n; ++i) {
        for (int j = 0; j < 2*n; ++j) {
            marginal_sds(i*2*n + j) = 1.0 / std::sqrt(eig_C(i, j).real());
        }
    }

    // Step 4: Update spectral decomposition
    eig_C = eig_C.array().sqrt();

    // Step 5: Generate samples
    Eigen::MatrixXd samples(N, n_samples);
    Eigen::FFT<double> fft;

    for (int s = 0; s < n_samples; ++s) {
        Eigen::VectorXcd z(4*N);
        for (int i = 0; i < 4*N; ++i) {
            double re = R::rnorm(0, 1);
            double im = R::rnorm(0, 1);
            z(i) = std::complex<double>(re, im);
        }

        for (int i = 0; i < 2*n; ++i) {
            for (int j = 0; j < 2*n; ++j) {
                int idx = i*2*n + j;
                z(idx) *= std::sqrt(eig_C(i, j).real()) * marginal_sds(idx);
            }
        }

        Eigen::VectorXd x = fft.inv(z).real();
        samples.col(s) = x.head(N);
    }

    return samples;
}