#include <RcppEigen.h>
#include <omp.h>
#include <random>


// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::plugins(openmp)]]

using namespace Rcpp;
using namespace Eigen;


// Function to create a 1-dimensional AR(1) precision matrix
// [[Rcpp::export]]
Eigen::SparseMatrix<double> make_AR_prec_matrix(int dim, double rho) {
    double scaling = 1.0 / (1.0 - rho * rho);
    double off_diags = -rho * scaling;
    double diag = (1.0 + rho * rho) * scaling;

    Eigen::SparseMatrix<double> Q(dim, dim);
    Q.reserve(Eigen::VectorXi::Constant(dim, 3));  // Reserve space for 3 non-zeros per column

    for (int i = 0; i < dim; ++i) {
        Q.insert(i, i) = (i == 0 || i == dim - 1) ? scaling : diag;
        if (i > 0) Q.insert(i, i-1) = off_diags;
        if (i < dim - 1) Q.insert(i, i+1) = off_diags;
    }

    Q.makeCompressed();
    return Q;
}

// [[Rcpp::export]]
Eigen::VectorXd marginal_sd_eigen(const Eigen::VectorXd& A1, const Eigen::MatrixXd& V1, int dim, int nu) {
    Eigen::VectorXd marginal_sds = Eigen::VectorXd::Zero(dim * dim);
    
    for (int i = 0; i < dim; ++i) {
        for (int j = 0; j < dim; ++j) {
            Eigen::VectorXd v = Eigen::kroneckerProduct(V1.col(i), V1.col(j));
            double lambda = (nu == 0) ? (A1(i) + A1(j)) : std::pow(A1(i) + A1(j), nu + 1);
            marginal_sds += (v.array().square() / lambda).matrix();
        }
    }

    return marginal_sds.array().sqrt();
}

// Function to create the standardized Matérn precision matrix with eigendecomposition method
// [[Rcpp::export]]
Eigen::SparseMatrix<double> make_standardized_matern_eigen(int dim, double rho, int nu) {
    // Step 1: Create unstandardized Matérn precision matrix
    Eigen::SparseMatrix<double> Q1 = make_AR_prec_matrix(dim, rho);
    
    // Perform eigendecomposition of Q1
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(Q1);
    Eigen::VectorXd A1 = solver.eigenvalues();
    Eigen::MatrixXd V1 = solver.eigenvectors();

    // Step 2: Calculate marginal standard deviations
    Eigen::VectorXd marginal_sds = marginal_sd_eigen(A1, V1, dim, nu);

    // Step 3: Create diagonal matrix of inverse standard deviations
    Eigen::SparseMatrix<double> D(dim * dim, dim * dim);
    D.reserve(Eigen::VectorXi::Constant(dim * dim, 1));
    for (int i = 0; i < dim * dim; ++i) {
        D.insert(i, i) = marginal_sds(i);
    }
    D.makeCompressed();

    // Create full Q matrix using Kronecker sum
    Eigen::SparseMatrix<double> I(dim, dim);
    I.setIdentity();
    Eigen::SparseMatrix<double> Q = Eigen::kroneckerProduct(Q1, I) + Eigen::kroneckerProduct(I, Q1);

    // Apply matrix multiplication nu times
    if (nu > 0) {
        Eigen::SparseMatrix<double> temp = Q;
        for (int i = 0; i < nu; ++i) {
            Q = Q * temp;
            Q.makeCompressed();
        }
    }

    // Standardize Q
    Eigen::SparseMatrix<double> Q_standardized = D * Q * D;

    return Q_standardized;
}

// [[Rcpp::export]]
Eigen::VectorXd matern_mvn_density_eigen(const Eigen::MatrixXd& X, int dim, double rho, int nu) {
    int n_obs = X.cols();
    int N = dim * dim;
    Eigen::VectorXd log_densities = Eigen::VectorXd::Zero(n_obs);
    Eigen::VectorXd quadform_sums = Eigen::VectorXd::Zero(n_obs);
    const double C = N * std::log(2 * M_PI);

    // Create precision matrix
    Eigen::MatrixXd Q1 = make_AR_prec_matrix(dim, rho);

    // Perform eigendecomposition
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(Q1);
    Eigen::VectorXd A1 = solver.eigenvalues();
    Eigen::MatrixXd V1 = solver.eigenvectors();

    // Compute marginal standard deviations
    Eigen::VectorXd marginal_sds = marginal_sd_eigen(A1, V1, dim, nu);

    // Create diagonal matrix D
    Eigen::DiagonalMatrix<double, Eigen::Dynamic> D(N);
    D.diagonal() = marginal_sds;

    double log_det = 0;

    #pragma omp parallel
    {
        Eigen::VectorXd local_quadform_sums = Eigen::VectorXd::Zero(n_obs);

        #pragma omp for reduction(+:log_det)
        for (int i = 0; i < dim; ++i) {
            for (int j = 0; j < dim; ++j) {
                // First calculate the Kronecker product
                Eigen::VectorXd v = Eigen::kroneckerProduct(V1.col(j), V1.col(i));
                
                // Then scale it
                v = D * v;
                
                // Now calculate the norm
                double norm_v = v.norm();
                
                // Normalize v
                v /= norm_v;
                
                // Compute the eigenvalue
                double A = (nu == 0) ? (A1(i) + A1(j)) : std::pow(A1(i) + A1(j), nu + 1);
                
                // Standardize the eigenvalue
                double lambda = A * (norm_v * norm_v);

                log_det += std::log(lambda);
                
                // Loop over observations
                for (int obs = 0; obs < n_obs; ++obs) {
                    double u = v.dot(X.col(obs));
                    local_quadform_sums(obs) += u * u * lambda;
                }
            }
        }

        #pragma omp critical
        {
            quadform_sums += local_quadform_sums;
        }
    }
    
    // Compute final log densities
    #pragma omp parallel for
    for (int obs = 0; obs < n_obs; ++obs) {
        log_densities(obs) = -0.5 * (C - log_det + quadform_sums(obs));
    }

    return log_densities;
}

// [[Rcpp::export]]
Eigen::MatrixXd sample_standardized_matern(int dim, double rho, int nu, int n_samples) {
    // Step 1: Create 1D AR(1) precision matrix
    Eigen::SparseMatrix<double> Q1 = make_AR_prec_matrix(dim, rho);
    
    // Step 2: Perform eigendecomposition of Q1
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(Q1);
    Eigen::VectorXd A1 = solver.eigenvalues();
    Eigen::MatrixXd V1 = solver.eigenvectors();

    // Step 3: Compute marginal standard deviations
    Eigen::VectorXd marginal_sds = marginal_sd_eigen(A1, V1, dim, nu);

    // Step 4: Prepare for sampling
    int D = dim * dim;
    Eigen::MatrixXd samples(D, n_samples);
    
    // Random number generation
    std::vector<std::mt19937> generators(omp_get_max_threads());
    for (auto& gen : generators) {
        std::random_device rd;
        gen.seed(rd());
    }
    std::normal_distribution<> d(0, 1);

    // Step 5: Generate samples
    #pragma omp parallel for
    for (int s = 0; s < n_samples; ++s) {
        int thread_id = omp_get_thread_num();
        Eigen::VectorXd x = Eigen::VectorXd::Zero(D);

        for (int i = 0; i < dim; ++i) {
            for (int j = 0; j < dim; ++j) {
                Eigen::VectorXd v = Eigen::kroneckerProduct(V1.col(j), V1.col(i));
                double lambda = std::pow(A1(i) + A1(j), -(nu + 1.0) / 2.0);
                x += lambda * d(generators[thread_id]) * v;
            }
        }

        // Standardize the sample
        x = x.array() / marginal_sds.array();

        samples.col(s) = x;
    }

    return samples;
}