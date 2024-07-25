#include <RcppEigen.h>
#include <omp.h>
#include <random>
#include "ar_matrix.h"


// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::plugins(openmp)]]

using namespace Rcpp;
using namespace Eigen;


// [[Rcpp::export]]
Eigen::VectorXd dmatern_eigen(const Eigen::MatrixXd& X, int dim, double rho, int nu) {
    int n_obs = X.cols();
    int N = dim * dim;
    Eigen::VectorXd quadform_sums = Eigen::VectorXd::Zero(n_obs);
    const double C = N * std::log(2 * M_PI);

    // Create precision matrix
    Eigen::MatrixXd Q1 = make_AR_prec_matrix(dim, rho);

    // Perform eigendecomposition
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(Q1);
    Eigen::VectorXd A1 = solver.eigenvalues();
    Eigen::MatrixXd V1 = solver.eigenvectors();


    double log_det = 0;

    #pragma omp parallel
    {
        Eigen::VectorXd local_quadform_sums = Eigen::VectorXd::Zero(n_obs);

        #pragma omp for reduction(+:log_det)
        for (int i = 0; i < dim; ++i) {
            for (int j = 0; j < dim; ++j) {
                // First calculate the Kronecker product
                Eigen::VectorXd v = Eigen::kroneckerProduct(V1.col(j), V1.col(i));
                
                // Compute the eigenvalue
                double lambda = (nu == 0) ? (A1(i) + A1(j)) : std::pow(A1(i) + A1(j), nu + 1);

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

    Eigen::VectorXd log_densities = -0.5 * (C - log_det + quadform_sums.array());

    return log_densities;
}

// [[Rcpp::export]]
Eigen::MatrixXd rmatern(int n, int dim, double rho, int nu) {
    // Step 1: Create 1D AR(1) precision matrix
    Eigen::SparseMatrix<double> Q1 = make_AR_prec_matrix(dim, rho);
    
    // Step 2: Perform eigendecomposition of Q1
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(Q1);
    Eigen::VectorXd A1 = solver.eigenvalues();
    Eigen::MatrixXd V1 = solver.eigenvectors();

    // Step 3: Prepare for sampling
    int D = dim * dim;
    Eigen::MatrixXd samples(D, n);

    // Random number generation
    std::vector<std::mt19937> generators(omp_get_max_threads());
    for (auto& gen : generators) {
        std::random_device rd;
        gen.seed(rd());
    }
    std::normal_distribution<> d(0, 1);

    // Step 5: Generate samples
    #pragma omp parallel for
    for (int s = 0; s < n; ++s) {
        int thread_id = omp_get_thread_num();
        Eigen::VectorXd x = Eigen::VectorXd::Zero(D);

        for (int i = 0; i < dim; ++i) {
            for (int j = 0; j < dim; ++j) {
                Eigen::VectorXd v = Eigen::kroneckerProduct(V1.col(j), V1.col(i));
                double lambda = std::pow(A1(i) + A1(j), -(nu + 1.0) / 2.0);
                x += lambda * d(generators[thread_id]) * v;
            }
        }

        samples.col(s) = x;
    }

    return samples;
}