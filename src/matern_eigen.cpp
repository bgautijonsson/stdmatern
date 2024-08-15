#include <RcppEigen.h>
#include <omp.h>
#include <random>
#include "ar_matrix.h"


// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::plugins(openmp)]]

using namespace Rcpp;
using namespace Eigen;

// [[Rcpp::export]]
Eigen::VectorXd dmatern_eigen(const Eigen::MatrixXd& X, int dim_x, int dim_y, double rho1, double rho2, int nu) {
    int n_obs = X.cols();
    int N = dim_x * dim_y;
    Eigen::VectorXd quadform_sums = Eigen::VectorXd::Zero(n_obs);
    const double C = N * std::log(2 * M_PI);

    // Create precision matrices
    Eigen::MatrixXd Q1 = make_AR_prec_matrix(dim_x, rho1);
    Eigen::MatrixXd Q2 = make_AR_prec_matrix(dim_y, rho2);

    // Perform eigendecompositions
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver1(Q1);
    Eigen::VectorXd A1 = solver1.eigenvalues();
    Eigen::MatrixXd V1 = solver1.eigenvectors();

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver2(Q2);
    Eigen::VectorXd A2 = solver2.eigenvalues();
    Eigen::MatrixXd V2 = solver2.eigenvectors();


    double log_det = 0;

    #pragma omp parallel
    {
        Eigen::VectorXd local_quadform_sums = Eigen::VectorXd::Zero(n_obs);

        #pragma omp for reduction(+:log_det)
        for (int i = 0; i < dim_y; ++i) {
            for (int j = 0; j < dim_x; ++j) {
                // First calculate the Kronecker product
                Eigen::VectorXd v = Eigen::kroneckerProduct(V1.col(j), V2.col(i));
                
                // Compute the eigenvalue
                double lambda = std::pow(A2(i) + A1(j), nu + 1);

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
Eigen::MatrixXd rmatern_eigen(int n, int dim_x, int dim_y, double rho1, double rho2, int nu) {
    // Step 1: Create 1D AR(1) precision matrices
    Eigen::SparseMatrix<double> Q1 = make_AR_prec_matrix(dim_x, rho1);
    Eigen::SparseMatrix<double> Q2 = make_AR_prec_matrix(dim_y, rho2);
    
    // Step 2: Perform eigendecomposition of Q1 and Q2
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver1(Q1);
    Eigen::VectorXd A1 = solver1.eigenvalues();
    Eigen::MatrixXd V1 = solver1.eigenvectors();

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver2(Q2);
    Eigen::VectorXd A2 = solver2.eigenvalues();
    Eigen::MatrixXd V2 = solver2.eigenvectors();

    // Step 3: Prepare for sampling
    int D = dim_x * dim_y;
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

        for (int i = 0; i < dim_x; ++i) {
            for (int j = 0; j < dim_y; ++j) {
                Eigen::VectorXd v = Eigen::kroneckerProduct(V2.col(j), V1.col(i));
                double lambda = std::pow(A1(i) + A2(j), -(nu + 1.0) / 2.0);
                x += lambda * d(generators[thread_id]) * v;
            }
        }

        samples.col(s) = x;
    }

    return samples;
}