#include <RcppEigen.h>
#include <omp.h>
#include <random>
#include "ar_matrix.h"
#include <Rcpp/Benchmark/Timer.h>


// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::plugins(openmp)]]


using namespace Rcpp;
using namespace Eigen;



// [[Rcpp::export]]
Eigen::VectorXd marginal_sd_eigen(const Eigen::VectorXd& A1, const Eigen::MatrixXd& V1, int dim_x, const Eigen::VectorXd& A2, const Eigen::MatrixXd& V2, int dim_y, int nu) {
    Eigen::ArrayXd marginal_sds = Eigen::ArrayXd::Zero(dim_x * dim_y);
    
    #pragma omp parallel
    {
        Eigen::ArrayXd local_marginal_sds = Eigen::ArrayXd::Zero(dim_x * dim_y);
        Eigen::VectorXd v(dim_x * dim_y);

        #pragma omp for collapse(2) nowait
        for (int i = 0; i < dim_x; ++i) {
            for (int j = 0; j < dim_y; ++j) {
                v = Eigen::kroneckerProduct(V1.col(i), V2.col(j));
                local_marginal_sds += v.array().square() / std::pow(A1(i) + A2(j), nu + 1);
            }
        }

        #pragma omp critical
        {
            marginal_sds += local_marginal_sds;
        }
    }

    return marginal_sds.sqrt();
}

// Function to create the standardized Matérn precision matrix with eigendecomposition method
// [[Rcpp::export]]
Eigen::SparseMatrix<double> make_standardized_matern_eigen(int dim_x, int dim_y, double rho1, double rho2, int nu) {
    int N = dim_x * dim_y;
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

    // Compute marginal standard deviations
    Eigen::ArrayXd marginal_sds = marginal_sd_eigen(A1, V1, dim_x, A2, V2, dim_y, nu);

    // Create diagonal matrix D
    Eigen::DiagonalMatrix<double, Eigen::Dynamic> D(N);
    D.diagonal() = marginal_sds;

    // Create full Q matrix using Kronecker sum
    Eigen::SparseMatrix<double> I1(dim_x, dim_x);
    I1.setIdentity();

    Eigen::SparseMatrix<double> I2(dim_y, dim_y);
    I2.setIdentity();


    Eigen::SparseMatrix<double> Q = Eigen::kroneckerProduct(Q1, I2) + Eigen::kroneckerProduct(I1, Q2);

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
Rcpp::NumericVector dmatern_copula_eigen(const Eigen::MatrixXd& X, int dim_x, int dim_y, double rho1, double rho2, int nu) {
    int n_obs = X.cols();
    int N = dim_x * dim_y;
    ArrayXd quadform_sums = ArrayXd::Zero(n_obs);
    Rcpp::Timer timer;

    // Create precision matrices
    MatrixXd Q1 = make_AR_prec_matrix(dim_x, rho1);
    MatrixXd Q2 = make_AR_prec_matrix(dim_y, rho2);

    // Perform eigendecompositions
    SelfAdjointEigenSolver<MatrixXd> solver1(Q1);
    VectorXd A1 = solver1.eigenvalues();
    MatrixXd V1 = solver1.eigenvectors();

    SelfAdjointEigenSolver<MatrixXd> solver2(Q2);
    VectorXd A2 = solver2.eigenvalues();
    MatrixXd V2 = solver2.eigenvectors();

    timer.step("Eigendecomposition");

    // Compute marginal standard deviations
    ArrayXd marginal_sds = ArrayXd::Zero(dim_x * dim_y);
    
    #pragma omp parallel
    {
        ArrayXd local_marginal_sds = ArrayXd::Zero(dim_x * dim_y);
        VectorXd v(dim_x * dim_y);

        #pragma omp for collapse(2) nowait
        for (int i = 0; i < dim_x; ++i) {
            for (int j = 0; j < dim_y; ++j) {
                v = kroneckerProduct(V1.col(i), V2.col(j));
                local_marginal_sds += v.array().square() / std::pow(A1(i) + A2(j), nu + 1);
            }
        }

        #pragma omp critical
        {
            marginal_sds += local_marginal_sds;
        }
    }

    timer.step("Marginal SDs");

    double log_det = 0;

    // Transpose X for better memory access if it's stored in column-major order
    MatrixXd X_t = X.transpose();

    #pragma omp parallel reduction(+:log_det)
    {
        ArrayXd local_quadform_sums = ArrayXd::Zero(n_obs);
        VectorXd v = VectorXd::Zero(dim_x * dim_y);
        VectorXd u = VectorXd::Zero(dim_x * dim_y);

        #pragma omp for collapse(2) nowait
        for (int i = 0; i < dim_x; ++i) {
            for (int j = 0; j < dim_y; ++j) {
                v = kroneckerProduct(V1.col(i), V2.col(j)).array() * marginal_sds;
                double norm_v = v.norm();
                v /= norm_v;
                
                double A = std::pow(A1(i) + A2(j), nu + 1);
                double lambda = A * (norm_v * norm_v);

                log_det += std::log(lambda);
                
                u = X_t * v;
                local_quadform_sums += (u.array().square() * lambda);
            }
        }

        #pragma omp critical
        {
            quadform_sums += local_quadform_sums;
        }
    }

    timer.step("Quadform Sums");

    ArrayXd x_squared = X.colwise().squaredNorm().array();
    ArrayXd log_densities = -0.5 * (quadform_sums - log_det - x_squared);
    timer.step("Log Densities");

    Rcpp::NumericVector res(timer);
    

    return res;
}

// [[Rcpp::export]]
Eigen::MatrixXd rmatern_copula_eigen(int n, int dim_x, int dim_y, double rho1, double rho2, int nu) {
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

    // Step 3: Compute marginal standard deviations
    Eigen::VectorXd marginal_sds = marginal_sd_eigen(A1, V1, dim_x, A2, V2, dim_y, nu);

    // Step 4: Prepare for sampling
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
                double lambda = std::pow(A1(i) + A2(j), -(nu + 1.0) / 2.0);
                Eigen::VectorXd v = Eigen::kroneckerProduct(V1.col(i), V2.col(j));
                x += lambda * d(generators[thread_id]) * v;
            }
        }

        // Standardize the sample
        x = x.array() / marginal_sds.array();

        samples.col(s) = x;
    }

    return samples;
}