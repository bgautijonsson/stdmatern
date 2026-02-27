#include <RcppEigen.h>
#include <omp.h>
#include <random>
#include "ar_matrix.h"
#include <nlopt.hpp>



// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::plugins(openmp)]]
// [[Rcpp::depends(nloptr)]]

using namespace Rcpp;
using namespace Eigen;


// Function to compute the marginal standard deviations of a Matérn-like field
// [[Rcpp::export]]
Eigen::ArrayXd marginal_sd_eigen(
    const Eigen::VectorXd& A1, const Eigen::MatrixXd& V1, int dim_x, 
    const Eigen::VectorXd& A2, const Eigen::MatrixXd& V2, int dim_y, 
    int nu
) {
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

    return marginal_sds.sqrt();
}

// Function to create the standardized Matérn precision matrix using eigendecomposition
// [[Rcpp::export]]
Eigen::SparseMatrix<double> make_standardized_matern_eigen(int dim_x, int dim_y, double rho1, double rho2, int nu) {
    int N = dim_x * dim_y;
    // Create precision matrices for each dimension
    SparseMatrix<double> Q1 = make_AR_prec_matrix(dim_x, rho1);
    SparseMatrix<double> Q2 = make_AR_prec_matrix(dim_y, rho2);

    // Perform eigendecompositions of Q1 and Q2
    SelfAdjointEigenSolver<SparseMatrix<double>> solver1(Q1);
    VectorXd A1 = solver1.eigenvalues();
    MatrixXd V1 = solver1.eigenvectors();

    SelfAdjointEigenSolver<SparseMatrix<double>> solver2(Q2);
    VectorXd A2 = solver2.eigenvalues();
    MatrixXd V2 = solver2.eigenvectors();

    // Compute marginal standard deviations
    ArrayXd marginal_sds = marginal_sd_eigen(A1, V1, dim_x, A2, V2, dim_y, nu);

    // Create diagonal matrix D with marginal standard deviations
    DiagonalMatrix<double, Dynamic> D(N);
    D.diagonal() = marginal_sds;

    // Create identity matrices for Kronecker product
    SparseMatrix<double> I1(dim_x, dim_x);
    I1.setIdentity();

    SparseMatrix<double> I2(dim_y, dim_y);
    I2.setIdentity();

    // Create full Q matrix using Kronecker sum
    SparseMatrix<double> Q = kroneckerProduct(Q1, I2) + kroneckerProduct(I1, Q2);

    // Apply matrix multiplication nu times for higher order Matérn fields
    if (nu > 0) {
        SparseMatrix<double> temp = Q;
        for (int i = 0; i < nu; ++i) {
            Q = Q * temp;
            Q.makeCompressed();
        }
    }

    // Standardize Q to ensure its inverse is a correlation matrix
    SparseMatrix<double> Q_standardized = D * Q * D;

    return Q_standardized;
}

// Function to compute the log-density of a Matérn copula
// [[Rcpp::export]]
Eigen::ArrayXd dmatern_copula_eigen(const Eigen::MatrixXd& X, int dim_x, int dim_y, double rho1, double rho2, int nu) {
    int n_obs = X.cols();
    int N = X.rows();
    
    // Create precision matrices for each dimension
    MatrixXd Q1 = make_AR_prec_matrix(dim_x, rho1);
    MatrixXd Q2 = make_AR_prec_matrix(dim_y, rho2);

    // Perform eigendecompositions of Q1 and Q2
    SelfAdjointEigenSolver<MatrixXd> solver1(Q1);
    VectorXd A1 = solver1.eigenvalues();
    MatrixXd V1 = solver1.eigenvectors();

    SelfAdjointEigenSolver<MatrixXd> solver2(Q2);
    VectorXd A2 = solver2.eigenvalues();
    MatrixXd V2 = solver2.eigenvectors();
    
    // Compute marginal standard deviations in parallel
    ArrayXd marginal_sds = marginal_sd_eigen(A1, V1, dim_x, A2, V2, dim_y, nu);

    MatrixXd X_scaled = (X.array().colwise() * marginal_sds).matrix();
    
    // Since Q_tilde = D * Q * D, we have log(det(Q_tilde)) = log(det(D)) + log(det(Q)) + log(det(D))
    // Since D is diagonal, we have log(det(D)) = sum(log(diag(D))) = sum(log(marginal_sds))
    double log_det_D = marginal_sds.log().sum();
    double log_det_Q = 0.0;

    // Compute quadratic forms in parallel
    ArrayXd quadform_sums = ArrayXd::Zero(n_obs);
    
    // The log-density for a normal distribution with precision matrix Q is
    // -0.5 * (x^T Q x + log(det(Q)) + log(2 * pi))
    #pragma omp parallel reduction(+:log_det_Q)
    {
        ArrayXd local_quadform_sums = ArrayXd::Zero(n_obs);
        VectorXd v(dim_x * dim_y);
        
        #pragma omp for collapse(2) nowait
        for (int i = 0; i < dim_x; ++i) {
            for (int j = 0; j < dim_y; ++j) {
                v = kroneckerProduct(V1.col(i), V2.col(j));
                
                double lambda = std::pow(A1(i) + A2(j), nu + 1.0);
                log_det_Q += std::log(lambda);
                
                // Compute for all observations
                for (int k = 0; k < n_obs; ++k) {
                    double dot_product = v.dot(X_scaled.col(k));
                    local_quadform_sums(k) += dot_product * dot_product * lambda;
                }
            }
        }

        #pragma omp critical
        {
            quadform_sums += local_quadform_sums;
        }
    }

    // Compute log determinant of Q_tilde = DQD as log(det(Q_tilde)) = log(det(D)) + log(det(Q)) + log(det(D))
    double log_det = log_det_Q + 2 * log_det_D;

    // Compute log densities
    ArrayXd x_squared = X.colwise().squaredNorm().array();
    ArrayXd log_densities = -0.5 * (quadform_sums - log_det - x_squared);
    
    return log_densities;
}

// [[Rcpp::export]]
Eigen::MatrixXd rmatern_copula_eigen(int n, int dim_x, int dim_y, double rho1, double rho2, int nu) {
    // Step 1: Create 1D AR(1) precision matrices
    SparseMatrix<double> Q1 = make_AR_prec_matrix(dim_x, rho1);
    SparseMatrix<double> Q2 = make_AR_prec_matrix(dim_y, rho2);
    
    // Step 2: Perform eigendecomposition of Q1 and Q2
    SelfAdjointEigenSolver<SparseMatrix<double>> solver1(Q1);
    VectorXd A1 = solver1.eigenvalues();
    MatrixXd V1 = solver1.eigenvectors();

    SelfAdjointEigenSolver<SparseMatrix<double>> solver2(Q2);
    VectorXd A2 = solver2.eigenvalues();
    MatrixXd V2 = solver2.eigenvectors();

    // Step 3: Compute marginal standard deviations
    ArrayXd marginal_sds = marginal_sd_eigen(A1, V1, dim_x, A2, V2, dim_y, nu);

    // Step 4: Prepare for sampling
    int D = dim_x * dim_y;
    MatrixXd samples = MatrixXd::Zero(D, n);
    
    // Random number generation with per-thread generators
    std::vector<std::mt19937> generators(omp_get_max_threads());
    for (auto& gen : generators) {
        std::random_device rd;
        gen.seed(rd());
    }
    std::normal_distribution<> d(0, 1);

    // Step 5: Generate samples
    #pragma omp parallel
    {
        int thread_id = omp_get_thread_num();
        auto& gen = generators[thread_id];
        MatrixXd local_samples = MatrixXd::Zero(D, n);
        VectorXd v(D);

        #pragma omp for collapse(2) nowait
        for (int i = 0; i < dim_x; ++i) {
            for (int j = 0; j < dim_y; ++j) {
                double lambda = std::pow(A1(i) + A2(j), -(nu + 1.0) / 2.0);
                v = kroneckerProduct(V1.col(i), V2.col(j));
                
                for (int s = 0; s < n; ++s) {
                    local_samples.col(s) += lambda * d(gen) * v;
                }
            }
        }

        #pragma omp critical
        {
            samples += local_samples;
        }
    }

    // Step 6: Apply marginal standard deviations
    samples = samples.array().colwise() / marginal_sds;

    return samples;
}


// ============================================================================
// MLE Estimation for Matérn Copula
// ============================================================================

// Struct to pass data to the NLopt objective function
struct OptData {
    const MatrixXd* Y;
    int dim_x;
    int dim_y;
    int nu;
};

// Objective function for NLopt: computes negative log-likelihood
double nll_objective(const std::vector<double>& x, std::vector<double>& grad, void* data) {
    OptData* opt_data = static_cast<OptData*>(data);
    
    double rho1 = x[0];
    double rho2 = x[1];
    
    // Compute log-densities using existing function
    ArrayXd log_densities = dmatern_copula_eigen(
        *(opt_data->Y),
        opt_data->dim_x,
        opt_data->dim_y,
        rho1,
        rho2,
        opt_data->nu
    );
    
    // Return negative sum of log-densities (we minimize, so negate for MLE)
    double nll = -log_densities.sum();
    if (!std::isfinite(nll)) return 1e15;
    return nll;
}

// Function to fit the Matérn copula parameters using MLE
// [[Rcpp::export]]
Rcpp::List fit_matern_copula_eigen(
    const Eigen::MatrixXd& Y,
    int dim_x, int dim_y, int nu,
    Rcpp::Nullable<Eigen::VectorXd> start = R_NilValue
) {
    // Set up optimization data
    OptData opt_data;
    opt_data.Y = &Y;
    opt_data.dim_x = dim_x;
    opt_data.dim_y = dim_y;
    opt_data.nu = nu;
    
    // Initialize starting values
    std::vector<double> x(2);
    if (start.isNotNull()) {
        VectorXd start_vec = Rcpp::as<VectorXd>(start);
        x[0] = start_vec(0);
        x[1] = start_vec(1);
    } else {
        // Default starting values
        x[0] = 0.5;
        x[1] = 0.5;
    }
    
    // Set up NLopt optimizer with BOBYQA algorithm (derivative-free, bounded)
    nlopt::opt optimizer(nlopt::LN_BOBYQA, 2);
    
    // Set bounds: rho1 and rho2 in [0.001, 0.999]
    std::vector<double> lb = {0.001, 0.001};
    std::vector<double> ub = {0.999, 0.999};
    optimizer.set_lower_bounds(lb);
    optimizer.set_upper_bounds(ub);
    
    // Set objective function
    optimizer.set_min_objective(nll_objective, &opt_data);
    
    // Set convergence tolerances
    optimizer.set_xtol_rel(1e-6);
    optimizer.set_ftol_rel(1e-8);
    optimizer.set_maxeval(1000);
    optimizer.set_maxtime(120.0);

    // Run optimization
    double min_nll;
    nlopt::result result;
    std::string message;

    try {
        result = optimizer.optimize(x, min_nll);

        // Convert result code to message
        switch (result) {
            case nlopt::SUCCESS:
                message = "Optimization succeeded";
                break;
            case nlopt::STOPVAL_REACHED:
                message = "Stopval reached";
                break;
            case nlopt::FTOL_REACHED:
                message = "Function tolerance reached";
                break;
            case nlopt::XTOL_REACHED:
                message = "Parameter tolerance reached";
                break;
            case nlopt::MAXEVAL_REACHED:
                message = "Maximum evaluations reached";
                break;
            case nlopt::MAXTIME_REACHED:
                message = "Maximum time reached";
                break;
            default:
                message = "Unknown result";
        }
    } catch (std::exception& e) {
        return Rcpp::List::create(
            Rcpp::Named("rho1") = NA_REAL,
            Rcpp::Named("rho2") = NA_REAL,
            Rcpp::Named("loglik") = NA_REAL,
            Rcpp::Named("convergence") = -1,
            Rcpp::Named("message") = std::string("Optimization failed: ") + e.what()
        );
    }
    
    return Rcpp::List::create(
        Rcpp::Named("rho1") = x[0],
        Rcpp::Named("rho2") = x[1],
        Rcpp::Named("loglik") = -min_nll,
        Rcpp::Named("convergence") = static_cast<int>(result),
        Rcpp::Named("message") = message
    );
}

// Objective function for isotropic case (single rho for both directions)
double nll_objective_isotropic(const std::vector<double>& x, std::vector<double>& grad, void* data) {
    OptData* opt_data = static_cast<OptData*>(data);
    
    double rho = x[0];  // Single rho for both directions
    
    // Compute log-densities using existing function with rho1 = rho2 = rho
    ArrayXd log_densities = dmatern_copula_eigen(
        *(opt_data->Y),
        opt_data->dim_x,
        opt_data->dim_y,
        rho,
        rho,
        opt_data->nu
    );
    
    // Return negative sum of log-densities
    double nll = -log_densities.sum();
    if (!std::isfinite(nll)) return 1e15;
    return nll;
}

// Function to fit the Matérn copula with isotropic correlation (single rho)
// [[Rcpp::export]]
Rcpp::List fit_matern_copula_eigen_isotropic(
    const Eigen::MatrixXd& Y,
    int dim_x, int dim_y, int nu,
    Rcpp::Nullable<double> start = R_NilValue
) {
    // Set up optimization data
    OptData opt_data;
    opt_data.Y = &Y;
    opt_data.dim_x = dim_x;
    opt_data.dim_y = dim_y;
    opt_data.nu = nu;
    
    // Initialize starting value
    std::vector<double> x(1);
    if (start.isNotNull()) {
        x[0] = Rcpp::as<double>(start);
    } else {
        x[0] = 0.5;
    }
    
    // Set up NLopt optimizer with BOBYQA algorithm (1D optimization)
    nlopt::opt optimizer(nlopt::LN_BOBYQA, 1);
    
    // Set bounds: rho in [0.001, 0.999]
    std::vector<double> lb = {0.001};
    std::vector<double> ub = {0.999};
    optimizer.set_lower_bounds(lb);
    optimizer.set_upper_bounds(ub);
    
    // Set objective function
    optimizer.set_min_objective(nll_objective_isotropic, &opt_data);
    
    // Set convergence tolerances
    optimizer.set_xtol_rel(1e-6);
    optimizer.set_ftol_rel(1e-8);
    optimizer.set_maxeval(1000);
    optimizer.set_maxtime(120.0);

    // Run optimization
    double min_nll;
    nlopt::result result;
    std::string message;

    try {
        result = optimizer.optimize(x, min_nll);

        // Convert result code to message
        switch (result) {
            case nlopt::SUCCESS:
                message = "Optimization succeeded";
                break;
            case nlopt::STOPVAL_REACHED:
                message = "Stopval reached";
                break;
            case nlopt::FTOL_REACHED:
                message = "Function tolerance reached";
                break;
            case nlopt::XTOL_REACHED:
                message = "Parameter tolerance reached";
                break;
            case nlopt::MAXEVAL_REACHED:
                message = "Maximum evaluations reached";
                break;
            case nlopt::MAXTIME_REACHED:
                message = "Maximum time reached";
                break;
            default:
                message = "Unknown result";
        }
    } catch (std::exception& e) {
        return Rcpp::List::create(
            Rcpp::Named("rho") = NA_REAL,
            Rcpp::Named("loglik") = NA_REAL,
            Rcpp::Named("convergence") = -1,
            Rcpp::Named("message") = std::string("Optimization failed: ") + e.what()
        );
    }
    
    return Rcpp::List::create(
        Rcpp::Named("rho") = x[0],
        Rcpp::Named("loglik") = -min_nll,
        Rcpp::Named("convergence") = static_cast<int>(result),
        Rcpp::Named("message") = message
    );
}
