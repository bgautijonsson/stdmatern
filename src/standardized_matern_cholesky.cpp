#include <RcppEigen.h>
#include <omp.h>
#include <random>


// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::plugins(openmp)]]

using namespace Rcpp;
using namespace Eigen;


// Function to create a 1-dimensional AR(1) precision matrix
Eigen::SparseMatrix<double> make_AR_prec_matrix_chol(int dim, double rho) {
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

// Function to create the 2-dimensional Matérn precision matrix using Kronecker products
// [[Rcpp::export]]
Eigen::SparseMatrix<double> make_matern_prec_matrix(int dim, double rho, int nu) {
    Eigen::SparseMatrix<double> Q = make_AR_prec_matrix_chol(dim, rho);
    Eigen::SparseMatrix<double> I = Eigen::SparseMatrix<double>(dim, dim);
    I.setIdentity();

    // Compute Kronecker products
    Eigen::SparseMatrix<double> QI = Eigen::kroneckerProduct(Q, I);
    Eigen::SparseMatrix<double> IQ = Eigen::kroneckerProduct(I, Q);

    // Sum the Kronecker products
    Eigen::SparseMatrix<double> result = QI + IQ;
    result.makeCompressed();

    // Apply matrix multiplication nu times
    if (nu > 0) {
        Eigen::SparseMatrix<double> temp = result;
        for (int i = 0; i < nu; ++i) {
            result = result * temp;
            result.makeCompressed();
        }
    }

    return result;
}

// Function to compute marginal variances
// [[Rcpp::export]]
Eigen::SparseMatrix<double> marginal_sd_cholesky(const Eigen::SparseMatrix<double>& Q) {
    Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> llt(Q);
    int n = Q.rows();
    Eigen::SparseMatrix<double> D(Q.rows(), Q.cols());
    D.reserve(Eigen::VectorXi::Constant(n, 1));

    for (int i = 0; i < n; ++i) {
        Eigen::VectorXd ei = Eigen::VectorXd::Unit(n, i);
        double marginal_variance = llt.solve(ei)(i);
        D.insert(i, i) = std::sqrt(marginal_variance);
    }

    D.makeCompressed();
    return D;
}

// Function to create the standardized Matérn precision matrix
// [[Rcpp::export]]
Eigen::SparseMatrix<double> make_standardized_matern_cholesky(int dim, double rho, int nu) {
    // Create the Matérn precision matrix
    Eigen::SparseMatrix<double> Q = make_matern_prec_matrix(dim, rho, nu);

    // Compute marginal variances
    Eigen::SparseMatrix<double> D = marginal_sd_cholesky(Q);

    Eigen::SparseMatrix<double> Q_standardized = D * Q * D;

    return Q_standardized;
}

std::pair<Eigen::SparseMatrix<double>, double> compute_normalized_cholesky(const Eigen::SparseMatrix<double>& Q) {
    Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> llt(Q);
    int n = Q.rows();
    Eigen::SparseMatrix<double> L = llt.matrixL();
    
    double log_det = 0.0;
    Eigen::VectorXd D_inv(n);

    #pragma omp parallel
    {
        #pragma omp for reduction(+:log_det)
        for (int i = 0; i < n; ++i) {
            Eigen::VectorXd ei = Eigen::VectorXd::Unit(n, i);
            double marginal_variance = llt.solve(ei)(i);
            double marginal_sd = std::sqrt(marginal_variance);
            D_inv(i) = 1.0 / marginal_sd;
            log_det += std::log(L.coeff(i, i));
        }

        // Normalize L
        #pragma omp for
        for (int k = 0; k < L.outerSize(); ++k) {
            for (Eigen::SparseMatrix<double>::InnerIterator it(L, k); it; ++it) {
                it.valueRef() *= D_inv(it.row());
            }
        }
    }

    L.makeCompressed();
    return std::make_pair(L, 2.0 * log_det);
}

// [[Rcpp::export]]
Eigen::VectorXd matern_mvn_density_cholesky(const Eigen::MatrixXd& X, int dim, double rho, int nu) {
    int n = X.rows();  // number of elements in each observation
    int n_obs = X.cols();  // number of observations
    const double log_2pi = std::log(2.0 * M_PI);
    
    // Create the Matérn precision matrix
    Eigen::SparseMatrix<double> Q = make_matern_prec_matrix(dim, rho, nu);
    
    // Compute normalized Cholesky factor and log determinant
    auto [L, log_det] = compute_normalized_cholesky(Q);
    
    // Prepare output vector
    Eigen::VectorXd log_densities(n_obs);
    
    // Compute log densities for each observation
    #pragma omp parallel for
    for (int i = 0; i < n_obs; ++i) {
        // Compute quadratic form (x^T * Q * x)
        Eigen::VectorXd z = L.triangularView<Eigen::Lower>().solve(X.col(i));
        double quadform = z.squaredNorm();
        
        // Compute log density
        log_densities(i) = -0.5 * (n * log_2pi + log_det + quadform);
    }
    
    return log_densities;
}