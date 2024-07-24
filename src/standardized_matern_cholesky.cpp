#include <RcppEigen.h>
#include <omp.h>
#include <random>
#include <queue>
#include <unordered_set>


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

double compute_sigma_ii(const Eigen::MatrixXd& L, int i) {
    int n = L.rows();
    Eigen::VectorXd ei = Eigen::VectorXd::Unit(n, i);
    
    // Solve Ly = e_i
    Eigen::VectorXd y = L.triangularView<Eigen::Lower>().solve(ei);
    
    // Solve L^T x = y
    Eigen::VectorXd x = L.triangularView<Eigen::Lower>().adjoint().solve(y);
    
    // Return x_i
    return x(i);
}

// Helper function to compute I(i)
std::vector<std::set<int>> compute_I(const Eigen::SparseMatrix<double>& L) {
    int n = L.rows();
    std::vector<std::set<int>> I(n);
    
    for (int k = 0; k < L.outerSize(); ++k) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(L, k); it; ++it) {
            if (it.row() > it.col()) {
                I[it.col()].insert(it.row());
            }
        }
    }
    
    return I;
}

// Rue's algorithm for computing marginal variances
Eigen::VectorXd compute_marginal_variances_rue(const Eigen::SparseMatrix<double>& Q) {
    int n = Q.rows();
    Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> llt(Q);
    Eigen::SparseMatrix<double> L = llt.matrixL();
    
    std::vector<std::set<int>> I = compute_I(L);
    
    Eigen::MatrixXd Sigma(n, n);
    
    for (int i = n - 1; i >= 0; --i) {
        for (int j = n - 1; j >= i; --j) {
            if (i == j) {
                Sigma(i, i) = 1.0 / (L.coeff(i, i) * L.coeff(i, i));
            } else {
                double sum = 0.0;
                for (int k : I[i]) {
                    sum += L.coeff(k, i) * Sigma(k, j);
                }
                Sigma(i, j) = Sigma(j, i) = -sum / L.coeff(i, i);
            }
        }
    }
    
    return Sigma.diagonal();
}

// [[Rcpp::export]]
Eigen::SparseMatrix<double> compute_normalized_cholesky(const Eigen::SparseMatrix<double>& Q) {
    Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> llt(Q);
    int n = Q.rows();
    Eigen::SparseMatrix<double> L = llt.matrixL();
    
    // Compute marginal variances using Rue's algorithm
    Eigen::VectorXd marginal_vars = compute_marginal_variances_rue(L);


    // Create diagonal matrix D
    Eigen::DiagonalMatrix<double, Eigen::Dynamic> D = marginal_vars.array().sqrt().matrix().asDiagonal();
    
    // Compute standardized Cholesky factor
    L = D * L;
    L.makeCompressed();
    return L;
}

// [[Rcpp::export]]
Eigen::VectorXd matern_mvn_density_cholesky(const Eigen::MatrixXd& X, int dim, double rho, int nu) {
    int N = X.rows();  // number of elements in each observation
    int n_obs = X.cols();  // number of observations
    const double C = N * std::log(2.0 * M_PI);
    
    // Create the Matérn precision matrix
    Eigen::SparseMatrix<double> Q = make_matern_prec_matrix(dim, rho, nu);
    
    // Compute normalized Cholesky factor
    Eigen::SparseMatrix<double> L = compute_normalized_cholesky(Q);

    // Compute log determinant
    double log_det = 2 * L.diagonal().array().log().sum();
    
    // Compute all quadratic forms at once
    Eigen::MatrixXd Z = L.transpose() * X;
    Eigen::VectorXd quadforms = Z.colwise().squaredNorm();
    
    // Compute log densities for all observations
    Eigen::VectorXd log_densities = -0.5 * (C - log_det + quadforms.array());
    
    return log_densities;
}