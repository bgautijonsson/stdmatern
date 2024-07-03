#include <RcppEigen.h>
#include <omp.h>
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

// Function to create the 2-dimensional Matérn precision matrix using Kronecker products
// [[Rcpp::export]]
Eigen::SparseMatrix<double> make_matern_prec_matrix(int dim, double rho, int nu) {
    Eigen::SparseMatrix<double> Q = make_AR_prec_matrix(dim, rho);
    Eigen::SparseMatrix<double> I = Eigen::SparseMatrix<double>(dim, dim);
    I.setIdentity();

    // Compute Kronecker products
    Eigen::SparseMatrix<double> QI = Eigen::kroneckerProduct(Q, I);
    Eigen::SparseMatrix<double> IQ = Eigen::kroneckerProduct(I, Q);

    // Sum the Kronecker products
    Eigen::SparseMatrix<double> result = QI + IQ;
    result.makeCompressed();

    // Apply matrix multiplication nu times
    if (nu > 1) {
        Eigen::SparseMatrix<double> temp = result;
        for (int i = 1; i < nu; ++i) {
            temp = temp * result;
            temp.makeCompressed();
        }
        result = temp;
    }

    return result;
}

// Function to compute marginal variances
// [[Rcpp::export]]
Eigen::SparseMatrix<double> compute_marginal_variances(const Eigen::SparseMatrix<double>& Q) {
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
Eigen::SparseMatrix<double> make_standardized_matern(int dim, double rho, int nu) {
    // Create the Matérn precision matrix
    Eigen::SparseMatrix<double> Q = make_matern_prec_matrix(dim, rho, nu);

    // Compute marginal variances
    Eigen::SparseMatrix<double> D = compute_marginal_variances(Q);

    Eigen::SparseMatrix<double> Q_standardized = D * Q * D;

    return Q_standardized;
}


// Function to create the Cholesky factor of the standardized Matérn precision matrix
// [[Rcpp::export]]
Eigen::SparseMatrix<double> make_standardized_matern_cholesky(int dim, double rho, int nu) {
    Eigen::SparseMatrix<double> Q = make_matern_prec_matrix(dim, rho, nu);
    Eigen::SimplicialLLT<Eigen::SparseMatrix<double>, Eigen::Lower, Eigen::NaturalOrdering<int>> llt(Q);
    Eigen::SparseMatrix<double> L = llt.matrixU();
    Eigen::VectorXd marginal_std(Q.rows());
    

    // Compute marginal standard deviations
    #pragma omp parallel for
    for (int i = 0; i < Q.rows(); ++i) {
        Eigen::VectorXd ei = Eigen::VectorXd::Unit(Q.rows(), i);
        marginal_std(i) = std::sqrt(llt.solve(ei)(i));
    }
    
    // Standardize L
    #pragma omp parallel for
    for (int k = 0; k < Q.rows(); ++k) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(L, k); it; ++it) {
            it.valueRef() *= marginal_std(it.row());
        }
    }
    
    // Convert to lower triangular sparse matrix
    return L;
}

// [[Rcpp::export]]
double dmvn_chol(const Eigen::VectorXd& x, const Eigen::SparseMatrix<double>& L) {
    
    int n = x.size();
    double C = -0.918938533204672669541 * n;
    
    Eigen::VectorXd q = L * x;
    
    // Compute log-density
    double quadform = q.squaredNorm();
    double log_det = L.diagonal().array().log().sum();
    
    return C + log_det - quadform/2;
}


// [[Rcpp::export]]
double matern_mvn_density(const Eigen::VectorXd& x, int dim, double rho, int nu) {
    
    // Create the standardized Matérn Cholesky matrix
    Eigen::SparseMatrix<double> L = make_standardized_matern_cholesky(dim, rho, nu);
    
    // Calculate and return the density
    return dmvn_chol(x, L);
}