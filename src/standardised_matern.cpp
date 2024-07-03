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
Eigen::VectorXd dmvn_chol_cpp(const Eigen::MatrixXd& X, const Eigen::SparseMatrix<double>& L) {
    int n_obs = X.cols();
    int n = X.rows();
    double C = -0.918938533204672669541 * n;
    double log_det = L.diagonal().array().log().sum();
    
    Eigen::VectorXd log_densities(n_obs);
    
    #pragma omp parallel for
    for (int i = 0; i < n_obs; ++i) {
        Eigen::VectorXd q = L * X.col(i);
        double quadform = q.squaredNorm();
        log_densities(i) = C + log_det - quadform/2;
    }
    
    return log_densities;
}


// [[Rcpp::export]]
Eigen::VectorXd matern_mvn_density_cpp(const Eigen::MatrixXd& X, int dim, double rho, int nu) {
    
    // Create the standardized Matérn Cholesky matrix
    Eigen::SparseMatrix<double> L = make_standardized_matern_cholesky(dim, rho, nu);
    
    // Calculate and return the densities
    return dmvn_chol_cpp(X, L);
}

// [[Rcpp::export]]
Eigen::VectorXd fast_marginal_variances(const Eigen::SparseMatrix<double>& Q1) {
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eig(Q1);
    Eigen::VectorXd A1 = eig.eigenvalues();
    Eigen::MatrixXd V1 = eig.eigenvectors();

    int dim = Q1.cols();
    
    Eigen::VectorXd msd = Eigen::VectorXd::Zero(dim * dim);
    
    for (int i = 0; i < dim; ++i) {
        for (int j = 0; j < dim; ++j) {
            Eigen::VectorXd V_temp = Eigen::kroneckerProduct(V1.col(j), V1.col(i));
            double A_temp = A1(i) + A1(j);
            msd += (V_temp.array().square() / A_temp).matrix();
        }
    }
    
    msd = msd.array().sqrt();
    
    return msd;
}

// Simplified function to compute log-density using eigendecomposition (μ = 0)
// [[Rcpp::export]]
double log_density_eigen(const Eigen::VectorXd& x,
                         double eigenvalue, const Eigen::VectorXd& eigenvector) {
    int n = x.size();
    double transformed = eigenvector.dot(x);
    double quadform = transformed * transformed * eigenvalue;
    double logdet = std::log(eigenvalue);
    
    return -0.5 * (n * std::log(2 * M_PI) + logdet + quadform);
}


// [[Rcpp::export]]
Eigen::VectorXd matern_mvn_density_eigen(const Eigen::MatrixXd& X, int dim, double rho, int nu) {
    Eigen::SparseMatrix<double> Q1 = make_AR_prec_matrix(dim, rho);
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eig(Q1);
    Eigen::VectorXd A1 = eig.eigenvalues();
    Eigen::MatrixXd V1 = eig.eigenvectors();
    
    int n_obs = X.cols();
    Eigen::VectorXd log_densities(n_obs);
    
    // Compute marginal standard deviations
    Eigen::VectorXd marginal_sds = fast_marginal_variances(Q1);
    
    double log_det_sum = A1.array().log().sum() * dim;  // Sum of log of eigenvalues of Q
    
    #pragma omp parallel for
    for (int obs = 0; obs < n_obs; ++obs) {
        double quadform_sum = 0;
        Eigen::VectorXd X_slice = X.col(obs).array() * marginal_sds.array();
        
        for (int i = 0; i < dim; ++i) {
            for (int j = 0; j < dim; ++j) {
                Eigen::VectorXd v = Eigen::kroneckerProduct(V1.col(j), V1.col(i));
                double lambda = A1(i) + A1(j);
                
                double transformed = v.dot(X_slice);
                quadform_sum += transformed * transformed * lambda;
            }
        }
        
        int n = dim * dim;
        log_densities(obs) = -0.5 * (n * std::log(2 * M_PI) + log_det_sum + quadform_sum);
    }
    
    return log_densities;
}