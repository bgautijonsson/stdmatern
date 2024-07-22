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
    if (nu > 0) {
        Eigen::SparseMatrix<double> temp = result;
        for (int i = 0; i < nu; ++i) {
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





// [[Rcpp::export]]
Eigen::VectorXd fast_marginal_standard_deviations(const Eigen::VectorXd& A1, const Eigen::MatrixXd& V1, int dim, int nu) {
    Eigen::VectorXd marginal_sds = Eigen::VectorXd::Zero(dim * dim);
    
    for (int i = 0; i < dim; ++i) {
        for (int j = 0; j < dim; ++j) {
            Eigen::VectorXd v = Eigen::kroneckerProduct(V1.col(i), V1.col(j));
            double lambda = std::pow(A1(i) + A1(j), nu + 1);
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
    Eigen::VectorXd marginal_sds = fast_marginal_standard_deviations(A1, V1, dim, nu);

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
    int D = dim * dim;
    Eigen::VectorXd log_densities(n_obs);
    const double C = D * std::log(2 * M_PI);

    // Create precision matrix
    Eigen::MatrixXd Q1 = make_AR_prec_matrix(dim, rho);

    // Perform eigendecomposition
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(Q1);
    Eigen::VectorXd A1 = solver.eigenvalues();
    Eigen::MatrixXd V1 = solver.eigenvectors();

    // Compute marginal standard deviations
    Eigen::VectorXd marginal_sds = fast_marginal_standard_deviations(A1, V1, dim, nu);


    
    for (int obs = 0; obs < n_obs; ++obs) {
        double quadform_sum = 0;
        double log_det_sum = 0;
        
        // Scale X_slice
        Eigen::VectorXd X_slice = X.col(obs).array() * marginal_sds.array();
        
        for (int i = 0; i < dim; ++i) {
            for (int j = 0; j < dim; ++j) {
                Eigen::VectorXd v = Eigen::kroneckerProduct(V1.col(j), V1.col(i));
                double lambda = std::pow(A1(i) + A1(j), nu + 1);
                log_det_sum -= std::log(lambda);
                
                double u = v.dot(X_slice);
                quadform_sum += u * u * lambda;
            }
        }
        
        log_densities(obs) = -0.5 * (C + log_det_sum + quadform_sum);
    }

    return log_densities;
}

// Function to generate samples from a standardized Matérn field
// [[Rcpp::export]]
Eigen::MatrixXd sample_standardized_matern(int dim, double rho, int nu, int n_samples) {
    // Step 1: Create 1D AR(1) precision matrix
    Eigen::SparseMatrix<double> Q1 = make_AR_prec_matrix(dim, rho);
    
    // Step 2: Perform eigendecomposition of Q1
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(Q1);
    Eigen::VectorXd A1 = solver.eigenvalues();
    Eigen::MatrixXd V1 = solver.eigenvectors();

    // Step 3: Compute marginal standard deviations
    Eigen::VectorXd marginal_sds = fast_marginal_standard_deviations(A1, V1, dim, nu);

    // Step 4: Prepare for sampling
    int D = dim * dim;
    Eigen::MatrixXd samples(D, n_samples);
    
    // Random number generation
    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<> d(0, 1);

    // Step 5: Generate samples
    for (int s = 0; s < n_samples; ++s) {

        Eigen::VectorXd x = Eigen::VectorXd::Zero(D);

        for (int i = 0; i < dim; ++i) {
            for (int j = 0; j < dim; ++j) {
                Eigen::VectorXd v = Eigen::kroneckerProduct(V1.col(j), V1.col(i));
                double lambda = std::pow(A1(i) + A1(j), -(nu + 1.0) / 2.0);
                x += lambda * d(gen) * v;
            }
        }

        // Standardize the sample
        x = x.array() / marginal_sds.array();

        samples.col(s) = x;
    }

    return samples;
}