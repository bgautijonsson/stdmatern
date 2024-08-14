// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// make_AR_prec_matrix
Eigen::SparseMatrix<double> make_AR_prec_matrix(int dim, double rho);
RcppExport SEXP _stdmatern_make_AR_prec_matrix(SEXP dimSEXP, SEXP rhoSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type dim(dimSEXP);
    Rcpp::traits::input_parameter< double >::type rho(rhoSEXP);
    rcpp_result_gen = Rcpp::wrap(make_AR_prec_matrix(dim, rho));
    return rcpp_result_gen;
END_RCPP
}
// create_base_matrix
Eigen::MatrixXd create_base_matrix(int dim1, int dim2, double rho1, double rho2);
RcppExport SEXP _stdmatern_create_base_matrix(SEXP dim1SEXP, SEXP dim2SEXP, SEXP rho1SEXP, SEXP rho2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type dim1(dim1SEXP);
    Rcpp::traits::input_parameter< int >::type dim2(dim2SEXP);
    Rcpp::traits::input_parameter< double >::type rho1(rho1SEXP);
    Rcpp::traits::input_parameter< double >::type rho2(rho2SEXP);
    rcpp_result_gen = Rcpp::wrap(create_base_matrix(dim1, dim2, rho1, rho2));
    return rcpp_result_gen;
END_RCPP
}
// compute_and_rescale_eigenvalues
Eigen::MatrixXcd compute_and_rescale_eigenvalues(const Eigen::MatrixXd& c, int nu);
RcppExport SEXP _stdmatern_compute_and_rescale_eigenvalues(SEXP cSEXP, SEXP nuSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type c(cSEXP);
    Rcpp::traits::input_parameter< int >::type nu(nuSEXP);
    rcpp_result_gen = Rcpp::wrap(compute_and_rescale_eigenvalues(c, nu));
    return rcpp_result_gen;
END_RCPP
}
// matrix_vector_product
Eigen::VectorXd matrix_vector_product(const Eigen::MatrixXcd& eigenvalues, const Eigen::VectorXd& v);
RcppExport SEXP _stdmatern_matrix_vector_product(SEXP eigenvaluesSEXP, SEXP vSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::MatrixXcd& >::type eigenvalues(eigenvaluesSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type v(vSEXP);
    rcpp_result_gen = Rcpp::wrap(matrix_vector_product(eigenvalues, v));
    return rcpp_result_gen;
END_RCPP
}
// dmatern_copula_circulant
Eigen::VectorXd dmatern_copula_circulant(const Eigen::MatrixXd& X, int dim1, int dim2, double rho1, double rho2, int nu);
RcppExport SEXP _stdmatern_dmatern_copula_circulant(SEXP XSEXP, SEXP dim1SEXP, SEXP dim2SEXP, SEXP rho1SEXP, SEXP rho2SEXP, SEXP nuSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type X(XSEXP);
    Rcpp::traits::input_parameter< int >::type dim1(dim1SEXP);
    Rcpp::traits::input_parameter< int >::type dim2(dim2SEXP);
    Rcpp::traits::input_parameter< double >::type rho1(rho1SEXP);
    Rcpp::traits::input_parameter< double >::type rho2(rho2SEXP);
    Rcpp::traits::input_parameter< int >::type nu(nuSEXP);
    rcpp_result_gen = Rcpp::wrap(dmatern_copula_circulant(X, dim1, dim2, rho1, rho2, nu));
    return rcpp_result_gen;
END_RCPP
}
// rmatern_copula_circulant
Eigen::MatrixXd rmatern_copula_circulant(int n_samples, int dim1, int dim2, double rho1, double rho2, int nu);
RcppExport SEXP _stdmatern_rmatern_copula_circulant(SEXP n_samplesSEXP, SEXP dim1SEXP, SEXP dim2SEXP, SEXP rho1SEXP, SEXP rho2SEXP, SEXP nuSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n_samples(n_samplesSEXP);
    Rcpp::traits::input_parameter< int >::type dim1(dim1SEXP);
    Rcpp::traits::input_parameter< int >::type dim2(dim2SEXP);
    Rcpp::traits::input_parameter< double >::type rho1(rho1SEXP);
    Rcpp::traits::input_parameter< double >::type rho2(rho2SEXP);
    Rcpp::traits::input_parameter< int >::type nu(nuSEXP);
    rcpp_result_gen = Rcpp::wrap(rmatern_copula_circulant(n_samples, dim1, dim2, rho1, rho2, nu));
    return rcpp_result_gen;
END_RCPP
}
// construct_circulant_precision
Eigen::SparseMatrix<double> construct_circulant_precision(int dim1, int dim2, double rho1, double rho2, int nu);
RcppExport SEXP _stdmatern_construct_circulant_precision(SEXP dim1SEXP, SEXP dim2SEXP, SEXP rho1SEXP, SEXP rho2SEXP, SEXP nuSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type dim1(dim1SEXP);
    Rcpp::traits::input_parameter< int >::type dim2(dim2SEXP);
    Rcpp::traits::input_parameter< double >::type rho1(rho1SEXP);
    Rcpp::traits::input_parameter< double >::type rho2(rho2SEXP);
    Rcpp::traits::input_parameter< int >::type nu(nuSEXP);
    rcpp_result_gen = Rcpp::wrap(construct_circulant_precision(dim1, dim2, rho1, rho2, nu));
    return rcpp_result_gen;
END_RCPP
}
// marginal_sd_eigen
Eigen::VectorXd marginal_sd_eigen(const Eigen::VectorXd& A1, const Eigen::MatrixXd& V1, int dim_x, const Eigen::VectorXd& A2, const Eigen::MatrixXd& V2, int dim_y, int nu);
RcppExport SEXP _stdmatern_marginal_sd_eigen(SEXP A1SEXP, SEXP V1SEXP, SEXP dim_xSEXP, SEXP A2SEXP, SEXP V2SEXP, SEXP dim_ySEXP, SEXP nuSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type A1(A1SEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type V1(V1SEXP);
    Rcpp::traits::input_parameter< int >::type dim_x(dim_xSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type A2(A2SEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type V2(V2SEXP);
    Rcpp::traits::input_parameter< int >::type dim_y(dim_ySEXP);
    Rcpp::traits::input_parameter< int >::type nu(nuSEXP);
    rcpp_result_gen = Rcpp::wrap(marginal_sd_eigen(A1, V1, dim_x, A2, V2, dim_y, nu));
    return rcpp_result_gen;
END_RCPP
}
// make_standardized_matern_eigen
Eigen::SparseMatrix<double> make_standardized_matern_eigen(int dim_x, int dim_y, double rho1, double rho2, int nu);
RcppExport SEXP _stdmatern_make_standardized_matern_eigen(SEXP dim_xSEXP, SEXP dim_ySEXP, SEXP rho1SEXP, SEXP rho2SEXP, SEXP nuSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type dim_x(dim_xSEXP);
    Rcpp::traits::input_parameter< int >::type dim_y(dim_ySEXP);
    Rcpp::traits::input_parameter< double >::type rho1(rho1SEXP);
    Rcpp::traits::input_parameter< double >::type rho2(rho2SEXP);
    Rcpp::traits::input_parameter< int >::type nu(nuSEXP);
    rcpp_result_gen = Rcpp::wrap(make_standardized_matern_eigen(dim_x, dim_y, rho1, rho2, nu));
    return rcpp_result_gen;
END_RCPP
}
// dmatern_copula_eigen
Eigen::VectorXd dmatern_copula_eigen(const Eigen::MatrixXd& X, int dim_x, int dim_y, double rho1, double rho2, int nu);
RcppExport SEXP _stdmatern_dmatern_copula_eigen(SEXP XSEXP, SEXP dim_xSEXP, SEXP dim_ySEXP, SEXP rho1SEXP, SEXP rho2SEXP, SEXP nuSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type X(XSEXP);
    Rcpp::traits::input_parameter< int >::type dim_x(dim_xSEXP);
    Rcpp::traits::input_parameter< int >::type dim_y(dim_ySEXP);
    Rcpp::traits::input_parameter< double >::type rho1(rho1SEXP);
    Rcpp::traits::input_parameter< double >::type rho2(rho2SEXP);
    Rcpp::traits::input_parameter< int >::type nu(nuSEXP);
    rcpp_result_gen = Rcpp::wrap(dmatern_copula_eigen(X, dim_x, dim_y, rho1, rho2, nu));
    return rcpp_result_gen;
END_RCPP
}
// rmatern_copula_eigen
Eigen::MatrixXd rmatern_copula_eigen(int n, int dim_x, int dim_y, double rho1, double rho2, int nu);
RcppExport SEXP _stdmatern_rmatern_copula_eigen(SEXP nSEXP, SEXP dim_xSEXP, SEXP dim_ySEXP, SEXP rho1SEXP, SEXP rho2SEXP, SEXP nuSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< int >::type dim_x(dim_xSEXP);
    Rcpp::traits::input_parameter< int >::type dim_y(dim_ySEXP);
    Rcpp::traits::input_parameter< double >::type rho1(rho1SEXP);
    Rcpp::traits::input_parameter< double >::type rho2(rho2SEXP);
    Rcpp::traits::input_parameter< int >::type nu(nuSEXP);
    rcpp_result_gen = Rcpp::wrap(rmatern_copula_eigen(n, dim_x, dim_y, rho1, rho2, nu));
    return rcpp_result_gen;
END_RCPP
}
// fold_data
Eigen::VectorXd fold_data(const Eigen::VectorXd& X, int n1, int n2);
RcppExport SEXP _stdmatern_fold_data(SEXP XSEXP, SEXP n1SEXP, SEXP n2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type X(XSEXP);
    Rcpp::traits::input_parameter< int >::type n1(n1SEXP);
    Rcpp::traits::input_parameter< int >::type n2(n2SEXP);
    rcpp_result_gen = Rcpp::wrap(fold_data(X, n1, n2));
    return rcpp_result_gen;
END_RCPP
}
// dmatern_copula_folded
Eigen::VectorXd dmatern_copula_folded(const Eigen::MatrixXd& X, int dim1, int dim2, double rho1, double rho2, int nu);
RcppExport SEXP _stdmatern_dmatern_copula_folded(SEXP XSEXP, SEXP dim1SEXP, SEXP dim2SEXP, SEXP rho1SEXP, SEXP rho2SEXP, SEXP nuSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type X(XSEXP);
    Rcpp::traits::input_parameter< int >::type dim1(dim1SEXP);
    Rcpp::traits::input_parameter< int >::type dim2(dim2SEXP);
    Rcpp::traits::input_parameter< double >::type rho1(rho1SEXP);
    Rcpp::traits::input_parameter< double >::type rho2(rho2SEXP);
    Rcpp::traits::input_parameter< int >::type nu(nuSEXP);
    rcpp_result_gen = Rcpp::wrap(dmatern_copula_folded(X, dim1, dim2, rho1, rho2, nu));
    return rcpp_result_gen;
END_RCPP
}
// rmatern_copula_folded_full
Eigen::MatrixXd rmatern_copula_folded_full(int n_samples, int dim1, int dim2, double rho1, double rho2, int nu);
RcppExport SEXP _stdmatern_rmatern_copula_folded_full(SEXP n_samplesSEXP, SEXP dim1SEXP, SEXP dim2SEXP, SEXP rho1SEXP, SEXP rho2SEXP, SEXP nuSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n_samples(n_samplesSEXP);
    Rcpp::traits::input_parameter< int >::type dim1(dim1SEXP);
    Rcpp::traits::input_parameter< int >::type dim2(dim2SEXP);
    Rcpp::traits::input_parameter< double >::type rho1(rho1SEXP);
    Rcpp::traits::input_parameter< double >::type rho2(rho2SEXP);
    Rcpp::traits::input_parameter< int >::type nu(nuSEXP);
    rcpp_result_gen = Rcpp::wrap(rmatern_copula_folded_full(n_samples, dim1, dim2, rho1, rho2, nu));
    return rcpp_result_gen;
END_RCPP
}
// make_matern_prec_matrix
Eigen::SparseMatrix<double> make_matern_prec_matrix(int dim_x, int dim_y, double rho1, double rho2, int nu);
RcppExport SEXP _stdmatern_make_matern_prec_matrix(SEXP dim_xSEXP, SEXP dim_ySEXP, SEXP rho1SEXP, SEXP rho2SEXP, SEXP nuSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type dim_x(dim_xSEXP);
    Rcpp::traits::input_parameter< int >::type dim_y(dim_ySEXP);
    Rcpp::traits::input_parameter< double >::type rho1(rho1SEXP);
    Rcpp::traits::input_parameter< double >::type rho2(rho2SEXP);
    Rcpp::traits::input_parameter< int >::type nu(nuSEXP);
    rcpp_result_gen = Rcpp::wrap(make_matern_prec_matrix(dim_x, dim_y, rho1, rho2, nu));
    return rcpp_result_gen;
END_RCPP
}
// dmatern_eigen
Eigen::VectorXd dmatern_eigen(const Eigen::MatrixXd& X, int dim_x, int dim_y, double rho1, double rho2, int nu);
RcppExport SEXP _stdmatern_dmatern_eigen(SEXP XSEXP, SEXP dim_xSEXP, SEXP dim_ySEXP, SEXP rho1SEXP, SEXP rho2SEXP, SEXP nuSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type X(XSEXP);
    Rcpp::traits::input_parameter< int >::type dim_x(dim_xSEXP);
    Rcpp::traits::input_parameter< int >::type dim_y(dim_ySEXP);
    Rcpp::traits::input_parameter< double >::type rho1(rho1SEXP);
    Rcpp::traits::input_parameter< double >::type rho2(rho2SEXP);
    Rcpp::traits::input_parameter< int >::type nu(nuSEXP);
    rcpp_result_gen = Rcpp::wrap(dmatern_eigen(X, dim_x, dim_y, rho1, rho2, nu));
    return rcpp_result_gen;
END_RCPP
}
// rmatern_eigen
Eigen::MatrixXd rmatern_eigen(int n, int dim_x, int dim_y, double rho1, double rho2, int nu);
RcppExport SEXP _stdmatern_rmatern_eigen(SEXP nSEXP, SEXP dim_xSEXP, SEXP dim_ySEXP, SEXP rho1SEXP, SEXP rho2SEXP, SEXP nuSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< int >::type dim_x(dim_xSEXP);
    Rcpp::traits::input_parameter< int >::type dim_y(dim_ySEXP);
    Rcpp::traits::input_parameter< double >::type rho1(rho1SEXP);
    Rcpp::traits::input_parameter< double >::type rho2(rho2SEXP);
    Rcpp::traits::input_parameter< int >::type nu(nuSEXP);
    rcpp_result_gen = Rcpp::wrap(rmatern_eigen(n, dim_x, dim_y, rho1, rho2, nu));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_stdmatern_make_AR_prec_matrix", (DL_FUNC) &_stdmatern_make_AR_prec_matrix, 2},
    {"_stdmatern_create_base_matrix", (DL_FUNC) &_stdmatern_create_base_matrix, 4},
    {"_stdmatern_compute_and_rescale_eigenvalues", (DL_FUNC) &_stdmatern_compute_and_rescale_eigenvalues, 2},
    {"_stdmatern_matrix_vector_product", (DL_FUNC) &_stdmatern_matrix_vector_product, 2},
    {"_stdmatern_dmatern_copula_circulant", (DL_FUNC) &_stdmatern_dmatern_copula_circulant, 6},
    {"_stdmatern_rmatern_copula_circulant", (DL_FUNC) &_stdmatern_rmatern_copula_circulant, 6},
    {"_stdmatern_construct_circulant_precision", (DL_FUNC) &_stdmatern_construct_circulant_precision, 5},
    {"_stdmatern_marginal_sd_eigen", (DL_FUNC) &_stdmatern_marginal_sd_eigen, 7},
    {"_stdmatern_make_standardized_matern_eigen", (DL_FUNC) &_stdmatern_make_standardized_matern_eigen, 5},
    {"_stdmatern_dmatern_copula_eigen", (DL_FUNC) &_stdmatern_dmatern_copula_eigen, 6},
    {"_stdmatern_rmatern_copula_eigen", (DL_FUNC) &_stdmatern_rmatern_copula_eigen, 6},
    {"_stdmatern_fold_data", (DL_FUNC) &_stdmatern_fold_data, 3},
    {"_stdmatern_dmatern_copula_folded", (DL_FUNC) &_stdmatern_dmatern_copula_folded, 6},
    {"_stdmatern_rmatern_copula_folded_full", (DL_FUNC) &_stdmatern_rmatern_copula_folded_full, 6},
    {"_stdmatern_make_matern_prec_matrix", (DL_FUNC) &_stdmatern_make_matern_prec_matrix, 5},
    {"_stdmatern_dmatern_eigen", (DL_FUNC) &_stdmatern_dmatern_eigen, 6},
    {"_stdmatern_rmatern_eigen", (DL_FUNC) &_stdmatern_rmatern_eigen, 6},
    {NULL, NULL, 0}
};

RcppExport void R_init_stdmatern(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
