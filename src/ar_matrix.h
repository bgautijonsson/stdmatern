#ifndef AR_MATRIX_H
#define AR_MATRIX_H

#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]

Eigen::SparseMatrix<double> make_AR_prec_matrix_cpp(int dim, double rho);

#endif
