#ifndef AR_MATRIX_H
#define AR_MATRIX_H

#include <RcppEigen.h>

Eigen::SparseMatrix<double> make_AR_prec_matrix(int dim, double rho);

#endif