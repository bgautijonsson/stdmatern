# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

make_AR_prec_matrix <- function(dim, rho) {
    .Call(`_stdmatern_make_AR_prec_matrix`, dim, rho)
}

create_base_matrix <- function(dim1, dim2, rho1, rho2) {
    .Call(`_stdmatern_create_base_matrix`, dim1, dim2, rho1, rho2)
}

compute_and_rescale_eigenvalues <- function(c, nu) {
    .Call(`_stdmatern_compute_and_rescale_eigenvalues`, c, nu)
}

matrix_vector_product <- function(eigenvalues, v) {
    .Call(`_stdmatern_matrix_vector_product`, eigenvalues, v)
}

make_matern_prec_matrix <- function(dim_x, dim_y, rho1, rho2, nu) {
    .Call(`_stdmatern_make_matern_prec_matrix`, dim_x, dim_y, rho1, rho2, nu)
}

dmatern_cholesky <- function(X, dim_x, dim_y, rho1, rho2, nu) {
    .Call(`_stdmatern_dmatern_cholesky`, X, dim_x, dim_y, rho1, rho2, nu)
}

marginal_sd_cholesky <- function(L, nu) {
    .Call(`_stdmatern_marginal_sd_cholesky`, L, nu)
}

dmatern_copula_cholesky <- function(X, dim_x, dim_y, rho1, rho2, nu) {
    .Call(`_stdmatern_dmatern_copula_cholesky`, X, dim_x, dim_y, rho1, rho2, nu)
}

dmatern_copula_circulant <- function(X, dim1, dim2, rho1, rho2, nu) {
    .Call(`_stdmatern_dmatern_copula_circulant`, X, dim1, dim2, rho1, rho2, nu)
}

rmatern_copula_circulant <- function(n_samples, dim1, dim2, rho1, rho2, nu) {
    .Call(`_stdmatern_rmatern_copula_circulant`, n_samples, dim1, dim2, rho1, rho2, nu)
}

construct_circulant_precision <- function(dim1, dim2, rho1, rho2, nu) {
    .Call(`_stdmatern_construct_circulant_precision`, dim1, dim2, rho1, rho2, nu)
}

marginal_sd_eigen <- function(A1, V1, dim_x, A2, V2, dim_y, nu) {
    .Call(`_stdmatern_marginal_sd_eigen`, A1, V1, dim_x, A2, V2, dim_y, nu)
}

make_standardized_matern_eigen <- function(dim_x, dim_y, rho1, rho2, nu) {
    .Call(`_stdmatern_make_standardized_matern_eigen`, dim_x, dim_y, rho1, rho2, nu)
}

dmatern_copula_eigen <- function(X, dim_x, dim_y, rho1, rho2, nu) {
    .Call(`_stdmatern_dmatern_copula_eigen`, X, dim_x, dim_y, rho1, rho2, nu)
}

rmatern_copula_eigen <- function(n, dim_x, dim_y, rho1, rho2, nu) {
    .Call(`_stdmatern_rmatern_copula_eigen`, n, dim_x, dim_y, rho1, rho2, nu)
}

fold_data <- function(X, dim1, dim2) {
    .Call(`_stdmatern_fold_data`, X, dim1, dim2)
}

dmatern_copula_folded <- function(X, dim1, dim2, rho1, rho2, nu) {
    .Call(`_stdmatern_dmatern_copula_folded`, X, dim1, dim2, rho1, rho2, nu)
}

rmatern_copula_folded_full <- function(n_samples, dim1, dim2, rho1, rho2, nu) {
    .Call(`_stdmatern_rmatern_copula_folded_full`, n_samples, dim1, dim2, rho1, rho2, nu)
}

dmatern_eigen <- function(X, dim_x, dim_y, rho1, rho2, nu) {
    .Call(`_stdmatern_dmatern_eigen`, X, dim_x, dim_y, rho1, rho2, nu)
}

rmatern_eigen <- function(n, dim_x, dim_y, rho1, rho2, nu) {
    .Call(`_stdmatern_rmatern_eigen`, n, dim_x, dim_y, rho1, rho2, nu)
}

