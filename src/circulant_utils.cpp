#include "circulant_utils.h"
#include <RcppEigen.h>
#include <complex>
#include <random>
#include <fftw3.h>

// [[Rcpp::depends(RcppEigen, fftwtools)]]

using namespace Rcpp;
using namespace Eigen;

Eigen::MatrixXcd fft2(const Eigen::MatrixXcd& input, double scale) {
    int dim = input.rows();
    Eigen::MatrixXcd output(dim, dim);

    fftw_complex *in, *out;
    in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * dim * dim);
    out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * dim * dim);

    // Copy data from Eigen matrix to FFTW input
    for (int i = 0; i < dim; ++i) {
        for (int j = 0; j < dim; ++j) {
            in[i*dim + j][0] = input(i, j).real();
            in[i*dim + j][1] = input(i, j).imag();
        }
    }

    // Forward FFT
    fftw_plan plan = fftw_plan_dft_2d(dim, dim, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(plan);

    // Copy result back to Eigen matrix and scale by 1/(dim*dim)
    for (int i = 0; i < dim; ++i) {
        for (int j = 0; j < dim; ++j) {
            output(i, j) = std::complex<double>(out[i*dim + j][0], out[i*dim + j][1]) * scale;
        }
    }

    fftw_destroy_plan(plan);
    fftw_free(in);
    fftw_free(out);

    return output;
}

// Function to create the base matrix c
// [[Rcpp::export]]
Eigen::MatrixXd create_base_matrix(int dim, double rho) {
    Eigen::MatrixXd c = Eigen::MatrixXd::Zero(dim, dim);
    double scale = 1.0 / (1.0 + std::pow(rho, 2));
    
    // Set the first row
    c(0, 0) = (2 + 2 * rho * rho) * scale;
    c(0, 1) = -rho * scale;
    c(0, dim - 1) = -rho * scale;
    
    // Set the second and last row
    c(1, 0) = -rho * scale;
    c(dim - 1, 0) = -rho * scale;
    
    return c;
}

// Function to compute FFT of c, apply nu, and rescale eigenvalues using marginal variance
// [[Rcpp::export]]
Eigen::MatrixXcd compute_and_rescale_eigenvalues(const Eigen::MatrixXd& c, int nu) {
    int dim = c.rows();

    fftw_complex *in, *out;
    in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * dim * dim);
    out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * dim * dim);

    // Copy data from Eigen matrix to FFTW input
    for (int i = 0; i < dim; ++i) {
        for (int j = 0; j < dim; ++j) {
            in[i*dim + j][0] = c(i, j);
            in[i*dim + j][1] = 0.0;
        }
    }

    // Forward FFT
    fftw_plan plan_forward = fftw_plan_dft_2d(dim, dim, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(plan_forward);

    // Copy eigenvalues to Eigen matrix and apply nu
    Eigen::MatrixXcd eigenvalues(dim, dim);
    for (int i = 0; i < dim; ++i) {
        for (int j = 0; j < dim; ++j) {
            std::complex<double> eig(out[i*dim + j][0], out[i*dim + j][1]);
            eigenvalues(i, j) = std::pow(eig, nu + 1);
        }
    }

    // Compute inverse of eigenvalues
    Eigen::MatrixXcd inv_eigenvalues = eigenvalues.cwiseInverse();

    // Inverse FFT of inverse eigenvalues
    for (int i = 0; i < dim * dim; ++i) {
        in[i][0] = inv_eigenvalues(i / dim, i % dim).real();
        in[i][1] = inv_eigenvalues(i / dim, i % dim).imag();
    }

    fftw_plan plan_backward = fftw_plan_dft_2d(dim, dim, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_execute(plan_backward);

    // Get marginal variance (first element of inverse FFT result)
    double mvar = out[0][0] / (dim * dim);

    // Scale eigenvalues
    Eigen::MatrixXcd updated_eigenvalues = mvar * eigenvalues;

    fftw_destroy_plan(plan_forward);
    fftw_destroy_plan(plan_backward);
    fftw_free(in);
    fftw_free(out);

    return updated_eigenvalues;
}

// Function for matrix-vector product using FFT and eigenvalues
// [[Rcpp::export]]
Eigen::VectorXd matrix_vector_product(const Eigen::MatrixXcd& eigenvalues, const Eigen::VectorXd& v) {
    int dim = eigenvalues.rows();

    fftw_complex *in, *out;
    in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * dim * dim);
    out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * dim * dim);

    // Copy data from Eigen vector to FFTW input
    for (int i = 0; i < dim * dim; ++i) {
        in[i][0] = v(i);
        in[i][1] = 0.0;
    }

    // Forward FFT
    fftw_plan plan_forward = fftw_plan_dft_2d(dim, dim, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(plan_forward);

    // Multiply with eigenvalues in frequency domain
    for (int i = 0; i < dim * dim; ++i) {
        std::complex<double> temp(out[i][0], out[i][1]);
        temp *= eigenvalues(i / dim, i % dim);
        out[i][0] = temp.real();
        out[i][1] = temp.imag();
    }

    // Inverse FFT
    fftw_plan plan_backward = fftw_plan_dft_2d(dim, dim, out, in, FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_execute(plan_backward);

    // Copy result back to Eigen vector and normalize
    Eigen::VectorXd result(dim * dim);
    for (int i = 0; i < dim * dim; ++i) {
        result(i) = in[i][0] / (dim * dim);
    }

    fftw_destroy_plan(plan_forward);
    fftw_destroy_plan(plan_backward);
    fftw_free(in);
    fftw_free(out);

    return result;
}