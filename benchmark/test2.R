library(stdmatern)
library(sparseMVN)
library(purrr)
library(tidyverse)

eval_fun <- function(rho) {
    
    Q1 <- make_AR_prec_matrix(grid_dim, rho)
    eig <- eigen(Q1)
    A1 <- eig$values
    V1 <- eig$vectors

    V <- kronecker(V1, V1)

    A <- kronecker(diag(A1), diag(rep(1, grid_dim))) + 
        kronecker(diag(rep(1, grid_dim)), diag(A1))
    
    A <- A^(nu + 1)

    
    msd <- fast_marginal_standard_deviations(A1, V1, grid_dim, nu)
    
    V <- diag(msd) %*% V
    A <- (diag(colSums(V^2)) * A)
    V <- diag(1/sqrt(colSums(V^2))) %*% V
    
    log_dens <- numeric(n_replicate)

    logdet <- 0

    for (i in seq_len(grid_dim)) {
        for (j in seq_len(grid_dim)) {
            lambda <- (A1[i] + A1[j])
            logdet <- logdet + log(lambda)
        }
    }

    logdet <- logdet * (nu + 1)
    logdet <- diag(A) |> log() |> sum()


    Y <- t(V) %*% Z

    quadform <- colSums(A %*% Y^2)

    
    # for (obs in seq_len(n_replicate)) {
    #     quadform <- 0
        
        
    #     for (i in seq_len(grid_dim)) {
    #         for (j in seq_len(grid_dim)) {
    #             v <- kronecker(V1[ , j], V1[ , i])
    #             lambda <- (A1[i] + A1[j])^(nu + 1)
                
    #             u <- sum(v * (Z[ , obs] * msd))
    #             quadform <- quadform + u^2 * lambda
    #         }
    #     }
    #     log_dens[obs] <-  - 0.5 * (n_replicate * log(2 * pi) - logdet + quadform)
    # }
    
    log_dens <- - 0.5 * (n_replicate * log(2 * pi) - logdet + quadform)

    log_dens |> sum()
}

fun2 <- function(rho) {
    Q <- make_standardized_matern_eigen(grid_dim, rho, nu)
    
    L <- chol(Q)
    q <- L %*% Z
    quadform <- colSums(q^2)
    detsum <- L |> diag() |> log() |> sum()
    C <- log(2 * pi)
    
    sum(- n_replicate / 2 * C + detsum - quadform / 2)
}

neg_ll_function <- function(params) {
    rho <- params[1]
    # - sum(matern_mvn_density(Z, grid_dim, rho, nu))
    
    #-fun2(rho)
    -eval_fun(rho)
}

grid_dim <- 10
rho <- 0.9
n_replicate <- 10
nu <- 2


Z <- sample_standardized_matern(grid_dim, rho, nu, n_replicate)
res <- optimize(
    f = neg_ll_function,
    interval = c(0, 1)
)

res



grid_dim <- 10
rho <- 0.5
n_replicate <- 10
nu <- 0


Z <- sample_standardized_matern(grid_dim, rho, nu, n_replicate)

Z <- matrix(
    rnorm(grid_dim^2 * n_replicate),
    nrow = grid_dim^2,
    ncol = n_replicate
)

Z <- matrix(
    0,
    nrow = grid_dim^2,
    ncol = n_replicate
)

Q <- make_standardized_matern_eigen(grid_dim, rho, nu)

Z <- rmvn.sparse(
    n_replicate,
    mu = rep(0, grid_dim^2),
    CH = Cholesky(Q)
) |> t()

map(
    seq(0.1, 0.9, length.out = 100),
    \(rho) {
        tibble(
            rho = rho,
            chol = -fun2(rho),
            #eigen = -sum(matern_mvn_density(Z, grid_dim, rho, nu))
            eigen = -eval_fun(rho)
        )
    }
) |> 
list_rbind() |> 
mutate(
    min_chol = rho[chol == min(chol)],
    min_eig = rho[eigen == min(eigen)]
) |> 
ggplot(aes(rho, chol)) +
geom_line(aes(col = "Cholesky")) +
geom_line(aes(y = eigen, col = "Eigen"))+
geom_vline(aes(xintercept = min_chol, col = "Cholesky"), lty = 2) +
geom_vline(aes(xintercept = min_eig, col = "Eigen"), lty = 2) 
scale_y_log10()






get_eigs <- function(dim, rho) {
    Q1 <- make_AR_prec_matrix(dim, rho)
    eig <- eigen(Q1)
    A1 <- eig$values
    V1 <- eig$vectors

    A <- c()

    for (i in seq_len(dim)) {
        for (j in seq_len(dim)) {
            A <- c(A, (A1[i] + A1[j])^(nu + 1))
        }
    }

    A
}


get_eigs(10, 0.99)

eigen(make_matern_prec_matrix(3, 0.9, 0))$values

dim <- 30
map(
    seq(0.01, 0.99, length.out = 10),
    \(rho) {
        tibble(
            rho = rho,
            eig = get_eigs(dim, rho) |> log() |> sum(),
            det = make_matern_prec_matrix(dim, rho, 0) |> 
                chol() |> 
                diag() |> 
                log() |> 
                sum()
        )
    }
) |> 
    list_rbind() |> 
    ggplot(aes(eig, 2 * det)) +
    geom_point() +
    geom_abline(intercept = 0, slope = 1) +
    scale_x_continuous(expand = expansion(0.2)) +
    scale_y_continuous(expand = expansion(0.2))


get_eigs(10, 0.5) |> log() |> sum()

make_matern_prec_matrix(10, 0.5, 0) |> det() |> log()


make_matern_prec_matrix(15, 0.95, 0) |> chol() |> diag() |> log() |> sum()



