library(stdmatern)
library(dplyr)

dim <- 5
rho <- 0.5
nu <- 2

Q1 <- make_standardized_matern(dim, rho, nu) 

Q2 <- make_standardized_matern_eigen(dim, rho, nu)

Q1 |> solve() |> diag()
Q2 |> solve() |> diag()

mean(abs(Q1[Q1 != 0] - Q2[Q2 != 0]))

fun1 <- function(dim, rho, nu) {
  out <-  make_standardized_matern(dim, rho, nu) 
  return(1)
}

fun2 <- function(dim, rho, nu) {
  out <-  make_standardized_matern_eigen(dim, rho, nu) 
  return(1)
}

dim <- 30
rho <- 0.5
nu <- 2

bench::mark(
  "Cholesky" = fun1(dim, rho, nu),
  "Eigen" = fun2(dim, rho, nu)
) |> 
  mutate(
    relative = as.numeric(median / min(median)),
    .before = min
  )
