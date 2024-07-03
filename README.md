
<!-- README.md is generated from README.Rmd. Please edit that file -->

# stdmatern

This directory is for development of fast and memory-efficient code that
creates Matérn precision matrices that have been standardized so that
their inverse is a correlation matrix. The code is written in C++ and
made available inside R with the `{Rcpp}` packages.

The package can be installed with

``` r
pak::pak("bgautijonsson/stdmatern")
```

``` r
library(stdmatern)
#> Loading required package: Matrix
```

``` r
Q <- make_standardized_matern(dim = 2, rho = 0.5, nu = 0)
```

``` r
Q
#> 4 x 4 sparse Matrix of class "dgCMatrix"
#>                                                 
#> [1,]  1.1666667 -0.2916667 -0.2916667  .        
#> [2,] -0.2916667  1.1666667  .         -0.2916667
#> [3,] -0.2916667  .          1.1666667 -0.2916667
#> [4,]  .         -0.2916667 -0.2916667  1.1666667
```

``` r
Q |> solve()
#> 4 x 4 sparse Matrix of class "dgCMatrix"
#>                                             
#> [1,] 1.0000000 0.2857143 0.2857143 0.1428571
#> [2,] 0.2857143 1.0000000 0.1428571 0.2857143
#> [3,] 0.2857143 0.1428571 1.0000000 0.2857143
#> [4,] 0.1428571 0.2857143 0.2857143 1.0000000
```

Creating and standardizing a 1600x1600 precision matrix

``` r
bench::mark(
  make_standardized_matern(dim = 40, rho = 0.5, nu = 0)
)
#> # A tibble: 1 × 6
#>   expression                             min median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr>                          <bch:> <bch:>     <dbl> <bch:byt>    <dbl>
#> 1 make_standardized_matern(dim = 40,… 55.1ms 55.5ms      18.0    98.3KB        0
```

# Sampling spatial data

``` r
start <- tictoc::tic()
grid_dim <- 80
rho <- 0.9
nu <- 2
Q <- make_standardized_matern(grid_dim, rho, nu = nu)


Z <- rmvn.sparse(
  n = 1,
  mu = rep(0, nrow(Q)),
  CH = Matrix::Cholesky(Q)
)

tibble(
  Z = as.numeric(Z)
) |> 
  mutate(
    id = row_number(),
    lat = (id - 1) %% grid_dim,
    lon = cumsum(lat == 0),
  ) |> 
  ggplot(aes(lat, lon, fill = Z)) +
  geom_raster() +
  scale_fill_viridis_c() +
  coord_fixed(expand = FALSE)
```

<img src="man/figures/README-unnamed-chunk-8-1.png" width="100%" />

``` r
stop <- tictoc::toc()
#> 3.385 sec elapsed
```

# Normal density

The package also implements a method for calculating the log-density of
a multivariate normal by creating and scaling an appropriate precision
matrix. This function is meant for use inside MCMC samplers or
optimization algorithms. If you plan to use the same precision matrix
often, it’s better to create and store the precision matrix instead of
calculating it again.

``` r
grid_dim <- 50
rho <- 0.5
nu <- 0
x <- rnorm(grid_dim^2)
bench::mark(
  matern_mvn_density(x, grid_dim, rho, nu)
)
#> # A tibble: 1 × 6
#>   expression                             min median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr>                          <bch:> <bch:>     <dbl> <bch:byt>    <dbl>
#> 1 matern_mvn_density(x, grid_dim, rh… 5.66ms 5.71ms      175.    40.6KB        0
```

The function can also take in a matrix, with ncol = n_replicates and
nrow = grim_dim^2

``` r
grid_dim <- 50
n_replicates <- 10
rho <- 0.5
nu <- 0
Q <- make_standardized_matern(grid_dim, rho, nu = nu)
L <- make_standardized_matern_cholesky(grid_dim, rho, nu = nu)
X <- rmvn.sparse(
  n = n_replicates,
  mu = rep(0, nrow(Q)),
  CH = Matrix::Cholesky(Q)
) |> t()

X |> dim()
#> [1] 2500   10
```

``` r

bench::mark(
  matern_mvn_density(X, grid_dim, rho, nu)
)
#> # A tibble: 1 × 6
#>   expression                             min median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr>                          <bch:> <bch:>     <dbl> <bch:byt>    <dbl>
#> 1 matern_mvn_density(X, grid_dim, rh… 6.91ms 7.93ms      125.        0B        0
```

``` r
library(dplyr)
library(sparseMVN)
library(ggplot2)
library(stdmatern)

grid_dim <- 100
n_replicates <- 1
rho <- 0.5
nu <- 0
Q <- make_standardized_matern(grid_dim, rho, nu = nu)
L <- make_standardized_matern_cholesky(grid_dim, rho, nu = nu)
X <- rmvn.sparse(
  n = n_replicates,
  mu = rep(0, nrow(Q)),
  CH = Matrix::Cholesky(Q)
) |> t()

bench::mark(
  matern_mvn_density(X, grid_dim, rho, nu)
)
#> # A tibble: 1 × 6
#>   expression                             min median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr>                          <bch:> <bch:>     <dbl> <bch:byt>    <dbl>
#> 1 matern_mvn_density(X, grid_dim, rh… 96.7ms 99.7ms      9.44    78.2KB        0
```
