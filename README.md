
<!-- README.md is generated from README.Rmd. Please edit that file -->

# cholesky

This directory is for development of fast and memory-efficient code that
creates Matérn precision matrices that have been standardized so that
their inverse is a correlation matrix. The code is written in C++ and
made available inside R with the `{Rcpp}` packages.

The package can be installed with

``` r
pak::pak("bgautijonsson/cholesky")
```

``` r
library(cholesky)
#> Loading required package: Matrix
```

``` r
Q <- make_standardized_matern(dim = 2, rho = 0.5)
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
  make_standardized_matern(dim = 40, rho = 0.5)
)
#> # A tibble: 1 × 6
#>   expression                             min median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr>                          <bch:> <bch:>     <dbl> <bch:byt>    <dbl>
#> 1 make_standardized_matern(dim = 40,… 55.5ms 58.3ms      15.9    98.3KB        0
```
