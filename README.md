
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
Q <- make_standardized_matern(dim = 3, rho = 0.5, nu = 0)
```

``` r
Q
#> 9 x 9 sparse Matrix of class "dgCMatrix"
#>                                                                                
#>  [1,]  1.155093 -0.2808520  .        -0.2808520  .          .          .       
#>  [2,] -0.280852  1.2291667 -0.280852  .         -0.2661131  .          .       
#>  [3,]  .        -0.2808520  1.155093  .          .         -0.2808520  .       
#>  [4,] -0.280852  .          .         1.2291667 -0.2661131  .         -0.280852
#>  [5,]  .        -0.2661131  .        -0.2661131  1.2962963 -0.2661131  .       
#>  [6,]  .         .         -0.280852  .         -0.2661131  1.2291667  .       
#>  [7,]  .         .          .        -0.2808520  .          .          1.155093
#>  [8,]  .         .          .         .         -0.2661131  .         -0.280852
#>  [9,]  .         .          .         .          .         -0.2808520  .       
#>                           
#>  [1,]  .          .       
#>  [2,]  .          .       
#>  [3,]  .          .       
#>  [4,]  .          .       
#>  [5,] -0.2661131  .       
#>  [6,]  .         -0.280852
#>  [7,] -0.2808520  .       
#>  [8,]  1.2291667 -0.280852
#>  [9,] -0.2808520  1.155093
```

``` r
Q |> solve()
#> 9 x 9 sparse Matrix of class "dgCMatrix"
#>                                                                       
#>  [1,] 1.00000000 0.27611088 0.08016032 0.27611088 0.1353601 0.05357375
#>  [2,] 0.27611088 1.00000000 0.27611088 0.13559322 0.2783556 0.13559322
#>  [3,] 0.08016032 0.27611088 1.00000000 0.05357375 0.1353601 0.27611088
#>  [4,] 0.27611088 0.13559322 0.05357375 1.00000000 0.2783556 0.08474576
#>  [5,] 0.13536011 0.27835560 0.13536011 0.27835560 1.0000000 0.27835560
#>  [6,] 0.05357375 0.13559322 0.27611088 0.08474576 0.2783556 1.00000000
#>  [7,] 0.08016032 0.05357375 0.02605210 0.27611088 0.1353601 0.05357375
#>  [8,] 0.05357375 0.08474576 0.05357375 0.13559322 0.2783556 0.13559322
#>  [9,] 0.02605210 0.05357375 0.08016032 0.05357375 0.1353601 0.27611088
#>                                       
#>  [1,] 0.08016032 0.05357375 0.02605210
#>  [2,] 0.05357375 0.08474576 0.05357375
#>  [3,] 0.02605210 0.05357375 0.08016032
#>  [4,] 0.27611088 0.13559322 0.05357375
#>  [5,] 0.13536011 0.27835560 0.13536011
#>  [6,] 0.05357375 0.13559322 0.27611088
#>  [7,] 1.00000000 0.27611088 0.08016032
#>  [8,] 0.27611088 1.00000000 0.27611088
#>  [9,] 0.08016032 0.27611088 1.00000000
```

Creating and standardizing a 1600x1600 precision matrix

``` r
bench::mark(
  make_standardized_matern(dim = 40, rho = 0.5, nu = 0)
)
#> # A tibble: 1 × 6
#>   expression                             min median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr>                           <bch> <bch:>     <dbl> <bch:byt>    <dbl>
#> 1 make_standardized_matern(dim = 40, … 1.3ms 1.36ms      709.    98.3KB        0
```

# Sampling spatial data

Here we sample highly dependent spatial data on a 100x100 grid,
i.e. there’s 10.000 observational locations.

``` r
start <- tictoc::tic()
grid_dim <- 100
rho <- 0.9
nu <- 2
Z <- rmatern_copula(1, grid_dim, rho, nu)

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
#> 0.209 sec elapsed
```

# Normal density

The package also implements a method for calculating the log-density of
a multivariate normal with appropriate precision matrix. The function
avoids creating the precision matrix Q by using known results about
kroncker sums and eigendecompositions. This causes the density
evaluation to be blazingly fast, even for very large spatial fields.

The package also contains functions for sampling from and calculating
densities for regular Matérn-like fields where the marginal variances
don’t have to be equal to one. Those functions are even faster than the
copula versions.

``` r
library(purrr)
library(glue)
library(tinytable)
library(tidyr)
#> 
#> Attaching package: 'tidyr'
#> The following objects are masked from 'package:Matrix':
#> 
#>     expand, pack, unpack
```

``` r

my_fun <- function(dim) {
  x_copula <- rmatern_copula(1, dim, 0.5, 0)
  x <- rmatern(1, dim, 0.5, 0)
  bench::mark(
    "Copula" = dmatern_copula(x, dim, 0.5, 0),
    "Regular" = dmatern(x, dim, 0.5, 0),
    filter_gc = FALSE,
    iterations = 10,
    check = FALSE
  ) |> 
    mutate(
      dim = dim
    )
}

results <- map(seq(10, 200, by = 10), my_fun)


results |> 
  list_rbind() |> 
  select(Q_size = dim, type = expression, time = median, memory = mem_alloc) |> 
  mutate(
    Q_size = glue("{Q_size^2}x{Q_size^2}"),
    type = as.character(type),
    time = as.character(time)
  ) |> 
  select(-memory) |> 
  pivot_wider(names_from = type, values_from = time) |> 
  tt()
```

| Q_size      | Copula   | Regular  |
|-------------|----------|----------|
| 100x100     | 103.55µs | 94.71µs  |
| 400x400     | 226.94µs | 115.39µs |
| 900x900     | 701.98µs | 191.06µs |
| 1600x1600   | 1.84ms   | 407.64µs |
| 2500x2500   | 4.17ms   | 769.35µs |
| 3600x3600   | 8.45ms   | 1.48ms   |
| 4900x4900   | 17.1ms   | 2.52ms   |
| 6400x6400   | 29.95ms  | 3.91ms   |
| 8100x8100   | 50.48ms  | 14.15ms  |
| 10000x10000 | 77.45ms  | 8.47ms   |
| 12100x12100 | 99.32ms  | 11.99ms  |
| 14400x14400 | 138.19ms | 17.57ms  |
| 16900x16900 | 182.92ms | 24.35ms  |
| 19600x19600 | 237.26ms | 27.68ms  |
| 22500x22500 | 302.48ms | 40.63ms  |
| 25600x25600 | 407.55ms | 44.19ms  |
| 28900x28900 | 538.36ms | 58.09ms  |
| 32400x32400 | 654.52ms | 77.96ms  |
| 36100x36100 | 788.99ms | 89.31ms  |
| 40000x40000 | 1.03s    | 118.81ms |
