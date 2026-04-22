# Pareto Scaling Leaves variables unscaled. Optional centering.

Pareto Scaling Leaves variables unscaled. Optional centering.

## Usage

``` r
pareto_scaling(center = FALSE)
```

## Arguments

- center:

  Logical. If TRUE, variables are mean-centered before scaling.

## Value

A `list` with elements:

- X:

  Numeric matrix containing the scaled data.

- prep:

  List describing the preprocessing, including centering and scaling
  parameters (`center`, `scale`, `X_mean`, `X_sd`).

## Details

Scales variables by the square root of their standard deviation.

## See also

Other scaling_strategies:
[`unscaled()`](https://tkimhofer.github.io/metabom8/reference/unscaled.md),
[`uv_scaling()`](https://tkimhofer.github.io/metabom8/reference/uv_scaling.md)

## Examples

``` r
paritalUV <- pareto_scaling(center=TRUE)
X <- matrix(c(10,10, 0,0, 0, 10, 0, 1000), ncol=4)
X_scaled <- prep_X(paritalUV, X)
str(X_scaled)
#> List of 2
#>  $ X   : num [1:2, 1:4] 0 0 0 0 -1.88 ...
#>  $ prep:List of 4
#>   ..$ center: logi TRUE
#>   ..$ scale : chr "Pareto"
#>   ..$ X_mean: num [1:4] 10 0 5 500
#>   ..$ X_sd  : num [1:4] 0 0 7.07 707.11
X_scaled$X
#>      [,1] [,2]      [,3]      [,4]
#> [1,]    0    0 -1.880302 -18.80302
#> [2,]    0    0  1.880302  18.80302
```
