# Unit Variance Scaling This function defines a preprocessing strategy that is applied via [`prep_X`](https://tkimhofer.github.io/metabom8/reference/prep_X.md).

Unit Variance Scaling This function defines a preprocessing strategy
that is applied via
[`prep_X`](https://tkimhofer.github.io/metabom8/reference/prep_X.md).

## Usage

``` r
uv_scaling(center = TRUE)
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

Centers variables and scales each feature to unit variance. UV scaling
divides each variable by its standard deviation. This is the default
scaling in many multivariate methods such as PCA and PLS.

## See also

Other scaling_strategies:
[`pareto_scaling()`](https://tkimhofer.github.io/metabom8/reference/pareto_scaling.md),
[`unscaled()`](https://tkimhofer.github.io/metabom8/reference/unscaled.md)

## Examples

``` r
autoscale <- uv_scaling(center=TRUE)
X <- matrix(c(10,10, 0,0, 0, 10), ncol=3)
X_scaled <- prep_X(autoscale, X)
str(X_scaled)
#> List of 2
#>  $ X   : num [1:2, 1:3] 0 0 0 0 -0.707 ...
#>  $ prep:List of 4
#>   ..$ center: logi TRUE
#>   ..$ scale : chr "UV"
#>   ..$ X_mean: num [1:3] 10 0 5
#>   ..$ X_sd  : num [1:3] 0 0 7.07
X_scaled$X
#>      [,1] [,2]       [,3]
#> [1,]    0    0 -0.7071068
#> [2,]    0    0  0.7071068
```
