# No Scaling This function defines a preprocessing strategy that is applied via [`prep_X`](https://tkimhofer.github.io/metabom8/reference/prep_X.md).

No Scaling This function defines a preprocessing strategy that is
applied via
[`prep_X`](https://tkimhofer.github.io/metabom8/reference/prep_X.md).

## Usage

``` r
unscaled(center = FALSE)
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

\#' @details Leaves variables unscaled. Optional centering.

## See also

Other scaling_strategies:
[`pareto_scaling()`](https://tkimhofer.github.io/metabom8/reference/pareto_scaling.md),
[`uv_scaling()`](https://tkimhofer.github.io/metabom8/reference/uv_scaling.md)

## Examples

``` r
just_centering <- unscaled(center=TRUE)
X <- matrix(c(10,10, 0,0, 0, 10), ncol=3)
X_centered <- prep_X(just_centering, X)
str(X_centered)
#> List of 2
#>  $ X   : num [1:2, 1:3] 0 0 0 0 -5 5
#>  $ prep:List of 4
#>   ..$ center: logi TRUE
#>   ..$ scale : chr "None"
#>   ..$ X_mean: num [1:3] 10 0 5
#>   ..$ X_sd  : num [1:3] 0 0 7.07
X_centered$X
#>      [,1] [,2] [,3]
#> [1,]    0    0   -5
#> [2,]    0    0    5
```
