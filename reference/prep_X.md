# Applies a preprocessing strategy to a numeric matrix.

Applies a preprocessing strategy to a numeric matrix.

## Usage

``` r
prep_X(preproc_strategy, X)
```

## Arguments

- preproc_strategy:

  list with elements 'center' (logical) and 'scale' (character).

- X:

  numeric matrix - processing is done column-wise.

## Value

A list containing the processed matrix and parameters.

## Examples

``` r
X <- matrix(c(10,10, 0,0, 0, 10), ncol=3)
autoscale <- uv_scaling(center=TRUE)
X_scaled <- prep_X(autoscale, X)
str(X_scaled)
#> List of 2
#>  $ X   : num [1:2, 1:3] 0 0 0 0 -0.707 ...
#>  $ prep:List of 4
#>   ..$ center: logi TRUE
#>   ..$ scale : chr "UV"
#>   ..$ X_mean: num [1:3] 10 0 5
#>   ..$ X_sd  : num [1:3] 0 0 7.07
```
