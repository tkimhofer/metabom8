# Extract fitted Y values

Extract fitted Y values

## Usage

``` r
fitted(object, ...)

# S4 method for class 'm8_model'
fitted(object)
```

## Arguments

- object:

  An object of class `m8_model`.

- ...:

  Additional arguments (currently ignored).

## Value

Numeric vector or matrix containing the fitted response values.

## Examples

``` r
data(covid)
cv <- balanced_mc(k=5, split=2/3)
scaling <- uv_scaling(center=TRUE)
model <-opls(X=covid$X, Y=covid$an$type, scaling, cv)
#> Performing discriminant analysis.
#> An O-PLS-DA model with 1 predictive and 1 orthogonal components was fitted.
show(model)
#> 
#> m8_model <opls>
#> ----------------------------------------
#> Dimensions  : 10 samples x 27819 variables
#> Mode        : classification 
#> Preprocess  : center | UV 
#> Components  : 2 (3 tested)
#> Validation  : BalancedMonteCarlo (k = 5) 
#> Stop rule   : cv_improvement_negligible 
#> ----------------------------------------
#> Use summary() for performance metrics.
#> 
Y_hat_dummy <- fitted(model)
```
