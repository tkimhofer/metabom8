# Extract model weights

Extract model weights

## Usage

``` r
weights(object, ...)

# S4 method for class 'm8_model'
weights(object, orth = FALSE)
```

## Arguments

- object:

  An object of class `m8_model`.

- ...:

  Additional arguments (currently ignored).

- orth:

  Logical indicating whether orthogonal scores should be returned (only
  applicable for OPLS models).

## Value

Numeric vector or matrix containing model weights.

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
W <- weights(model)
Wo <- weights(model, orth = TRUE)
dim(W)
#> [1]     1 27819
dim(Wo) == dim(W)
#> [1] TRUE TRUE
```
