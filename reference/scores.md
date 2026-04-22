# PLS/OPLS model scores

PLS/OPLS model scores

## Usage

``` r
scores(object, ...)

# S4 method for class 'm8_model'
scores(object, orth = FALSE, cv = FALSE)
```

## Arguments

- object:

  An object of class `m8_model`.

- ...:

  Additional arguments (currently ignored).

- orth:

  Logical indicating whether orthogonal scores should be returned (only
  applicable for OPLS models).

- cv:

  Logical indicating whether cross-validated scores should be returned

## Value

Numeric vector or matrix containing scores.

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
scores(model, orth=FALSE)
#>            [,1]
#>  [1,]  64.35455
#>  [2,]  51.71614
#>  [3,]  79.38213
#>  [4,]  84.32423
#>  [5,]  81.29468
#>  [6,] -70.48437
#>  [7,] -69.99682
#>  [8,] -71.61648
#>  [9,] -71.89225
#> [10,] -77.08180
scores(model, orth=TRUE)
#>             [,1]
#>  [1,] -30.158164
#>  [2,]  -6.827912
#>  [3,] 197.529834
#>  [4,] -77.801381
#>  [5,] -82.742378
#>  [6,] -19.695453
#>  [7,] -14.447133
#>  [8,]  -5.123022
#>  [9,]  31.012710
#> [10,]   8.252897
scores(model, cv=TRUE)
#>            [,1]
#>  [1,]  15.00521
#>  [2,]  21.33020
#>  [3,]  85.42529
#>  [4,]  20.78544
#>  [5,]   0.00000
#>  [6,] -10.96025
#>  [7,] -22.67001
#>  [8,] -18.26880
#>  [9,] -16.62769
#> [10,] -22.07234
```
